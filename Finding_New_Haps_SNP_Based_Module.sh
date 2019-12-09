#!/bin/bash

rm SPECIMENS_TO_SEARCH_OTHER

/usr/local/bin/gsed -i 's@.*/@@' NOT_JUNC_REF

working_directory=`cat RUN_DIR`

non_junction_markers=`cat NOT_JUNC_REF`

complete_reference_database=`cat ALL_REF`



cd $working_directory/PROCESSED_READS/


##MAKE LIST OF SPECIMENS TO SEARCH FOR NEW JUNCTION SEQUENCES
for file in *
do
echo "$file" >> ../SPECIMENS_TO_SEARCH_OTHER
done
cd $working_directory

/usr/local/bin/gsed -i 's/.clean_merged.fastq//g' SPECIMENS_TO_SEARCH_OTHER



for SPECIMEN_NAME in `cat SPECIMENS_TO_SEARCH_OTHER`


do


#SPECIMEN_NAME=`cat SPECIMENS_TO_SEARCH_OTHER`


bwa index $working_directory/$non_junction_markers
bwa mem $non_junction_markers $working_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $SPECIMEN_NAME.mapped_only.fastq




##### CLUSTER USING CD-HIT --- STEP 4
#FLAGS

# c = % similarity of bins
# g =  1 or 0, default 0 -- explanation below:
#By cd-hit’s default algorithm, a sequence is clustered to the first
#cluster that meet the threshold (fast mode). If set to 1, the program
#will cluster it into the most similar cluster that meet the threshold
#(accurate but slow mode)
# d = length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space


$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $working_directory/$SPECIMEN_NAME.mapped_only.fastq -o $working_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -c 1 -g 1 -d 0 -T 10

##STEP_5
####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
seqret -sequence $working_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -outseq $working_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta

###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT
cd $working_directory
perl make_multi_seq.pl $SPECIMEN_NAME.clean_merged_CLUSTERS.fasta $SPECIMEN_NAME.clean_merged_CLUSTERS.fq.clstr multi-seq 10  ###ADJUSTED JOEL CHECK THIS



cd multi-seq
cat * >> ../$SPECIMEN_NAME.all_clusters.fasta
cd ..
#rm -rf multi-seq
#rm -rf clean_merged_CLUSTERS.*




#NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.

$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER/$SPECIMEN_NAME.all_clusters.fasta -input_type fasta -dbtype nucl -title junction_check
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $SPECIMEN_NAME.all_clusters.fasta -query $complete_reference_database -word_size 7 -evalue 0.001 -perc_identity 80 -qcov_hsp_perc 100 -out $working_directory/$SPECIMEN_NAME.blast_result -max_target_seqs 30 -outfmt "6 qseqid qseq"

/usr/local/bin/gsed -i 's/^/>/' $SPECIMEN_NAME.blast_result
/usr/local/bin/gsed -i 's/\s/\n/g' $SPECIMEN_NAME.blast_result

awk -F'\n' 'BEGIN { ORS=""; RS=">"; N=-1 ; OFS="\n" } { if(N++ >= 0) { $1=">"N; print $0; } }' < $SPECIMEN_NAME.blast_result > $SPECIMEN_NAME.blast_result.fasta

$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $working_directory/$SPECIMEN_NAME.blast_result.fasta -o $working_directory/$SPECIMEN_NAME.blast_result_CLUSTERS.fasta -c 1 -g 1 -d 0 -T 10




#Parameters: -N = number of mismatches -L = length of seed substring -i = specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length.
# -D = Up to <int> consecutive seed extension attempts can “fail” before Bowtie 2 moves on. -R the maximum number of times Bowtie 2 will “re-seed” reads with repetitive seeds.


$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $SPECIMEN_NAME.blast_result_CLUSTERS.fasta $SPECIMEN_NAME.blast_result_CLUSTERS.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $SPECIMEN_NAME.blast_result_CLUSTERS.fasta_BT_INDEX -U $working_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads 10 --local > OUTPUT.bam


samtools sort -T -n OUTPUT.bam -o aln.sorted.bam

samtools index aln.sorted.bam


cat $SPECIMEN_NAME.blast_result_CLUSTERS.fasta | /usr/local/bin/gsed '/>/!d' > LIST_TO_RUN.txt

samtools depth -aa -d 0 aln.sorted.bam > $SPECIMEN_NAME.coverage # exports a file containing the coverage.


for NEW_HAPS in `cat LIST_TO_RUN.txt` 

do

#NEW_HAPS=`cat LIST_TO_RUN.txt`

if [ `cat $SPECIMEN_NAME.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt 20 ];

then echo $NEW_HAPS >> haps.to.remove.txt

fi

done


/usr/local/bin/gsed 's/>/>_/g' $SPECIMEN_NAME.blast_result_CLUSTERS.fasta | /usr/local/bin/gsed '/^>_/ s/$/_/' > $SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta


for HAPS_TO_REMOVE in `cat haps.to.remove.txt`

do

/usr/local/bin/gsed -i "/$HAPS_TO_REMOVE/,+1 d" $SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta

done




#USE THIS TO SPLIT ALL THE CLUSTERS
mkdir new_potential_otherhaps

cd new_potential_otherhaps

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fasta
        echo $line > $outfile
    else
        echo $line >> $outfile
   fi
done < ../$SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta





#### 
##############################################################################################################360i2 A
cat $complete_reference_database > $working_directory/ALL_REFERENCE_SEQUENCES.fasta

file_check=`echo "$working_directory/CYCLO_REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ -f $file_check ]; 
then 
cat $working_directory/CYCLO_REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $working_directory/ALL_REFERENCE_SEQUENCES.fasta

else echo "NO NEW HAPLOTYPES FOUND PREVIOUSLY -- LETS CHECK THE LATEST SPECIMEN"
fi

#############################################################################################################


complete_reference_with_new_haplotypes_included



for OTHER_HAPS in `ls -1`
do

complete_reference_with_new_haplotypes_included=$working_directory/ALL_REFERENCE_SEQUENCES.fasta

#OTHER_HAPS=`ls -1 | head -1`


##############################################################################################################
#SPECIMEN_NAME=`echo "SNY12411-19"`

# NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in $working_directory/new_potential_otherhaps/$OTHER_HAPS -input_type fasta -dbtype nucl


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $OTHER_HAPS -query $complete_reference_with_new_haplotypes_included -word_size 7 \
-evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -out $working_directory/new_potential_otherhaps/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result \
-max_target_seqs 1 -outfmt "6 qseqid"

match_present=`cat $working_directory/new_potential_otherhaps/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result | wc -l`



if [ $match_present == 0 ];

then  $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $SPECIMEN_NAME.TEST_$OTHER_HAPS -query $complete_reference_with_new_haplotypes_included -word_size 7 -evalue 0.001 \
-perc_identity 80 -qcov_hsp_perc 100 -out $working_directory/new_potential_otherhaps/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY -max_target_seqs 1 -outfmt "6 qseqid"

which_marker=`/usr/local/bin/gsed 's/.$//' $SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY`
hap_count=`grep $which_marker $complete_reference_with_new_haplotypes_included | wc -l` ############################################################################################################################## YOU DONT WANT TO READ THIS FROM THE ORIGINAL REFERENCE DB BECAUSE THIS DOES NOT CONSIDER NEWLY DISCOVERED HAPLOTYPES. THIS IS IMPORTANT TO CORRECT.
hap_count_plus_one=`echo "$hap_count + 1" | bc`
cat $OTHER_HAPS | awk '/^>/{print ">'$which_marker$hap_count_plus_one'"; next}{print}' > $working_directory/CYCLO_REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.$which_marker$hap_count_plus_one.fasta


elif [ $match_present -gt 0 ];

then echo "not a new haplotype"

fi

done


cd ..

#rm -rf
















