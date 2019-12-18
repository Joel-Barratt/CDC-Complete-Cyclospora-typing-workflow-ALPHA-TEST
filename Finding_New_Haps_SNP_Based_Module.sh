#!/bin/bash
#####################################################################################################################################################################
#### SEARCH FOR NEW HAPLOTYPES (NON-REPEAT) #########################################################################################################################
#####################################################################################################################################################################



#####################################################################################################################################################################
#### DEFINE VARIABLES - THESE ARE TXT FILES GENERATED DURING THE IMPORT AND TRIME PHASE  ############################################################################
#####################################################################################################################################################################

working_directory=`cat RUN_DIR`
non_junction_markers=$working_directory/REF_SEQS/READ_RECOVERY/FOR_READ_RECOVERY_long_types_NON_JUNCTION_2019.fasta
complete_reference_database=`cat ALL_REF`
tmp_directory=`cat TMP_FOL`
REQUIRED_DEPTH_NEW_HAP=`cat DEPTH_NEW`
#####################################################################################################################################################################


#


##############################################################################################################
#### GENERATE LIST OF TRIMMED READ FILES TO ANALYZE
##############################################################################################################
cd $tmp_directory/PROCESSED_READS/
for file in *
             do
                echo "$file" >> $tmp_directory/SPECIMENS_TO_SEARCH_OTHER
             done
cd $tmp_directory
/usr/local/bin/gsed -i 's/.clean_merged.fastq//g' SPECIMENS_TO_SEARCH_OTHER






for SPECIMEN_NAME in `cat SPECIMENS_TO_SEARCH_OTHER`


do


#### FILE CHECK 1 -- ARE THERE ANY NEW HAPLOTYPES THAT YOU MIGHT WANT TO ADD TO THE MAPPING REFERENCES?
##############################################################################################################
rm READ_RECOVERY_REFERENCE_SEQUENCES.fasta
cat $working_directory/REF_SEQS/READ_RECOVERY/FOR_READ_RECOVERY_long_types_NON_JUNCTION_2019.fasta > $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta

file_check_1=`echo "$tmp_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ `cat $file_check_1 |wc -l` -gt 0 ]; 
then 
cat $tmp_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
/usr/local/bin/gsed -i '/Junction/,+1 d' $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
fi
#############################################################################################################






#SPECIMEN_NAME=`cat SPECIMENS_TO_SEARCH_OTHER`

#### MAP READS TO THE REFERENCE DATABASE OF NON-JUNCTION MARKERS TO OBTAIN CORRECT READS
bwa index $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
bwa mem -t 10 READ_RECOVERY_REFERENCE_SEQUENCES.fasta $tmp_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $SPECIMEN_NAME.mapped_only.fastq
#######################################################################################################################################

#CD-HIT
#FLAGS

# c = % similarity of bins
# g =  1 or 0, default 0 -- explanation below:
#By cd-hit’s default algorithm, a sequence is clustered to the first
#cluster that meet the threshold (fast mode). If set to 1, the program
#will cluster it into the most similar cluster that meet the threshold
#(accurate but slow mode)
# d = length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space

$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $tmp_directory/$SPECIMEN_NAME.mapped_only.fastq -o $tmp_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -c 1 -g 1 -d 0 -T 10


####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
seqret -sequence $tmp_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -outseq $tmp_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta



###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT
perl $working_directory/make_multi_seq.pl $tmp_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta $tmp_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq.clstr multi-seq 30  ###ADJUSTED JOEL CHECK THIS



cd multi-seq
cat * >> $tmp_directory/$SPECIMEN_NAME.all_clusters.fasta
cd $tmp_directory



##########################################################################################################################################################################################################################
##MAKE A BLAST DATABASE FROM THE CLUSTERS --- THIS IS WHERE WE TRY TO IDENTIFY SEQUENCES FROM THE CLUSTERED RAW READS THAT MAY REPRESENT NEW TYPES
##############################################################################################################
# make a database of known non-junction markers to see if there is a match
rm ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
cat $complete_reference_database > $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta

file_check_2=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ `cat $file_check_2 |wc -l` -gt 0 ]; 
then 
cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
/usr/local/bin/gsed -i '/Junction/,+1 d' $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
fi
##############################################################################################################


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $tmp_directory/$SPECIMEN_NAME.all_clusters.fasta \
-input_type fasta \
-dbtype nucl


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $SPECIMEN_NAME.all_clusters.fasta \
-query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-word_size 7 \
-evalue 0.001 \
-perc_identity 80 \
-qcov_hsp_perc 100 \
-max_target_seqs 1000 \
-num_threads 10 \
-outfmt "6 qseqid sseq" \
-out $tmp_directory/$SPECIMEN_NAME.blast_result 
/usr/local/bin/gsed -i 's/^/>/' $tmp_directory/$SPECIMEN_NAME.blast_result
/usr/local/bin/gsed -i 's/\s/\n/g' $tmp_directory/$SPECIMEN_NAME.blast_result
awk -F'\n' 'BEGIN { ORS=""; RS=">"; N=-1 ; OFS="\n" } { if(N++ >= 0) { $1=">"N; print $0; } }' < $SPECIMEN_NAME.blast_result > $SPECIMEN_NAME.blast_result.fasta
##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $tmp_directory/$SPECIMEN_NAME.blast_result.fasta -o $tmp_directory/$SPECIMEN_NAME.blast_result_CLUSTERS.fasta -c 1 -g 1 -d 0 -T 10
##########################################################################################################################################################################################################################


# MAP READS TO THE OUTPUT OF PREVIOUS CLUSTERING WITH BOWTIE - USE HIGHLY STRINGENT PARAMETERS
#Parameters: -N = number of mismatches -L = length of seed substring -i = specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length.
# -D = Up to <int> consecutive seed extension attempts can “fail” before Bowtie 2 moves on. -R the maximum number of times Bowtie 2 will “re-seed” reads with repetitive seeds.
$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $tmp_directory/$SPECIMEN_NAME.blast_result_CLUSTERS.fasta $tmp_directory/$SPECIMEN_NAME.blast_result_CLUSTERS.fasta_BT_INDEX
$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $tmp_directory/$SPECIMEN_NAME.blast_result_CLUSTERS.fasta_BT_INDEX -U $tmp_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads 10 --local > $SPECIMEN_NAME.OUTPUT.bam



##### EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
samtools sort -T -n $tmp_directory/$SPECIMEN_NAME.OUTPUT.bam -o $tmp_directory/$SPECIMEN_NAME.aln.sorted.bam
samtools index $tmp_directory/$SPECIMEN_NAME.aln.sorted.bam
cat $tmp_directory/$SPECIMEN_NAME.blast_result_CLUSTERS.fasta | /usr/local/bin/gsed '/>/!d' > $tmp_directory/GOOD_HAPS_AMONG_CLUSTERS.txt
samtools depth -aa -d 0 $SPECIMEN_NAME.aln.sorted.bam > $tmp_directory/$SPECIMEN_NAME.coverage # exports a file containing the coverage.
##########################################################################################################################################################################################################################

### find out which haplotypes contain bases with coverage less than X reads
for NEW_HAPS in `cat GOOD_HAPS_AMONG_CLUSTERS.txt` 
      do
          if [ `cat $SPECIMEN_NAME.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $REQUIRED_DEPTH_NEW_HAP ];
               then echo $NEW_HAPS >> clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt
          fi
      done


#####  THIS IS WHERE YOU SET THE REQUIRED DEPTH ---- SEE ABOVE


##########################################################################################################################################################################################################################
## this will put underscores at the begining and end of fasta names -- if they are kept solely numerical, this can create problems
/usr/local/bin/gsed -i 's/>/>_/g' clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt
/usr/local/bin/gsed -i '/^>_/ s/$/_/' clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt
/usr/local/bin/gsed 's/>/>_/g' $SPECIMEN_NAME.blast_result_CLUSTERS.fasta | /usr/local/bin/gsed '/^>_/ s/$/_/' > $SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta
##########################################################################################################################################################################################################################

## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST OF NEW ONES THAT POSSESS LOW COVERAGE
for HAPS_TO_REMOVE in `cat clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt`
do
/usr/local/bin/gsed -i "/$HAPS_TO_REMOVE/,+1 d" $SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta #removes the line containing the fasta sequence name. The '+1' removes the line after as well.
done



#split each of the remaining sequences that meet the coverage cutoff into separate fasta files.
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
done < $tmp_directory/$SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta





for OTHER_HAPS in `ls -1`
do


##############################################################################################################
# make a database of known non-junction markers to see if there is a match
rm ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
cat $complete_reference_database > $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta

file_check_3=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ `cat $file_check_3 |wc -l` -gt 0 ]; 
then 
cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
/usr/local/bin/gsed -i '/Junction/,+1 d' $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
fi



#############################################################################################################

# NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in $tmp_directory/new_potential_otherhaps/$OTHER_HAPS -input_type fasta -dbtype nucl


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $OTHER_HAPS \
-query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-word_size 7 \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_directory/new_potential_otherhaps/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result \
-max_target_seqs 1 \
-outfmt "6 qseqid"


match_present=`cat $tmp_directory/new_potential_otherhaps/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result | wc -l`


if [ $match_present == 0 ];

then  $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $OTHER_HAPS \
-query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta -word_size 7 -evalue 0.001 \
-perc_identity 80 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_directory/new_potential_otherhaps/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY \
-max_target_seqs 1 \
-outfmt "6 qseqid"

which_marker=`/usr/local/bin/gsed 's/_Hap_.*/_Hap_/' $SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY | head -1`
hap_count=`cat $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta | grep $which_marker | wc -l`
hap_count_plus_one=`echo "$hap_count + 1" | bc`
cat $OTHER_HAPS | awk '/^>/{print ">'$which_marker$hap_count_plus_one'"; next}{print}' > $working_directory/REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.$which_marker$hap_count_plus_one.fasta
#cat $OTHER_HAPS | awk '/^>/{print ">'$which_marker$hap_count_plus_one'"; next}{print}' >> $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta


elif [ $match_present -gt 0 ];

then echo "not a new haplotype"

fi

done

cd $tmp_directory

rm -rf multi-seq new_potential_otherhaps *.sam *.bam *BT_INDEX* *_CLUSTERS* *.amb *.ann *.bwt *.pac *.sa *blast_result* *.coverage *.blast_result *coverage_less_than* *all_clusters* *mapped_only* *.REMOVE_BAD_SEQUENCES.fasta SPECIMENS_TO_SEARCH_OTHER *.aln.sorted.*

done















