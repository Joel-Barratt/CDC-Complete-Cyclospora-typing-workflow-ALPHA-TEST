#!/bin/bash

rm SPECIMENS_TO_SEARCH

/usr/local/bin/gsed -i 's@.*/@@' JUNC_REF

running_directory=`cat RUN_DIR`

junction_with_primer=`cat JUNC_REF`

complete_reference_database=`cat ALL_REF`



cd $running_directory/PROCESSED_READS/


##MAKE LIST OF SPECIMENS TO GENOTYPE
for file in *
do
echo "$file" >> ../SPECIMENS_TO_SEARCH
done
cd $running_directory/

/usr/local/bin/gsed -i 's/.clean_merged.fastq//g' SPECIMENS_TO_SEARCH



for SPECIMEN_NAME in `cat SPECIMENS_TO_SEARCH`


do


#SPECIMEN_NAME=`cat SPECIMENS_TO_SEARCH`


bwa index $running_directory/$junction_with_primer
bwa mem $junction_with_primer $running_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $SPECIMEN_NAME.mapped_only.fastq



# convert the fastq sequence just created to a fasta
$running_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/miraconvert $SPECIMEN_NAME.mapped_only.fastq $SPECIMEN_NAME.mapped_only.fasta



#make those fasta sequences a single line - this is necessary for the string searches that will follow
cat $SPECIMEN_NAME.mapped_only.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' >> $SPECIMEN_NAME.mapped_only_single_line.fasta
awk 'NF > 0' $SPECIMEN_NAME.mapped_only_single_line.fasta > $SPECIMEN_NAME.mapped_only_single_line2.fasta

###now we collect the reads that were mapped to the junction sequences and then split them into their own fasta files in the directory "mapped_only"


mkdir mapped_only
cd mapped_only



while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fasta
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < ../$SPECIMEN_NAME.mapped_only_single_line2.fasta


##TIDY UP
rm -rf ../$SPECIMEN_NAME.*.sam
rm -rf ../$SPECIMEN_NAME.mapped_only.*



## now we need to find reads that have the primer between each end (or partial primer). Then we need to list the names of those reads.
#>junction_primer_forward
#TACCAAAGCATCCATCTACAGC
#>junction_primer_reverse-reverse_complement
#AACACGATCCGATTGCTTGGG
#>reverse_sequence_of_forward_junction_primer
#GCTGTAGATGGATGCTTTGGTA
#>junction_reverse_primer
#CCCAAGCAATCGGATCGTGTT


####find those sequences that mapped to the junction that contain a section of both primers

for my_fastas in *
do
grep -l 'TACCAAAGCATCCATCTACAGC.*AACACGATCCGATTGCTTGGG' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
grep -l 'CCCAAGCAATCGGATCGTGTT.*GCTGTAGATGGATGCTTTGGTA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
done



mkdir ../assemble_these

for READS in `cat ../$SPECIMEN_NAME.list_of_sequences_to_review.txt`
do
cat $READS >> $SPECIMEN_NAME.assemble_these_seqs
done
cp $SPECIMEN_NAME.assemble_these_seqs ../assemble_these
cd ../assemble_these

$running_directory/CAP3/cap3.macosx.intel64/cap3 $SPECIMEN_NAME.assemble_these_seqs -o 90 -p 95


length1=`cat *.contigs | wc -l`
if [ $length1 -gt 1 ];
then cat *.contigs >> ../$SPECIMEN_NAME.map_to_this.fasta

fi

length2=`cat *.singlets | wc -l`
if [ $length2 -gt 1 ];
then cat *.singlets >> ../$SPECIMEN_NAME.map_to_this.fasta

fi

cd ..


###TIDY UP

rm -rf mapped_only *mapped_only_single_line.fasta assemble_these *list_of_sequences_to_review


cat $SPECIMEN_NAME.map_to_this.fasta | awk '/^>/{print ">new_junction_" ++i "_sequence"; next}{print}' > $SPECIMEN_NAME.map_to_this_2.fasta

seqtk seq -F '#' $SPECIMEN_NAME.map_to_this_2.fasta > $SPECIMEN_NAME.map_to_this.fastq

cp $SPECIMEN_NAME.map_to_this.fastq map_to_this.fastq
cp ./PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq clean_merged.fastq

$running_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/mira Junction_manifest_file.txt

cp /$running_directory/novel_junction_finder_assembly/novel_junction_finder_d_results/novel_junction_finder_out_new_junc.unpadded.fas* .
$running_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/miraconvert novel_junction_finder_out_new_junc.unpadded.fasta $SPECIMEN_NAME.newly_found_junction.fastq

## TIDY UP
rm -rf novel_junction_finder_assembly novel_junction_finder_out* *map_to_this.fastq



#### CUTADAPT TO TRIM OFF THE PRIMERS
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.newly_found_junction.fastq -g TACCAAAGCATCCATCTACAGC --output $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq ### five prime adapter only
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACACGATCCGATTGCTTGGG --output $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq ### THREE prime adapter only
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g CCCAAGCAATCGGATCGTGTT --output $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq ### five prime adapter only   ----- IF CONTIG IS REVERSED
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGTAGATGGATGCTTTGGTA --output $SPECIMEN_NAME.REV_three_prime_trim.fasta  --fasta ### five prime adapter only   ----- IF CONTIG IS REVERSED
cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $SPECIMEN_NAME.novel_junction_sequence_found.fasta




rm -rf $SPECIMEN_NAME.CUTADAPT_* $SPECIMEN_NAME.newly_found_junction.fastq $SPECIMEN_NAME.REV_three_prime_trim.fasta clean_merged*



#This will split the newly discovered junction haplotypes				IS THERE A BETTER WAY TO DO THIS? SOMETIMES I GET AN ERROR AND IT IS NOT CLEAR WHY
mkdir new_potential_junctions

cd new_potential_junctions

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fasta
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < ../$SPECIMEN_NAME.novel_junction_sequence_found.fasta



for FOOBAR in `ls -1`
do 

junction_count=`grep 'Junction_Hap' $complete_reference_database | wc -l`
new_number=`echo "$junction_count + 1" | bc`

cat $FOOBAR | awk '/^>/{print ">Junction_Hap_'$new_number'"; next}{print}' > $SPECIMEN_NAME.TEST_$FOOBAR

# NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.
$running_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER/$junction_with_primer -input_type fasta -dbtype nucl -title junction_check
$running_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db ../$junction_with_primer -query $SPECIMEN_NAME.TEST_$FOOBAR -word_size 7 -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -out $running_directory/$SPECIMEN_NAME.RESULT_$FOOBAR.blast_result -max_target_seqs 1 -outfmt "6 qseqid pident mismatch gapopen gaps sseqid sseq evalue bitscore"

match_present=`cat $running_directory/$SPECIMEN_NAME.RESULT_$FOOBAR.blast_result | wc -l`

if [ $match_present == 0 ];
then echo "you have discovered a new junction type"
cat $SPECIMEN_NAME.TEST_$FOOBAR >> $complete_reference_database
cat $SPECIMEN_NAME.TEST_$FOOBAR >> $running_directory/$junction_with_primer

elif [ $match_present -gt 0 ];
then echo "no new junction sequences were found"
fi

#mv $SPECIMEN_NAME* TRASH

done
cd $running_directory
rm -rf new_potential_junctions *fasta.amb *fasta.ann *fasta.bwt *fasta.nhr *fasta.nin *fasta.nsq *fasta.pac *fasta.sa *TEST* *review.txt *RESULT* *novel* *mapped_only* *map_to_this*
done

rm SPECIMENS_TO_SEARCH



