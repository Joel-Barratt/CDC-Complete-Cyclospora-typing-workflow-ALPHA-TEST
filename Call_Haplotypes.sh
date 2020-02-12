###THIS MODULE WILL ACTUALLY CALL THE HAPLOTYPES
#!/bin/bash

#cd TMP #see what happens if you unhash this. I just hashed it.

working_directory=`cat RUN_DIR`
tmp_folder=`cat TMP_FOL`
complete_reference_database=`cat ALL_REF`
#non_junction_markers=$working_directory/REF_SEQS/READ_RECOVERY/FOR_READ_RECOVERY_long_types_NON_JUNCTION_2019.fasta # THIS MAY NOT BE RELEVANT FOR OTHER PARASITES.

minimum_depth_for_haplotype_assignment=`cat DEPTH_FOR_CALLING`


##MAKE LIST OF SPECIMENS TO GENOTYPE
cd $tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/
for file in *
do
echo "$file" >> $tmp_folder/SPECIMENS_TO_TYPE
done
cd $tmp_folder
/usr/local/bin/gsed -i 's/.FINAL_CLIPPED_FOR_HAP_CALLING.fastq//g' SPECIMENS_TO_TYPE


## COMPILE REFERENCE DATABASES

rm ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta READ_RECOVERY_REFERENCE_SEQUENCES.fasta
cat $complete_reference_database > $tmp_folder/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
file_check=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ `cat $file_check |wc -l` -gt 0 ]; 
then 
cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_folder/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
fi


#####STEP 7 IS TO TURN YOUR CLUSTERS INTO A BLAST DATABASE

$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $tmp_folder/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-input_type fasta \
-dbtype nucl


for SPECIMEN_NAME in `cat SPECIMENS_TO_TYPE`
do

#SPECIMEN_NAME=`cat SPECIMENS_TO_TYPE | head -1`


#### MAP READS TO THE REFERENCE DATABASE OF NON-JUNCTION MARKERS TO OBTAIN CORRECT READS
#bwa index $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
#bwa mem -t 10 $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta $tmp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
#samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
#samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $tmp_folder/$SPECIMEN_NAME.mapped_only.fastq
#######################################################################################################################################



##RENAME READS THAT ARE ALREADY TRIMMED FOR HAP CALLING

cat $tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAP_CALLING.fastq | awk '{if(NR%4==1) $0=sprintf("@Read_%d",(1+i++)); print;}' >  \
$tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/RENAMED.$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAP_CALLING.fastq


##### CLUSTER USING CD-HIT --- STEP 4
#FLAGS

# c = % similarity of bins
# g =  1 or 0, default 0 -- explanation below:
#By cd-hit’s default algorithm, a sequence is clustered to the first
#cluster that meet the threshold (fast mode). If set to 1, the program
#will cluster it into the most similar cluster that meet the threshold
#(accurate but slow mode)
# d = length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space
#-s   length difference cutoff, default 0.0     --- IT IS ESSENTIAL THAT YOU CHANGE THIS TO 1!!!!!!

$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/RENAMED.$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAP_CALLING.fastq \
-o $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -c 1 -g 1 -d 0 -T 10 -s 1


##STEP_5
####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
seqret -sequence $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -outseq $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta


###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT

perl $working_directory/make_multi_seq.pl $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq.clstr multi-seq $minimum_depth_for_haplotype_assignment  ###ADJUSTED JOEL CHECK THIS -- multi-seq determines how many reads is in a cluster before its written.
#cd multi-seq
#cat * >> $tmp_folder/$SPECIMEN_NAME.all_clusters.fasta
#cd ..
cp -r $tmp_folder/multi-seq $tmp_folder/$SPECIMEN_NAME.hap_call.multi-seq
rm -rf $tmp_folder/multi-seq clean_merged_CLUSTERS.*


cd $SPECIMEN_NAME.hap_call.multi-seq


####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA
seqret -sequence $tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/RENAMED.$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAP_CALLING.fastq \
-outseq $tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/RENAMED.$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAPLOTYPE_CALLING.fasta

mkdir $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP



for cluster_file in *
do


#cluster_file=26


#####STEP 8 IS TO TURN BLAST YOUR REFERENCES AGAINST THAT DATABASE OF CLUSTERS
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_folder/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-query $cluster_file \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/temp_result \
-max_target_seqs 1 \
-word_size 7 \
-dust no \
-soft_masking false \
-outfmt "6 sseqid"

test_for_a_decent_hit=`cat $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/temp_result | grep "_Hap_" | wc -l`

if [  $test_for_a_decent_hit -gt 0 ];

then

tmp_var=`cat $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/temp_result | /usr/local/bin/gsed 's/_Hap_.*//'`

cat $cluster_file | /usr/local/bin/gsed "1s/.*/>$tmp_var/" > $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$tmp_var.fasta
rm $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/temp_result

dynamic_coverage_cutoff=`cat $tmp_folder/$SPECIMEN_NAME.COVERAGE_CUTOFF_PER_LOCUS_FOR_HAP_CALLING/$SPECIMEN_NAME.$tmp_var.COVERAGE_CUTOFF`


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$tmp_var.fasta \
-input_type fasta \
-dbtype nucl

$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$tmp_var.fasta \
-query $tmp_folder/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/RENAMED.$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAPLOTYPE_CALLING.fasta \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$SPECIMEN_NAME.$tmp_var.raw_fasta_reads_with_hits \
-max_target_seqs 1 \
-word_size 7 \
-dust no \
-soft_masking false \
-outfmt "6 qseqid"

coverage_estimate_by_blast_hap_calling=`cat $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$SPECIMEN_NAME.$tmp_var.raw_fasta_reads_with_hits | grep "Read_" | wc -l`

if [ $coverage_estimate_by_blast_hap_calling -gt $dynamic_coverage_cutoff ] && [ $coverage_estimate_by_blast_hap_calling -gt $minimum_depth_for_haplotype_assignment ];

then cat $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$tmp_var.fasta >> $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$SPECIMEN_NAME.all_clusters.fasta

fi

elif [  $test_for_a_decent_hit == 0 ];

then

rm $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/temp_result

fi

done


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$SPECIMEN_NAME.all_clusters.fasta \
-input_type fasta \
-dbtype nucl


#####STEP 8 IS TO TURN BLAST YOUR REFERENCES AGAINST THAT DATABASE OF CLUSTERS
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP/$SPECIMEN_NAME.all_clusters.fasta \
-query $tmp_folder/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $working_directory/RESULTS/$SPECIMEN_NAME \
-max_target_seqs 1 \
-word_size 7 \
-dust no \
-soft_masking false \
-outfmt "6 qseqid pident mismatch gapopen gaps sseq evalue bitscore"

cd $tmp_folder

###############################################
#rm -r $tmp_folder/$SPECIMEN_NAME.BLASTING_TEMP#
###############################################

done
rm -rf *all_clusters.fast* *clean_merged*




###THIS IS ONLY FOR THE JUNCTION -- hash this part out if you dont want to run the junction.

#cd CONFIRMED_JUNCTIONS


for SPECIMEN_NAME in `cat SPECIMENS_TO_TYPE`
do


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $tmp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.validated_junction_sequences.fasta \
-input_type fasta \
-dbtype nucl


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.validated_junction_sequences.fasta \
-query $tmp_folder/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.junction_only \
-max_target_seqs 1 \
-word_size 7 \
-dust no \
-soft_masking false \
-outfmt "6 qseqid pident mismatch gapopen gaps sseq evalue bitscore"


cat $tmp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.junction_only >> $working_directory/RESULTS/$SPECIMEN_NAME

#cat $tmp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.junction_only >> $working_directory/RESULTS/$SPECIMEN_NAME.junction.txt


done

cd $tmp_folder



#rm SPECIMENS_TO_TYPE
cd $working_directory
