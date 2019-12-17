###THIS MODULE WILL ACTUALLY CALL THE HAPLOTYPES
#!/bin/bash

cd TMP

working_directory=`cat RUN_DIR`
tmp_folder=`cat TMP_FOL`
complete_reference_database=`cat ALL_REF`
non_junction_markers=$working_directory/REF_SEQS/READ_RECOVERY/FOR_READ_RECOVERY_long_types_NON_JUNCTION_2019.fasta # THIS MAY NOT BE RELEVANT FOR OTHER PARASITES.

##MAKE LIST OF SPECIMENS TO GENOTYPE
cd $tmp_folder/PROCESSED_READS/
for file in *
do
echo "$file" >> $tmp_folder/SPECIMENS_TO_TYPE
done
cd $tmp_folder
/usr/local/bin/gsed -i 's/.clean_merged.fastq//g' SPECIMENS_TO_TYPE



## COMPILE REFERENCE DATABASES

cat $complete_reference_database > $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
cat $non_junction_markers >> $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta

cat $complete_reference_database > $tmp_folder/ALL_REFERENCES_FOR_BLASTING.fasta

file_check=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`


if [ -f $file_check ]; 
then 
cat $working_directory/CYCLO_REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
cat $working_directory/CYCLO_REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_folder/ALL_REFERENCES_FOR_BLASTING.fasta
fi




for SPECIMEN_NAME in `cat SPECIMENS_TO_TYPE`
do

#SPECIMEN_NAME=`cat SPECIMENS_TO_TYPE`


#### MAP READS TO THE REFERENCE DATABASE OF NON-JUNCTION MARKERS TO OBTAIN CORRECT READS
bwa index $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
bwa mem -t 10 $tmp_folder/READ_RECOVERY_REFERENCE_SEQUENCES.fasta $tmp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $tmp_folder/$SPECIMEN_NAME.mapped_only.fastq
#######################################################################################################################################




##### CLUSTER USING CD-HIT --- STEP 4
#FLAGS

# c = % similarity of bins
# g =  1 or 0, default 0 -- explanation below:
#By cd-hit’s default algorithm, a sequence is clustered to the first
#cluster that meet the threshold (fast mode). If set to 1, the program
#will cluster it into the most similar cluster that meet the threshold
#(accurate but slow mode)
# d = length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space


$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $tmp_folder/$SPECIMEN_NAME.mapped_only.fastq -o $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -c 1 -g 1 -d 0 -T 10


##STEP_5
####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
seqret -sequence $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -outseq $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta


###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT

perl $working_directory/make_multi_seq.pl $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta $tmp_folder/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq.clstr multi-seq 30  ###ADJUSTED JOEL CHECK THIS -- multi-seq determines how many reads is in a cluster before its written.
cd multi-seq
cat * >> $tmp_folder/$SPECIMEN_NAME.all_clusters.fasta
cd ..
rm -rf multi-seq
rm -rf clean_merged_CLUSTERS.*

#####STEP 7 IS TO TURN YOUR CLUSTERS INTO A BLAST DATABASE

$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $tmp_folder/$SPECIMEN_NAME.all_clusters.fasta \
-input_type fasta \
-dbtype nucl

#####STEP 8 IS TO TURN BLAST YOUR REFERENCES AGAINST THAT DATABASE OF CLUSTERS
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $SPECIMEN_NAME.all_clusters.fasta \
-query $tmp_folder/ALL_REFERENCES_FOR_BLASTING.fasta \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $working_directory/RESULTS/$SPECIMEN_NAME \
-max_target_seqs 1 \
-word_size 7 \
-outfmt "6 qseqid pident mismatch gapopen gaps sseqid sseq evalue bitscore"

########################################
rm $SPECIMEN_NAME.all_clusters.fasta ### 	IMPORTANT #######  THIS WILL STOP THE DOUBLING UP OF THE SPECIMENS
########################################

done
rm -rf *all_clusters.fast* *clean_merged*
rm SPECIMENS_TO_TYPE
cd $working_directory
