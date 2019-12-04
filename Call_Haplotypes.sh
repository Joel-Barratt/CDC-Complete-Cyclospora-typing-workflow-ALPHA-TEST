###THIS MODULE WILL ACTUALLY CALL THE HAPLOTYPES

rm SPECIMENS_TO_TYPE

running_directory=`cat RUN_DIR`

complete_reference_database=`cat ALL_REF`


##MAKE LIST OF SPECIMENS TO GENOTYPE

cd $running_directory/PROCESSED_READS/

for file in *
do
echo "$file" >> ../SPECIMENS_TO_TYPE
done
cd $running_directory

/usr/local/bin/gsed -i 's/.clean_merged.fastq//g' SPECIMENS_TO_TYPE

mkdir TYPING_RESULTS


for SPECIMEN_NAME in `cat SPECIMENS_TO_TYPE`
do

#SPECIMEN_NAME=`cat SPECIMENS_TO_TYPE`

##### CLUSTER USING CD-HIT --- STEP 4
#FLAGS

# c = % similarity of bins
# g =  1 or 0, default 0 -- explanation below:
#By cd-hit’s default algorithm, a sequence is clustered to the first
#cluster that meet the threshold (fast mode). If set to 1, the program
#will cluster it into the most similar cluster that meet the threshold
#(accurate but slow mode)
# d = length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space


$running_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $running_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -o $running_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -c 1 -g 1 -d 0 -T 10




##STEP_5
####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
seqret -sequence $running_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fq -outseq $running_directory/$SPECIMEN_NAME.clean_merged_CLUSTERS.fasta




###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT
cd $running_directory
perl make_multi_seq.pl $SPECIMEN_NAME.clean_merged_CLUSTERS.fasta $SPECIMEN_NAME.clean_merged_CLUSTERS.fq.clstr multi-seq 10  ###ADJUSTED JOEL CHECK THIS
cd multi-seq
cat * >> ../$SPECIMEN_NAME.all_clusters.fasta
cd ..
rm -rf multi-seq
rm -rf clean_merged_CLUSTERS.*


#####STEP 7 IS TO TURN YOUR CLUSTERS INTO A BLAST DATABASE

$running_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in $running_directory/$SPECIMEN_NAME.all_clusters.fasta -input_type fasta -dbtype nucl


#####STEP 8 IS TO TURN BLAST YOUR REFERENCES AGAINST THAT DATABASE OF CLUSTERS
$running_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $SPECIMEN_NAME.all_clusters.fasta -query $complete_reference_database \
-evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -out $running_directory/TYPING_RESULTS/$SPECIMEN_NAME \
-max_target_seqs 1 -word_size 7 -outfmt "6 qseqid pident mismatch gapopen gaps sseqid sseq evalue bitscore"


########################################
rm $SPECIMEN_NAME.all_clusters.fasta ### 	IMPORTANT #######  THIS WILL STOP THE DOUBLING UP OF THE SPECIMENS
########################################


done
rm -rf *all_clusters.fast* *clean_merged*
rm SPECIMENS_TO_TYPE

