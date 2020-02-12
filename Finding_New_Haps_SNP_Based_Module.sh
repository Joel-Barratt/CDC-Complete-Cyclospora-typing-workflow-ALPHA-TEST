#!/bin/bash
#####################################################################################################################################################################
#### SEARCH FOR NEW HAPLOTYPES (NON-REPEAT) #########################################################################################################################
#####################################################################################################################################################################

### need pysam installed. pip install pysam

#####################################################################################################################################################################
#### DEFINE VARIABLES - THESE ARE TXT FILES GENERATED DURING THE IMPORT AND TRIME PHASE  ############################################################################
#####################################################################################################################################################################

working_directory=`cat RUN_DIR`
read_recovery_refs=`cat READ_REC` ######                  THIS LINE MIGHT NEED TO BE CHANGED FOR CYCLOSPORA.
complete_reference_database=`cat ALL_REF`
tmp_directory=`cat TMP_FOL`
REQUIRED_DEPTH_NEW_HAP=`cat DEPTH_NEW`
bed_references=`cat BED_LOCATION`
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






cat $bed_references | cut -f4 > $tmp_directory/LIST_OF_LOCI_UNIQUE
mkdir $tmp_directory/CLIPPED_READS_FOR_HAPLOTYPE_CALLING


for SPECIMEN_NAME in `cat SPECIMENS_TO_SEARCH_OTHER`


do

## hash out the line below.
#SPECIMEN_NAME=`cat SPECIMENS_TO_SEARCH_OTHER | head -1`

mkdir $tmp_directory/$SPECIMEN_NAME.COVERAGE_CUTOFF_PER_LOCUS_FOR_HAP_CALLING


#### FILE CHECK 1 -- ARE THERE ANY NEW HAPLOTYPES THAT YOU MIGHT WANT TO ADD TO THE MAPPING REFERENCES?
##############################################################################################################
rm READ_RECOVERY_REFERENCE_SEQUENCES.fasta
cat $read_recovery_refs > $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta

file_check_1=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ `cat $file_check_1 |wc -l` -gt 0 ]; 
then 
cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
/usr/local/bin/gsed -i '/Junction/,+1 d' $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta # THIS LINE IS HASHED OUT BECAUSE IT IS ONLY RELEVANT TO CYCLOSPORA.
fi
#############################################################################################################

#### MAP READS TO THE REFERENCE DATABASE OF NON-JUNCTION MARKERS TO OBTAIN CORRECT READS
bwa index $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
bwa mem -t 10 READ_RECOVERY_REFERENCE_SEQUENCES.fasta $tmp_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' | awk '{if(NR%4==1) $0=sprintf("@Read_%d",(1+i++)); print;}' > $SPECIMEN_NAME.mapped_only.fastq  # the last awk bit renames the reads 
#######################################################################################################################################


#### KEEP HERE -- MAY NEED TO REMOVE JUNCTION FROM THIS PART OF THE REFERENCE SEQEUENCE FOR CYCLOSPORA.

#### NOW WE NEED TO MAKE A CONSENSUS OF HAPLOTYPES TO MAP TO


#cat $complete_reference_database > $tmp_directory/MAPPING_REFS_TO_MAKE_CONSENSUS.fasta

#file_check_1=`echo "$tmp_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

#if [ `cat $file_check_1 |wc -l` -gt 0 ]; 
#then 
#cat $tmp_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/MAPPING_REFS_TO_MAKE_CONSENSUS.fasta
#/usr/local/bin/gsed -i '/Junction/,+1 d' $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta # THIS LINE IS HASHED OUT BECAUSE IT IS ONLY RELEVANT TO CYCLOSPORA. THIS LINE REMOVES THE JUNCTION FROM THIS STEP WHICH IS IMPORTANT BECAUSE JUNCTION HAP CALLING IS A NEW WORKFLOW.
#fi
#############################################################################################################



##### now we align to the sequences in order to calculate the depth of sequencing obtained for each marker, determine the cutoff and trim the reads


cd $tmp_directory

cat $read_recovery_refs > $tmp_directory/FULL_LENGTH_MARKERS.fasta

mkdir $SPECIMEN_NAME.MAPPING
cd $SPECIMEN_NAME.MAPPING


###### NOW WE MAP READS TO THE CONSENSUS SEQUENCE OF EVERY MARKER --- MAP TO EACH OF THE LONGER MARKERS INSTEAD - WE WILL USE THE BED FILE TO EXTRACT CORRECT REGION.

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $tmp_directory/FULL_LENGTH_MARKERS.fasta $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.BT_INDEX
$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.BT_INDEX \
-U $tmp_directory/$SPECIMEN_NAME.mapped_only.fastq --threads 10 --local > $SPECIMEN_NAME.alignment.bam




samtools sort $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.alignment.bam > $tmp_directory/$SPECIMEN_NAME.MAPPING/sorted.$SPECIMEN_NAME.alignment.bam
samtools index sorted.$SPECIMEN_NAME.alignment.bam


### THE LINE BELOW - YOU NEED TO USE THE FOR cat LIST_OF_LOCI_UNIQUE   -- PULL THESE CURRENT LINE FROM THE LARGE BED REFERENCE FILE AND MAKE A SMALLER BED FILE
## THEN EXTRACT READS OF THE CORRECT LOCATION FROM THIS.

##THIS LOOP IS CURRENTLY HAVING PROBLEMS.
for find_depth_and_alignments in `cat $tmp_directory/LIST_OF_LOCI_UNIQUE`
 
                do

#hash out line below
#find_depth_and_alignments=`cat $tmp_directory/LIST_OF_LOCI_UNIQUE | head -1`
#find_depth_and_alignments=Mt_MSR_PART_A


                         #line below is a problem -- find_depth_and_alignments should be column 4 of the bed file which is the contents of the LIST_OF_LOCI_UNIQUE text file.
                         marker_range=`echo $find_depth_and_alignments | /usr/local/bin/gsed 's/_PART_*.//g'`
                                          
                         
                         cat $bed_references | grep $find_depth_and_alignments > $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed
                         left_trim=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed | cut -f2`
                         right_trim=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed | cut -f3`
                         

                         ###THIS LINE IS 100%  NECESSARY TO GET THE TRIMMING DONE PROPERLY
                         #IF A READ DOES NOT COMPLETELY SPAN THE REGION, IT WILL NOT BE TRIMMED PROPERLY.



####THE PROBLEM IS SAMTOOLS VIEW HERE -- IT IS NOT RECOGNISING THE VARIABLE
                         samtools view  -bu -F 4  sorted.$SPECIMEN_NAME.alignment.bam $marker_range | \
                         java -jar /Users/joelbarratt/jvarkit/dist/samjs.jar  \
                         -e "record.alignmentStart <= $left_trim && record.alignmentEnd >= $right_trim"  --samoutputformat BAM \
                         -o $tmp_directory/$SPECIMEN_NAME.MAPPING/OVERLAP_ONLY.$SPECIMEN_NAME.$find_depth_and_alignments.bam


                         #SOFT CLIP FROM BED FILE
                         java -jar /Users/joelbarratt/jvarkit/dist/pcrclipreads.jar --interval $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed  \
                         $tmp_directory/$SPECIMEN_NAME.MAPPING/OVERLAP_ONLY.$SPECIMEN_NAME.$find_depth_and_alignments.bam --samoutputformat BAM \
                         -o $tmp_directory/$SPECIMEN_NAME.MAPPING/SOFT_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam

                         ##REMOVES SOFT CLIPPED BASES
                         java -jar /Users/joelbarratt/jvarkit/dist/biostar84452.jar $tmp_directory/$SPECIMEN_NAME.MAPPING/SOFT_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam \
                         > $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam
 
                         
                         #extract only mapped reads
                         ##### FIX THIS UP
                         samtools view -h -F 4 $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam \
                         | samtools view -bS > $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.sam  #### take mapped reads only
                         samtools view $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.sam \
                         | awk '{print("@"$1"\n"$10"\n+\n"$11)}'  > $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq


                        ##### EXTRACT COVERAGE INFORMATION --- THIS INFORMATION WILL BE USED TO DETERMINE DYNAMIC COVERAGE THRESHOLDS
                        samtools sort -T -n $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam \
                        -o $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only_aln_sorted.bam
                        samtools index $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only_aln_sorted.bam
                        samtools depth -aa -d 0 $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only_aln_sorted.bam \
                        | awk "/^$marker_range/" | head -$right_trim | awk "NR > $left_trim" \
                         > $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.coverage # exports a file containing the coverage
                        average_coverage=`awk -v N=3 '{ sum += $N } END { if (NR > 0) print sum / NR }' < $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.coverage`
                        echo "scale=0; ($average_coverage*0.10)/1" | bc > $tmp_directory/$SPECIMEN_NAME.COVERAGE_CUTOFF_PER_LOCUS_FOR_HAP_CALLING/$SPECIMEN_NAME.$find_depth_and_alignments.COVERAGE_CUTOFF
                        coverage_cutoff=`echo "scale=0; ($average_coverage*0.25)/1" | bc`
                         
#####MMASSIVE CD-HIT ERROR IF SOME SEQUENCES HAVE THE SAME NAME SO AT THIS STEP WE CAN RENAME THE FASTQ FILE HEADERS TO OVERCOME THIS.





#CD-HIT
#FLAGS

# c = % similarity of bins
# g =  1 or 0, default 0 -- explanation below:
#By cd-hit’s default algorithm, a sequence is clustered to the first
#cluster that meet the threshold (fast mode). If set to 1, the program
#will cluster it into the most similar cluster that meet the threshold
#(accurate but slow mode)
# d = length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space
#-p	1 or 0, default 0 --- if set to 1, print alignment overlap in .clstr file
# -s 1  <-- it is absolutely essential that this value is equal to 1.


$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq \
-o $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fq -p 1 -c 1 -g 1 -d 0 -T 10 -s 1


####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
seqret -sequence $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fq -outseq $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fasta



###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT
perl $working_directory/make_multi_seq.pl $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fasta \
     $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fq.clstr multi-seq $coverage_cutoff  ###dynamic cutoff.


cp -r $tmp_directory/$SPECIMEN_NAME.MAPPING/multi-seq $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq
rm -rf $tmp_directory/CONSENSUS_REFERENCES_FOR_READ_TRIMMING_AFTER_MAPPING.fasta $tmp_directory/$SPECIMEN_NAME.MAPPING/multi-seq


#### you are still in the $SPECIMEN_NAME.MAPPING directory



for OTHER_HAPS in `ls -1 $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq`
 do

# hASH the line below out when done
#OTHER_HAPS=`ls -1 $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq`


##############################################################################################################
# make a database of known non-junction markers to see if there is a match
#rm ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
cat $complete_reference_database > $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta

file_check_3=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

if [ `cat $file_check_3 |wc -l` -gt 0 ]; 
then 
cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
/usr/local/bin/gsed -i '/Junction/,+1 d' $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta # HASHED OUT BECAUSE ONLY RELEVANT TO CYCLOSPORA
fi



#############################################################################################################

# NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in /$tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS -input_type fasta -dbtype nucl


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
-query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
-word_size 7 \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result \
-max_target_seqs 1 \
-dust no \
-soft_masking false \
-outfmt "6 qseqid"


match_present=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result | wc -l`


## THE NEXT BLAST SEARCH JUST DOES A MORE STRINGENT CHECK OF COVERAGE FOR THE NEW HAPLOTYPE.

# convert the original clipped fastq into a fasta
seqret -sequence $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq -outseq $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fasta

# blast this to the potentially new haplotype with 100% identity and decent coverage
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
-query $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fasta \
-word_size 7 \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-num_threads 10 \
-out $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result_COVERAGE_GUESS \
-dust no \
-soft_masking false \
-outfmt "6 qseqid"

coverage_estimate_by_blast=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result_COVERAGE_GUESS | grep "Read_" | wc -l`

# count number of hits and this must be greater than 30 and greater than or equal to "REQUIRED_DEPTH_NEW_HAP"

# $REQUIRED_DEPTH_NEW_HAP <- this is the minimum requirement that is always user defined. Default is at 30
# $coverage_cutoff <- this is the sliding cutoff that changes with depth. We need to make sure that the depth is greater than the sliding cutoff and the minimum cutoff before a new haplotype is called.

if [ $match_present == 0 ] && [ $coverage_estimate_by_blast -gt $REQUIRED_DEPTH_NEW_HAP ] && [ $coverage_estimate_by_blast -gt $coverage_cutoff ];

then  $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
-db $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
-query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta -word_size 7 -evalue 0.001 \
-perc_identity 75 \
-qcov_hsp_perc 85 \
-num_threads 10 \
-out $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY \
-dust no \
-soft_masking false \
-max_target_seqs 1 \
-outfmt "6 qseqid"

which_marker=`/usr/local/bin/gsed 's/_Hap_.*/_Hap_/' $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY | head -1`
hap_count=`cat $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta | grep $which_marker | wc -l`
hap_count_plus_one=`echo "$hap_count + 1" | bc`

cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS | \
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | \
awk '/^>/{print ">'$which_marker$hap_count_plus_one'"; next}{print}' | awk NF > $working_directory/REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.$which_marker$hap_count_plus_one.fasta




#cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/19USxxxxx00D6PfD0000.chr5_mark1_part3_Hap_3.fasta | awk NF



elif [ $match_present -gt 0 ];

then echo "not a new haplotype"

fi


## cat all the FINAL_CLIPPED files to a new fastq file for mapping.
cat $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq >> $tmp_directory/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAP_CALLING.fastq

#

done ## completes blasting of the clusters mapping against a marker that meet the cutoff.    THIS LOOP SEEMS TO WORK PROPERLY NOW.









done ## completes clipping of reads and establishment of cutoffs for every individual marker -LOOP NOT WORKING.

cd $tmp_directory
rm -rf $SPECIMEN_NAME.MAPPING
## TRY TO HERE AND SEE IF IT DIES ----- THIS LOOP NOW WORKING PERFECTLY!








done ## completes running of the haplotype finding module for every specimen and haplotype

cd $tmp_directory

rm -rf multi-seq *.sam *.bam *BT_INDEX* *_CLUSTERS* *.amb *.ann *.bwt *.pac *.sa *blast_result* *.coverage *.blast_result *coverage_less_than* *all_clusters* *mapped_only* *.REMOVE_BAD_SEQUENCES.fasta SPECIMENS_TO_SEARCH_OTHER *.aln.sorted.*
















