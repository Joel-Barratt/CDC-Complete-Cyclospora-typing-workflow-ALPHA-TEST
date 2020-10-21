#!/bin/bash
#####################################################################################################################################################################

### need pysam installed. pip install pysam

working_directory=`cat RUN_DIR`
read_recovery_refs=`cat READ_REC` 
complete_reference_database=`cat ALL_REF`
tmp_directory=`cat TMP_FOL`
REQUIRED_DEPTH_NEW_HAP=`cat DEPTH_NEW`
bed_references=`cat BED_LOCATION`
number_of_threads=`cat THREADS_TO_USE`
RAM_NEEDED=`cat RAM_ALLOCATION`


cd $tmp_directory/PROCESSED_READS/
for file in *
             do
                echo "$file" >> $tmp_directory/SPECIMENS_TO_SEARCH_OTHER
             done
cd $tmp_directory

gsed -i 's/.clean_merged.fastq//g' SPECIMENS_TO_SEARCH_OTHER
cat $bed_references | cut -f4 > $tmp_directory/LIST_OF_LOCI_UNIQUE


mkdir $tmp_directory/CLIPPED_READS_FOR_HAPLOTYPE_CALLING


for SPECIMEN_NAME in `cat SPECIMENS_TO_SEARCH_OTHER`

                              do


                                           mkdir $tmp_directory/$SPECIMEN_NAME.COVERAGE_CUTOFF_PER_LOCUS_FOR_HAP_CALLING



                                           rm READ_RECOVERY_REFERENCE_SEQUENCES.fasta
                                           cat $read_recovery_refs > $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
                                           file_check_1=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

                                                      if [ `cat $file_check_1 |wc -l` -gt 0 ]; 
                                                      then 
                                                      cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
                                                      gsed -i '/Junction/,+1 d' $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta # Only relevant to Cyclospora.
                                                      fi


                                           #### MAP READS TO THE REFERENCE DATABASE OF NON-JUNCTION MARKERS TO OBTAIN CORRECT READS
                                           bwa index $tmp_directory/READ_RECOVERY_REFERENCE_SEQUENCES.fasta
                                           bwa mem -t $number_of_threads READ_RECOVERY_REFERENCE_SEQUENCES.fasta $tmp_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
                                           samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
                                           samtools view $SPECIMEN_NAME.mapped_only.sam \
                                           | awk '{print("@"$1"\n"$10"\n+\n"$11)}' \
                                           | awk '{if(NR%4==1) $0=sprintf("@Read_%d",(1+i++)); print;}' > $SPECIMEN_NAME.mapped_only.fastq  # the last awk section renames the reads 


                                           cd $tmp_directory
                                           cat $read_recovery_refs > $tmp_directory/FULL_LENGTH_MARKERS.fasta
                                           mkdir $SPECIMEN_NAME.MAPPING
                                           cd $SPECIMEN_NAME.MAPPING


                                           ###### MAP TO LONG (UNSPLIT) REFERENCES using Bowtie
                                           $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build \
                                           $tmp_directory/FULL_LENGTH_MARKERS.fasta \
                                           $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.BT_INDEX
                                           $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 \
                                           -x $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.BT_INDEX \
                                           -U $tmp_directory/$SPECIMEN_NAME.mapped_only.fastq \
                                           --threads $number_of_threads \
                                           --local > $SPECIMEN_NAME.alignment.bam


                                           samtools sort $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.alignment.bam \
                                           > $tmp_directory/$SPECIMEN_NAME.MAPPING/sorted.$SPECIMEN_NAME.alignment.bam

                                           samtools index sorted.$SPECIMEN_NAME.alignment.bam






                                           for find_depth_and_alignments in `cat $tmp_directory/LIST_OF_LOCI_UNIQUE`
 
                                                                                                      do

                                                           marker_range=`echo $find_depth_and_alignments | gsed 's/_PART_*.//g'`
                                          
                                                           cat $bed_references | grep $find_depth_and_alignments > $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed
                                                           left_trim=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed | cut -f2`
                                                           right_trim=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed | cut -f3`
                         

                                                           samtools view  -bu -F 4  sorted.$SPECIMEN_NAME.alignment.bam $marker_range | \
                                                           java -jar /Users/joelbarratt/jvarkit/dist/samjs.jar  \
                                                           -e "record.alignmentStart <= $left_trim && record.alignmentEnd >= $right_trim"  --samoutputformat BAM \
                                                           -o $tmp_directory/$SPECIMEN_NAME.MAPPING/OVERLAP_ONLY.$SPECIMEN_NAME.$find_depth_and_alignments.bam


                                                           #SOFT CLIP FROM BED FILE
                                                           java -jar /Users/joelbarratt/jvarkit/dist/pcrclipreads.jar --interval $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.bed  \
                                                           $tmp_directory/$SPECIMEN_NAME.MAPPING/OVERLAP_ONLY.$SPECIMEN_NAME.$find_depth_and_alignments.bam --samoutputformat BAM \
                                                           -o $tmp_directory/$SPECIMEN_NAME.MAPPING/SOFT_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam

                                                           ##REMOVES SOFT CLIPPED BASES
                                                           java -jar /Users/joelbarratt/jvarkit/dist/biostar84452.jar \
                                                           $tmp_directory/$SPECIMEN_NAME.MAPPING/SOFT_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam \
                                                           > $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.map_to_marker_only.bam
 
                         
                                                           #extract only mapped reads
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




                                                          ### The two lines below are used to determine the sliding cutoff for haplotype assignment to a specimen (are not used in the present workflow to decide if we want to write a new haplotype  
                                                          ### to our database. For the haplotype assignment workflow, want at least 20 reads OR greater than 10% of the entire coverage whichever of these two values is higher. I do this as part
                                                          ### of the current workflow to save time (saves me having to map reads again). So the results of the next two steps will be written to file now, just to be references during the haplotype assignment workflow.
                                                          average_coverage=`awk -v N=3 '{ sum += $N } END { if (NR > 0) print sum / NR }' < $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.coverage`
                                                          echo "scale=0; ($average_coverage*0.10)/1" | bc > $tmp_directory/$SPECIMEN_NAME.COVERAGE_CUTOFF_PER_LOCUS_FOR_HAP_CALLING/$SPECIMEN_NAME.$find_depth_and_alignments.COVERAGE_CUTOFF



                                                          ### The line below is for the present workflow. To write a new haplotype to the reference database, I want it to have at least 100 times coverage AND I want the new haplotype to be supported by 
                                                          ### at least 25% of the reads for that locus (whichever value is higher will be the final cutoff used).
                                                          coverage_cutoff=`echo "scale=0; ($average_coverage*0.25)/1" | bc`


                         
                                                          #####CD-HIT ERROR IF SOME SEQUENCES HAVE THE SAME NAME - SO I RENAME THE FASTQ FILE HEADERS TO OVERCOME THIS.
                                                          $working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est \
                                                          -i $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq \
                                                          -o $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fq \
                                                          -p 1 -c 1 -g 1 -d 0 -T $number_of_threads -s 1 -M $RAM_NEEDED


                                                          ####THIS IS AN EMBOSS TOOL THAT WILL CONVERT FASTQ TO FASTA - you need to do this in order to turn the clusters into fasta files
                                                          seqret -sequence $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fq \
                                                          -outseq $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fasta


                                                          ###### THIS NEXT STEP WILL FILTER CLUSTES BY SIZE - REMOVING CLUSTERS THAT DO NOT HAVE MORE THAN 10 SEQUENCES IN THEM --- STEP 6 OUTPUT
                                                          perl $working_directory/make_multi_seq.pl $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fasta \
                                                          $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.clean_merged_CLUSTERS.fq.clstr multi-seq $coverage_cutoff  ###dynamic cutoff.

                                                          cp -r $tmp_directory/$SPECIMEN_NAME.MAPPING/multi-seq $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq
                                                          rm -rf $tmp_directory/CONSENSUS_REFERENCES_FOR_READ_TRIMMING_AFTER_MAPPING.fasta $tmp_directory/$SPECIMEN_NAME.MAPPING/multi-seq







                                                                                              for OTHER_HAPS in `ls -1 $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq`
                                                                                                                                      do


                                                                                                                    # make a database of known non-junction markers to see if there is a match
                                                                                                                    cat $complete_reference_database > $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta

                                                                                                                    file_check_3=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta"`

                                                                                                                    if [ `cat $file_check_3 |wc -l` -gt 0 ]; 
                                                                                                                    then 
                                                                                                                    cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*.fasta >> $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta
                                                                                                                    gsed -i '/Junction/,+1 d' $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta # HASHED OUT BECAUSE ONLY RELEVANT TO CYCLOSPORA
                                                                                                                    fi



                                                                                                                    $working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
                                                                                                                    -in /$tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
                                                                                                                    -input_type fasta -dbtype nucl


                                                                                                                    $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
                                                                                                                    -db $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
                                                                                                                    -query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta \
                                                                                                                    -word_size 7 \
                                                                                                                    -evalue 0.001 \
                                                                                                                    -perc_identity 100 \
                                                                                                                    -qcov_hsp_perc 100 \
                                                                                                                    -num_threads $number_of_threads \
                                                          -out $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result \
                                                                                                                    -max_target_seqs 1 \
                                                                                                                    -dust no \
                                                                                                                    -soft_masking false \
                                                                                                                    -outfmt "6 qseqid"


match_present=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result | wc -l`



                                                                                                                    # convert the original clipped fastq into a fasta
                                                                                                                    seqret -sequence $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq \
                                                                                                                    -outseq $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fasta

                                                                                                                    # blast this to the potentially new haplotype with 100% identity and decent coverage
                                                                                                                    $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
                                                                                                                    -db $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
                                                                                                                    -query $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fasta \
                                                                                                                    -word_size 7 \
                                                                                                                    -evalue 0.001 \
                                                                                                                    -perc_identity 100 \
                                                                                                                    -qcov_hsp_perc 100 \
                                                                                                                    -num_threads $number_of_threads \
                                                          -out $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result_COVERAGE_GUESS \
                                                                                                                    -dust no \
                                                                                                                    -soft_masking false \
                                                                                                                    -outfmt "6 qseqid"

coverage_estimate_by_blast=`cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.$find_depth_and_alignments.RESULT_$OTHER_HAPS.blast_result_COVERAGE_GUESS | grep "Read_" | wc -l`

                                                          ## The BLAST search below is meant to be less than 100% - this is just for a last minute depth check.

                                                          if [ $match_present == 0 ] && [ $coverage_estimate_by_blast -gt $REQUIRED_DEPTH_NEW_HAP ] && [ $coverage_estimate_by_blast -gt $coverage_cutoff ];

                                                          then  $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn \
                                                          -db $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS \
                                                          -query $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta -word_size 7 -evalue 0.001 \
                                                          -perc_identity 75 \
                                                          -qcov_hsp_perc 85 \
                                                          -num_threads $number_of_threads \
                                                          -out $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY \
                                                          -dust no \
                                                          -soft_masking false \
                                                          -max_target_seqs 1 \
                                                          -outfmt "6 qseqid"

which_marker=`gsed 's/_Hap_.*/_Hap_/' $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$SPECIMEN_NAME.RESULT_$OTHER_HAPS.blast_result_REDUCED_SIMILARITY | head -1`
                                                          hap_count=`cat $tmp_directory/ALL_REFERENCE_SEQUENCES_FOR_BLASTING.fasta | grep $which_marker | wc -l`
                                                          hap_count_plus_one=`echo "$hap_count + 1" | bc`

                                                          cat $tmp_directory/$SPECIMEN_NAME.MAPPING/$SPECIMEN_NAME.$find_depth_and_alignments.multi-seq/$OTHER_HAPS | \
                                                          awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | \
                                                          awk '/^>/{print ">'$which_marker$hap_count_plus_one'"; next}{print}' \
                                                          | awk NF > $working_directory/REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.$which_marker$hap_count_plus_one.fasta

                                                          elif [ $match_present -gt 0 ];

                                                          then echo "not a new haplotype"

                                                          fi

                                                                                                                    ## cat all the FINAL_CLIPPED files to a new fastq file for mapping.
                                                                                                                    cat $tmp_directory/$SPECIMEN_NAME.MAPPING/FINAL_CLIPPED.$SPECIMEN_NAME.$find_depth_and_alignments.mapped_only.fastq >> \
                                                                                                                    $tmp_directory/CLIPPED_READS_FOR_HAPLOTYPE_CALLING/$SPECIMEN_NAME.FINAL_CLIPPED_FOR_HAP_CALLING.fastq

                                                                                                                                       done ## completes blasting of the clusters mapping against a marker that meet the cutoff.




                                                                                              done ## completes clipping of reads and establishment of cutoffs for every individual marker.

                                                          cd $tmp_directory
                                                          rm -rf $SPECIMEN_NAME.MAPPING


                                       done ## completes running of the haplotype finding module for every specimen and haplotype.

cd $tmp_directory

rm -rf multi-seq *.sam *.bam *BT_INDEX* *_CLUSTERS* *.amb *.ann *.bwt *.pac *.sa *blast_result* *.coverage *.blast_result *coverage_less_than* *all_clusters* *mapped_only* *.REMOVE_BAD_SEQUENCES.fasta SPECIMENS_TO_SEARCH_OTHER *.aln.sorted.*
















