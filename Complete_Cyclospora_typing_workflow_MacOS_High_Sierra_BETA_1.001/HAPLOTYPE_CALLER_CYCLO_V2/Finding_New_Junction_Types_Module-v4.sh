#!/bin/bash
#####################################################################################################################################################################

#UP November 19, 2020



working_directory=`cat RUN_DIR`
complete_reference_database=`cat ALL_REF`
temp_folder=`cat TMP_FOL`
junction_with_primer=`cat READ_REC_JUNCTION`
RAM_NEEDED=`cat RAM_ALLOCATION`
REQUIRED_DEPTH_NEW_HAP=`cat DEPTH_NEW`
number_of_threads=`cat THREADS_TO_USE`


cp ../Junction_manifest_file.txt ../MOD_Junction_manifest_file.txt


sed -i "" "s/parameters = -GE:not=12/parameters = -GE:not=$number_of_threads/g" $working_directory/MOD_Junction_manifest_file.txt

mkdir CONFIRMED_JUNCTIONS

cd $temp_folder/PROCESSED_READS/


     ##MAKE LIST OF SPECIMENS TO SEARCH FOR NEW JUNCTION SEQUENCES
     for file in *
     do
     echo "$file" >> $temp_folder/SPECIMENS_TO_SEARCH_JUNCTION
     done
     cd $temp_folder

gsed -i 's/.clean_merged.fastq//g' $temp_folder/SPECIMENS_TO_SEARCH_JUNCTION


##### This massive loop repeats for all specimens

             for SPECIMEN_NAME in `cat $temp_folder/SPECIMENS_TO_SEARCH_JUNCTION`

                         do

                                 rm JUNCTION_REFS.fasta

                                 cat $junction_with_primer > $temp_folder/JUNCTION_REFS.fasta # pull only junction sequences from original reference file.

                                 file_check=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*Junction*.fasta"`

                                       if [ `cat $file_check |wc -l` -gt 0 ]; 
                                       then 
                                       cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*Junction*.fasta >> $temp_folder/JUNCTION_REFS.fasta
                                       fi
                                 
                                 bwa index JUNCTION_REFS.fasta
                                 bwa mem -t $number_of_threads JUNCTION_REFS.fasta $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
                                 samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
                                 samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $SPECIMEN_NAME.mapped_only.fastq

                                 ## convert the fastq sequence just created to a fasta
                                 $working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/miraconvert $SPECIMEN_NAME.mapped_only.fastq $SPECIMEN_NAME.mapped_only.fasta

                                 ## make those fasta sequences a single line - this is necessary for the string searches that will follow
                                 cat $SPECIMEN_NAME.mapped_only.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' >> $SPECIMEN_NAME.mapped_only_single_line.fasta
                                 awk 'NF > 0' $SPECIMEN_NAME.mapped_only_single_line.fasta > $SPECIMEN_NAME.mapped_only_single_line2.fasta

               
                                ## collect the reads that were mapped to the junction sequences and then split them into their own fasta files in the directory "mapped_only"

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




                              rm -rf ../$SPECIMEN_NAME.*.sam
                              rm -rf ../$SPECIMEN_NAME.mapped_only.*


                                    READS_MAPPING_TO_YOUR_JUNCTION=`ls -1 | wc -l`



                                    ## Search for the forward and reverse primer sequences in every read that mapped to the junction sequence 
                                    ## so that these can be used as a refernce for a mapping assembly (MIRA).

                                    for my_fastas in *
                                                do

                                                      grep -l 'CCATCTACAGC.*AACAC' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
                                                      grep -l 'GTGTT.*GCTGTAGATGG' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt

                                                      grep -l 'TCTACAGC.*AACACGATC' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
                                                      grep -l 'GATCGTGTT.*GCTGTAGA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt

                                                done


                                    sort -u ../$SPECIMEN_NAME.list_of_sequences_to_review.txt > ../$SPECIMEN_NAME.2.list_of_sequences_to_review.txt






                                    mkdir ../assemble_these



                                      for READS in `cat ../$SPECIMEN_NAME.2.list_of_sequences_to_review.txt`
                                           do
                                              cat $READS >> $SPECIMEN_NAME.assemble_these_seqs
                                           done



                                    cp $SPECIMEN_NAME.assemble_these_seqs ../assemble_these
                                    cd ../assemble_these




                                 awk '/^>/{print ">Contig_" ++i; next}{print}' < $SPECIMEN_NAME.assemble_these_seqs > $SPECIMEN_NAME.cluster_these_seqs.fasta



                                 ## CLUSTER THE CLUSTERS PREVIOUSLY FOUND BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
                                 $working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est \
                                 -i $SPECIMEN_NAME.cluster_these_seqs.fasta \
                                 -o $temp_folder/assemble_these/$SPECIMEN_NAME.junction_clusters.fasta \
                                 -c 1 -g 1 -d 0 -T $number_of_threads -M $RAM_NEEDED


                                 $working_directory/CAP3/cap3.macosx.intel64/cap3 \
                                 $temp_folder/assemble_these/$SPECIMEN_NAME.junction_clusters.fasta \
                                 -o 95 -p 99.99

                                 cp  $SPECIMEN_NAME.junction_clusters.fasta.cap.contigs ../$SPECIMEN_NAME.map_to_this.fasta
                                 cat $SPECIMEN_NAME.junction_clusters.fasta.cap.singlets >> ../$SPECIMEN_NAME.map_to_this.fasta


                                 cd ..


                                 rm -rf mapped_only *mapped_only_single_line.fasta assemble_these *list_of_sequences_to_review


                                 cat $SPECIMEN_NAME.map_to_this.fasta | awk '/^>/{print ">new_junction_" ++i "_sequence"; next}{print}' > $SPECIMEN_NAME.map_to_this_2.fasta

                                 seqtk seq -F '#' $SPECIMEN_NAME.map_to_this_2.fasta > $SPECIMEN_NAME.map_to_this.fastq

                                 cp $SPECIMEN_NAME.map_to_this.fastq map_to_this.fastq
                                 cp $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq clean_merged.fastq

                                 #MIRA assembly

                                 $working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/mira /$working_directory/MOD_Junction_manifest_file.txt

                                 cp $temp_folder/novel_junction_finder_assembly/novel_junction_finder_d_results/novel_junction_finder_out_new_junc.unpadded.fas* .

                                 $working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/miraconvert novel_junction_finder_out_new_junc.unpadded.fasta $SPECIMEN_NAME.newly_found_junction.fastq

                                 rm -rf novel_junction_finder_assembly novel_junction_finder_out* *map_to_this.fastq




                                 ## CUTADAPT to trim all junction haplotypes to the primer

                                 cutadapt $SPECIMEN_NAME.newly_found_junction.fastq -g ACAGC \
                                 --output $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq
                                 cutadapt $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACAC \
                                 --output $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq
                                 cutadapt $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g GTGTT \
                                 --output $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq
                                 cutadapt $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT \
                                 --output $SPECIMEN_NAME.REV_three_prime_trim.fasta  --fasta


                                 cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_0_junction_sequences.fasta



                                 mkdir $temp_folder/CONFIRMED_JUNCTIONS/artefact_remove
                                 cd $temp_folder/CONFIRMED_JUNCTIONS/artefact_remove



                                                        while read line
                                                               do
                                                                   if [[ ${line:0:1} == '>' ]]
                                                                   then
                                                                       outfile=${line#>}.fasta
                                                                       echo $line > $outfile
                                                                   else
                                                                       echo $line >> $outfile
                                                                   fi
                                                               done < ../$SPECIMEN_NAME.DRAFT_0_junction_sequences.fasta

                                 ## Searches for forward and reverse primer sequences in each new contig, to make sure that the mapping assembly has not generated odd artifacts.

                                 for my_fastas in *
                                          do

                                               if 

                                                  [ `grep -il 'TGCGGAAACTGTATTTTTA.*TAAAAATTTAGTACACCTAGCC' $my_fastas | wc -l` -gt 0 ] || \
                                                  [ `grep -il 'GGCTAGGTGTACTAAATTTTTA.*TAAAAATACAGTTTCCGCA' $my_fastas | wc -l` -gt 0 ] || \

                                                  [ `grep -il 'TTATTTAATTTTACTATTTTAAAT.*ATTTAGTACACC' $my_fastas | wc -l` -gt 0 ] || \
                                                  [ `grep -il 'GGTGTACTAAAT.*ATTTAAAATAGTAAAATTAAATAA' $my_fastas | wc -l` -gt 0 ] || \

                                                  [ `grep -il 'GAAACTGTATTTTTA.*TAAAAATTTAGTACACCT' $my_fastas | wc -l` -gt 0 ] || \
                                                  [ `grep -il 'AGGTGTACTAAATTTTTA.*TAAAAATACAGTTTC' $my_fastas | wc -l` -gt 0 ] || \

                                                  [ `grep -il 'AATTTTACTATTTTAAAT.*ATTTAGTACA' $my_fastas | wc -l` -gt 0 ] || \
                                                  [ `grep -il 'TGTACTAAAT.*ATTTAAAATAGTAAAATT' $my_fastas | wc -l` -gt 0 ];

                                               then cat $my_fastas >> ../$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta

                                               fi

                                          done




                                               cd $temp_folder

                                               rm -rf $SPECIMEN_NAME.CUTADAPT_* $SPECIMEN_NAME.newly_found_junction.fastq $SPECIMEN_NAME.REV_three_prime_trim.fasta clean_merged*



####FIRST PASS - the following steps are repeated for several passes in to completely validate any junction sequences found.
##########################################################################################################################################
                                               $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta_BT_INDEX

                                               $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta_BT_INDEX \
                                               -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads \
                                               --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam

                                               ##### SORT BAM BEFORE MAKING CONSENSUS
                                               samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam \
                                               -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam
                                               samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam

                                               # call variants
                                               bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam | bcftools call -mv -Oz \
                                               -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

                                               bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

                                               cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
                                               | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz > \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.consensus.fasta

                                               # CLUSTER
                                               $working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
                                               -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads -M $RAM_NEEDED

                                               $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX

                                               $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX \
                                               -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 \
                                               --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam

                                               ##### SORT BAM BEFORE MAKING CONSENSUS
                                               samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam \
                                               -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam
                                               samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam

                                               # call variants
                                               bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam | bcftools call -mv -Oz \
                                               -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

                                               bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

                                               cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta \
                                               | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz > \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001.consensus.fasta


                                               seqtk seq -F '#' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001.consensus.fasta > \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fastq




                                               #### CUTADAPT TO TRIM OFF THE PRIMERS
                                               #### Note - this primer trimming step is intentionally repeated! When Illumina sequencing the junction, 
                                               #### a random artifact (a strange primer concatemer) is sometimes produced on either side of a true junction repeat.
                                               #### While it seems redundant, I run cutadapt multiple times in this workflow in order to remove these contcatemers. 
                                               #### Please do not remove these seemingly redundant steps. I did this on purpose!

                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fastq -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACAC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g GTGTT \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim.fasta  --fasta


                                               ## Do not remove these steps. They are here for a reason.
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim.fasta -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim2.fasta
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim2.fasta -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim3.fasta
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim3.fasta -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim4.fasta
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim4.fasta -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim5.fasta
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim5.fasta -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim6.fasta
                                               cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim6.fasta -g ACAGC \
                                               --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta




                                               ####Bowtie again
                                               ##########################################################################################################################################

                                               $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta_BT_INDEX


                                               $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX \
                                               -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 \
                                               --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT2.bam


                                               samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT2.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted2.bam
                                               samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted2.bam


                                               ## call variants
                                               bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted2.bam \
                                               | bcftools call -mv -Oz -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls2.vcf.gz
                                               bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls2.vcf.gz


                                               cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta \
                                               | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls2.vcf.gz > \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00002.consensus.fasta



                                               ## get your sequences on one line
                                               gsed -e 's/\(^>.*$\)/#\1#/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00002.consensus.fasta \
                                               | tr -d "\r" | tr -d "\n" | gsed -e 's/$/#/' | tr "#" "\n" | gsed -e '/^$/d' > \
                                               $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta





                                                          for CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00002.consensus.fasta | gsed 's/>//g'`
                                                                   do
                                                                            echo $CHECK_CONSENSUS_DRAFT_1 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt
                                                                            grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta > \
                                                                            $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta
                                                                            grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta > \
                                                                            $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta

                                                                            rm $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

                                                                                if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta | tail -1` == `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta | tail -1` ];
                                                                                    then gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                                                                fi
                                                                    done





                                   ##### The following steps take the putative new contigs through a series of logical routines that allow the 
                                   ##### workflow to decide if any contig produced in the first Bowtie pass is real or an artifact.
                                   ##### They essentially check for contigs that have N bases in their sequence

                                   for CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001.consensus.fasta | gsed 's/>//g'`

                                      do
                                         echo $CHECK_CONSENSUS_DRAFT_1 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

                                         grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta > \
                                         $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_CONSENSUS.fasta ## sequence obtained after mapping using bowtie (a consensus)
                                         grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta > \
                                         $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_BEFORE_consensus.fasta ### sequence obtained from the original cap3 assembly - before alignment.
                                         rm $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

                                         ## The following loops logically check whether the sequence of both of the above fasta files are sane, if one is sane, or neither are sane.

                                         ## does the consensus contain any N bases? If it does, then delete it from the pile of potential junction sequences.
                                         if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_CONSENSUS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];

                                              then gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                         fi


                                         ## does the original contig before bowtie mapping have any N bases? If it does, then delete it from the pile of potential junction sequences.
                                         if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_BEFORE_consensus.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];
                                              then gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                         fi


                                         ## Check results of the above two "if" loops, and if the two conditions below are met, scrap that cap3 contig and the consensus for it 
                                         ## completely from the list of potential junctions.
                                         if [ `gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta | tail -1 | wc -l` -eq 0 ] && \
                                            [ `gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];
                                              then gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta
                                         fi

                                         if [ `gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ] && \
                                            [ `gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];
                                         then gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta
                                         gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                         fi
                                       done






                                       ## Additional sanity check. Any sequences left over (before and after correction of cap3 contigs by mapping) are BLASTed to see if they
                                       ## had a decent BLAST hit before a consensus was generated and if the same hit is obtained after the consensus was generated.

                                    for SCND_CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta | gsed 's/>//g'`
                                         do
                                            echo $SCND_CHECK_CONSENSUS_DRAFT_1 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

                                            grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta
                                            grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta


                                            $working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
                                            -in $temp_folder/JUNCTION_REFS.fasta \
                                            -input_type fasta \
                                            -dbtype nucl

                                            ## Please do not turn on dust or soft masking. *TRUST ME*
                                            $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $temp_folder/JUNCTION_REFS.fasta \
                                            -query $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta \
                                            -word_size 7 \
                                            -evalue 0.001 \
                                            -perc_identity 100 \
                                            -qcov_hsp_perc 100 \
                                            -num_threads $number_of_threads \
                                            -out $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result \
                                            -max_target_seqs 1 \
                                            -dust no \
                                            -soft_masking false \
                                            -outfmt "6 qseqid qseq"


                                            $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $temp_folder/JUNCTION_REFS.fasta \
                                            -query $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta \
                                            -word_size 7 \
                                            -evalue 0.001 \
                                            -perc_identity 100 \
                                            -qcov_hsp_perc 100 \
                                            -num_threads $number_of_threads \
                                            -out $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result \
                                            -max_target_seqs 1 \
                                            -dust no \
                                            -soft_masking false \
                                            -outfmt "6 qseqid qseq"


                                            ## did they obtain the same blast hit?
                                            if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` == `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` ];
                                            then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                            fi

                                            ## did the consensus receive a blast hit but the original sequence not?
                                            if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` -eq 0 ] && \
                                               [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -gt 0 ];
                                            then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                            fi

                                            ## did neither obtain a blast hit? If this is true, you may have discovered a new junction haplotype.
                                            if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` -eq 0 ] \
                                            && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -eq 0 ];
                                            then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta
                                            fi

                                            ## did the original cap3 contig NOT obtain a blast hit, but the consensus did?
                                            if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -eq 0 ] \
                                            && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` -gt 0 ];
                                            then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta
                                            fi
                                         done


                                            
                                            ## you now have a set of sequences that passed the first sanity check. Write these to file.
                                            cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta >> \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.01_junction_sequences.fasta
                                            cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta >> \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.01_junction_sequences.fasta

                                            
                                            ## Tidy up the remaining junk files before moving on to the next stages.
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta_BT_*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta.fai
                                            rm -rf $temp_folder/CONFIRMED_JUNCTIONS/artefact_remove
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT2.bam
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1*




                                            ########################################################################################################################################################
                                            # Run a clustering step. This will check that none of the junction sequences remaining in your list are identical to each other. If they are, these will
                                            # merged in this clustering step before moving on to the next pass.
                                            $working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.01_junction_sequences.fasta \
                                            -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads -M $RAM_NEEDED
                                            ########################################################################################################################################################



                                            ## Convert the remaining fasta clusters (from the previous clustering step) to a fastq file with maximum qualtiy before I feed back into cutadapt.
                                            seqtk seq -F '#' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fasta > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fastq


                                            ## CUTADAPT TO TRIM OFF THE PRIMERS - please do this again. Trust me. Cutadapt is repeated here multiple times for a very good reason. This is not a mistake.

                                            cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fastq -g ACAGC \
                                            --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq
                                            cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACAC \
                                            --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq
                                            cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g GTGTT \
                                            --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq
                                            cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT \
                                            --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.1_junction_sequences.fasta  --fasta


                                            ######################################################################################################################################################
                                            # Cluster again hoping that after the cutadapt step, we may have made some fasta sequences that now look identical to each other and should be merged before
                                            # the next pass.
                                            $working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.1_junction_sequences.fasta \
                                            -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads -M $RAM_NEEDED
                                            ######################################################################################################################################################





##################### SECOND PASS
### this is where we generate draft 2 of the junction sequences. The previous pass had the first set of draft junction sequences.


                                            ## BOWTIE mapping to generate a consensus sequence from the contigs that survived the first pass.
                                            $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta_BT_INDEX

                                            $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta_BT_INDEX \
                                            -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.OUTPUT.bam


                                            ## EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
                                            samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.aln.sorted.bam
                                            samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.aln.sorted.bam
                                            samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.aln.sorted.bam > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.coverage # exports a file containing the coverage.

                                   
                                   ## Find contigs with coverage less than 10 at **ANY** base across the entire contig.         
                                   for JUNK_REMOVAL in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.coverage | cut -f 1 | uniq`
                                            do
         if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.coverage | gsed 's/^/>/' | gsed "s/$JUNK_REMOVAL\t/$JUNK_REMOVAL\t/g" | grep "$JUNK_REMOVAL" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt 10 ]
                 then echo $JUNK_REMOVAL >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.contigs_with_coverage_less_than.10.txt
         fi
                                            done


                                            cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta
                                            gsed -i 's/^/>/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.contigs_with_coverage_less_than.10.txt ## print the names of contigs that have a low coverage base.


                                            ## Delete these low coverage contigs from your set of potential junctions using the previously generated list as a reference.
                                            for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.contigs_with_coverage_less_than.10.txt`
                                                  do
                                                     gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta # removes line containing the fasta seq name. 
                                                                                                                                                              # The '+1' removes the line after as well.
                                                  done


                                            # Tidy up remaining junk files
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_0*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT* 
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_H*
                                            rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1*


                                            ## Yep, Bowtie again. This is also intentional - not a mistake. We slowly but surely discard contigs that seem less likely to be real with every iteration
                                            ## As contigs are steadily discarded and we map to the remaining contigs, the coverage obtained for the remaining contigs increases.
                                            ## We can therefore steadily increase our quality as we move towards a consensus for the true junction sequence in these specimens.
                                            $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_consensus.fasta_BT_INDEX

                                            $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_consensus.fasta_BT_INDEX \
                                            -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q --threads $number_of_threads --local > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED.OUTPUT.bam

                                            ## SORT BAM BEFORE MAKING CONSENSUS
                                            samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED.OUTPUT.bam \
                                            -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_sorted.bam
                                            samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_sorted.bam

                                            # call variants
                                            bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_sorted.bam | bcftools call -mv -Oz \
                                            -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.calls3.vcf.gz
                                            bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.calls3.vcf.gz
                                            cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta | bcftools consensus \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.calls3.vcf.gz > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.RELAXED.consensus.fasta



                                            ## Do the blasting again for a sanity check of the remaining contigs
                                            ## get your sequences on one line
                                            gsed -e 's/\(^>.*$\)/#\1#/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.RELAXED.consensus.fasta \
                                            | tr -d "\r" | tr -d "\n" | gsed -e 's/$/#/' | tr "#" "\n" | gsed -e '/^$/d' > \
                                            $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta




                                            for SCND_CHECK_CONSENSUS_DRAFT_2 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta | gsed 's/>//g'`
                                                  do
                                                      echo $SCND_CHECK_CONSENSUS_DRAFT_2 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt
                                                      grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta > \
                                                      $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.fasta
                                                      grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta > \
                                                      $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.fasta
                                                      rm $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

                                                      $working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
                                                      -in $temp_folder/JUNCTION_REFS.fasta \
                                                      -input_type fasta \
                                                      -dbtype nucl

                                                      $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $temp_folder/JUNCTION_REFS.fasta \
                                                      -query $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.fasta \
                                                      -word_size 7 \
                                                      -evalue 0.001 \
                                                      -perc_identity 100 \
                                                      -qcov_hsp_perc 100 \
                                                      -num_threads $number_of_threads \
                                                      -out $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result \
                                                      -max_target_seqs 1 \
                                                      -dust no \
                                                      -soft_masking false \
                                                      -outfmt "6 qseqid qseq"


                                                      $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db $temp_folder/JUNCTION_REFS.fasta \
                                                      -query $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.fasta \
                                                      -word_size 7 \
                                                      -evalue 0.001 \
                                                      -perc_identity 100 \
                                                      -qcov_hsp_perc 100 \
                                                      -num_threads $number_of_threads \
                                                      -out $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result \
                                                      -max_target_seqs 1 \
                                                      -dust no \
                                                      -soft_masking false \
                                                      -outfmt "6 qseqid qseq"



                                                      if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` == `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` ];
                                                              then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
                                                      fi

                                                      if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -eq 0 ] && \
                                                         [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -gt 0 ];
                                                            then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
                                                      fi


                                                      if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -eq 0 ] && \
                                                         [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -eq 0 ];
                                                                then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
                                                      fi


                                                      if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -eq 0 ] && \
                                                         [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -gt 0 ];
                                                                 then gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta
                                                      fi
                                                 done


                                                 ## Write remaining contigs to file for the next pass.
                                                 cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX.fasta
                                                 cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX.fasta
                                                 cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta




################################################## THIRD PASS FOR IDENTIFYING THE NEXT JUNCTION -- nearly done.... The junction repeat is a great marker, but hard to detect and validate (its a repeat).

                                                 $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta_BT_INDEX

                                                 $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta_BT_INDEX \
                                                 -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.OUTPUT.bam

                                                 ## EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
                                                 samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.OUTPUT.bam -o \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.aln.sorted.bam
                                                 samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.aln.sorted.bam

                                                 samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.aln.sorted.bam > \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.coverage # exports a file containing the coverage.

                                                 




                                                 ## remove any contig that has a base with coverage that is not above 1%
                                                 JUNCTION_rough_coverage_cutoff=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.010)/1" | bc`

                                                 for JUNK_REMOVAL in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.coverage | cut -f 1 | uniq`
                                                     do
if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.coverage | gsed 's/^/>/' | gsed "s/$JUNK_REMOVAL\t/$JUNK_REMOVAL\t/g" | grep "$JUNK_REMOVAL" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff ]
      then echo $JUNK_REMOVAL >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt
fi
                                                     done

                                                 gsed -i 's/^/>/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt
                                                 cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.consensus.fasta





                                                 ## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST OF NEW ONES THAT POSSESS LOW COVERAGE
                                                 for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt`
                                                      do
                                                 gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.consensus.fasta
                                                      done

                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2*





                                                 ## Cluster again. This is because the increased coverage obtained by discarding low-coverage contigs previously might result in some
                                                 ## previously identified contigs becoming identical to each other.
                                                 #########################################################################################################################################
                                                 # WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
                                                 $working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.consensus.fasta \
                                                 -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta -c 1 -g 1 -d 0 -T $number_of_threads -M $RAM_NEEDED
                                                 #########################################################################################################################################


                                                 ## Clean up
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX*


                                                 ## Bowtie again.
                                                 $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta_BT_INDEX

                                                 $working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta_BT_INDEX \
                                                 -U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.OUTPUT.bam

                                                 ## EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
                                                 samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.bam
                                                 samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.bam

                                                 samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.bam > \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage # exports a file containing the coverage.

                                                 ## slight increase in coverage cutoff to filter out any contigs that may still have low coverage
                                                 JUNCTION_rough_coverage_cutoff_2=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.015)/1" | bc` ####### DO NOT CHANGE THIS CUTOFF -- TRUST ME.




                                                 ### find out which haplotypes contain bases with coverage less than your default coverage cutoff and the coverage cutoff identified above.
                                                 for NEW_HAPS in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | cut -f 1 | uniq` 
                                                       do
          if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | gsed 's/^/>/' | gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $REQUIRED_DEPTH_NEW_HAP ] \
          && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | gsed 's/^/>/' | gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff_2 ];
               then echo $NEW_HAPS >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.$REQUIRED_DEPTH_NEW_HAP.txt
          fi
                                                       done


                                                 ## generates list of sequence names that failed to meet the cutoffs identified above.
                                                 gsed -i 's/^/>/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.$REQUIRED_DEPTH_NEW_HAP.txt





                                                 ## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST THAT POSSESS LOW COVERAGE
                                                 for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.$REQUIRED_DEPTH_NEW_HAP.txt`
                                                      do
                                                 gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta
                                                      done

                                                 cp $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.validated_junction_sequences.fasta 
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3*




                                            ###### Now we look for sharp drops and rises in coverage within each remaining contig. 
                                            ###### If there are sharp drops and rises of more than 2 standard deviations above the mean, the contig is deleted as it is probably false.
                                            ###### Recall that the complete sequence of the junction is pretty short (it should fit within a single pair of Illumina reads, 
                                            ###### so the coverage across a junction sequence should be pretty consistent.


                                                 for JUNCTION_TESTER in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | cut -f 1 | uniq` 
                                                       do

                                                 coverage_at_each_base=`cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | gsed -n "/$JUNCTION_TESTER/p" | cut -f 3`
                                                 number_of_bases=`echo $coverage_at_each_base | wc -w`
                                                 sum_of_bases=`echo $coverage_at_each_base | xargs | sed -e 's/\ /+/g' | bc`
                                                 average_coverage=`echo "($sum_of_bases/$number_of_bases)" | bc`

                                                 standardDeviation=$(
                                                     echo "$coverage_at_each_base" |
                                                         awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}'
                                                                     )

                                                 two_SDs=`echo "(2*$standardDeviation)" | bc`

                                                 ### now if the average coverage is less than 2 standard deviations from the average coverage, then delete the contig.
                                                 ### It means that there is a sharp drop in coverage somewhere, indicative of some erroneous bases - the Contig is probably not a true sequence

                                                 rounded_2SDs=`echo ${two_SDs%%.*}`

                                                 if [ $rounded_2SDs -gt $average_coverage ];

                                                 ### delete the contig in question
                                                 then gsed -i "/$JUNCTION_TESTER/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.validated_junction_sequences.fasta

                                                 fi

                                                        done


                                                 
                                                 ## Tidy up junk files left over from previous passes.

                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.OUTPUT.bam
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.ba*
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage
                                                 rm -rf $temp_folder/CONFIRMED_JUNCTIONS/artefact_remove



                                                 cd $temp_folder



                                                 ## you are now done confirming the sequence of the mt junctions that you have found. All remaining junction sequences are considered valid.
                                                 ## Write them to a final file.
                                                 cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.validated_junction_sequences.fasta > \
                                                 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FINAL.validated_junction_sequences.fasta



                                                 # This will split the newly discovered junction haplotypes into separate fasta files if there are multiple junctions found.

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
                                                     done < $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FINAL.validated_junction_sequences.fasta



                                                 ## Clean up remaining junk files
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4*
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage*
                                                 rm $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2*


                                                 #### Check if the remaining valid junction types are new haplotypes, and if they are, write them to our new haplotypes folder
                                                 #### This part of the code will make sure that their names are correct (name of our junction types is based on their length).
                                                 for FOOBAR in `ls -1`
                                                    do 


                                                    junction_count=`grep 'Junction_Hap' $temp_folder/JUNCTION_REFS.fasta | wc -l`
                                                    new_number=`echo "$junction_count + 1" | bc`
                                                    cat $FOOBAR | awk '/^>/{print ">Junction_Hap_'$new_number'"; next}{print}' > \
                                                    $SPECIMEN_NAME.junction_recount_$FOOBAR

                                                    # IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL TO AND NOT OBTAINING 100% 
                                                    # COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE WRITE IT TO OUR LIST OF NEW HAPS.
                                                    $working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
                                                    -in $temp_folder/JUNCTION_REFS.fasta \
                                                    -input_type fasta \
                                                    -dbtype nucl \
                                                    -title junction_check


                                                    $working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db ../JUNCTION_REFS.fasta \
                                                    -query $SPECIMEN_NAME.junction_recount_$FOOBAR \
                                                    -word_size 7 \
                                                    -evalue 0.001 \
                                                    -perc_identity 100 \
                                                    -qcov_hsp_perc 100 \
                                                    -num_threads $number_of_threads \
                                                    -out $temp_folder/new_potential_junctions/$SPECIMEN_NAME.RESULT_$FOOBAR.blast_result \
                                                    -max_target_seqs 1 \
                                                    -dust no \
                                                    -soft_masking false \
                                                    -outfmt "6 qseqid pident mismatch gapopen gaps sseqid sseq evalue bitscore"

                                                    match_present=`cat $temp_folder/new_potential_junctions/$SPECIMEN_NAME.RESULT_$FOOBAR.blast_result | wc -l`

                                                    if [ $match_present == 0 ];
                                                       then echo "you have discovered a new junction type"

                                                       repeat_length=`cat $SPECIMEN_NAME.junction_recount_$FOOBAR | sed "1d" | tr -cd '[:alpha:]' | wc -m`
                                                       rep_length_2=`echo "$repeat_length + 21 + 22" | bc` ## Number of bases in the new haplotype plus the length of each primer sequence. 
                                                                                                           ## These are the primers in the first Mt Junction paper by Nascimento et al.

                                                       gsed -i "1s/.*/>Mt_Cmt$rep_length_2.X_Junction_Hap_$new_number/" $SPECIMEN_NAME.junction_recount_$FOOBAR
                                                       cat $SPECIMEN_NAME.junction_recount_$FOOBAR >> $temp_folder/JUNCTION_REFS.fasta
                                                       cat $SPECIMEN_NAME.junction_recount_$FOOBAR >> \
                                                       $working_directory/REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.Mt_Cmt$rep_length_2.X_Junction_Hap_$new_number.fasta

                                                    elif [ $match_present -gt 0 ];
                                                    then echo "no new junction sequences were found"
                                                    fi

                                                       done

                                                    cd ..

                                                    ## tidy up of junk sequences
                                                    rm -rf new_potential_junctions
                                                    rm -rf *fasta.amb *fasta.ann *fasta.bwt *fasta.nhr *fasta.nin *fasta.nsq 
                                                    rm -rf *fasta.pac *fasta.sa *junction_recount_* *review.txt *RESULT* *novel* *mapped_only* *map_to_this*


            done   ## this last massive mega loop run for each speicmen is done.


     rm $working_directory/MOD_Junction_manifest_file.txt # remove the copied Mira manifest
     cd $temp_folder



