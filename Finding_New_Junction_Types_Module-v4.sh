#!/bin/bash
#####################################################################################################################################################################
####     DISCOVER NEW MT JUNCTION TYPES    ##########################################################################################################################
#####################################################################################################################################################################

#cd TMP

#rm SPECIMENS_TO_SEARCH_JUNCTION

working_directory=`cat RUN_DIR` # /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER
complete_reference_database=`cat ALL_REF` # /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER/REF_SEQS/BLASTING/ORIGINAL_REFS/2019_ORIGINAL_REFERENCES_SHORTENED.fasta
temp_folder=`cat TMP_FOL` #/Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER/TMP
junction_with_primer=`cat READ_REC_JUNCTION`



REQUIRED_DEPTH_NEW_HAP=`cat DEPTH_NEW`

#REQUIRED_DEPTH_NEW_HAP=50


number_of_threads=`cat CORES_TO_USE`


cp ../Junction_manifest_file.txt ../MOD_Junction_manifest_file.txt


sed -i "" "s/parameters = -GE:not=12/parameters = -GE:not=$number_of_threads/g" $working_directory/MOD_Junction_manifest_file.txt


### used to be here

mkdir CONFIRMED_JUNCTIONS



cd $temp_folder/PROCESSED_READS/


##MAKE LIST OF SPECIMENS TO SEARCH FOR NEW JUNCTION SEQUENCES
for file in *
do
echo "$file" >> $temp_folder/SPECIMENS_TO_SEARCH_JUNCTION
done
cd $temp_folder

/usr/local/bin/gsed -i 's/.clean_merged.fastq//g' $temp_folder/SPECIMENS_TO_SEARCH_JUNCTION



for SPECIMEN_NAME in `cat $temp_folder/SPECIMENS_TO_SEARCH_JUNCTION`


       do


#SPECIMEN_NAME=`cat SPECIMENS_TO_SEARCH_JUNCTION | head -1`
#SPECIMEN_NAME=S_TX019_18


## or does it work best here?
rm JUNCTION_REFS.fasta

############################################################################################################## get all known junction sequences for read recovery
cat $junction_with_primer > $temp_folder/JUNCTION_REFS.fasta # pull only junction sequences from original reference file.

file_check=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*Junction*.fasta"`

if [ `cat $file_check |wc -l` -gt 0 ]; 
then 
cat $working_directory/REF_SEQS/BLASTING/NEW_HAPS/*Junction*.fasta >> $temp_folder/JUNCTION_REFS.fasta
fi

#############################################################################################################





               bwa index JUNCTION_REFS.fasta
               bwa mem -t $number_of_threads JUNCTION_REFS.fasta $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
               samtools view -h -F 4 $SPECIMEN_NAME.alignment.sam | samtools view -bS > $SPECIMEN_NAME.mapped_only.sam  #### take mapped reads only
               samtools view $SPECIMEN_NAME.mapped_only.sam | awk '{print("@"$1"\n"$10"\n+\n"$11)}' > $SPECIMEN_NAME.mapped_only.fastq



               # convert the fastq sequence just created to a fasta
               $working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/miraconvert $SPECIMEN_NAME.mapped_only.fastq $SPECIMEN_NAME.mapped_only.fasta



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





READS_MAPPING_TO_YOUR_JUNCTION=`ls -1 | wc -l`






                              ## now we need to find reads that have the primer between each end (or partial primer). Then we need to list the names of those reads.
                              #>junction_primer_forward
                              #TACCAAAGCATCCATCTACAGC
                              #>junction_primer_reverse-reverse_complement
                              #AACACGATCCGATTGCTTGGG --- this is the full length primer, but shorten to this: AACACGATCCGATTG so that we detect China.
                              #>reverse_sequence_of_forward_junction_primer
                              #GCTGTAGATGGATGCTTTGGTA
                              #>junction_reverse_primer
                              #CAATCGGATCGTGTT --- shortened so we detect china.


####find those sequences that mapped to the junction that contain a section of both primers

for my_fastas in *
do
#grep -l 'TACCAAAGCATCCATCTACAGC.*AACACGATCCGATTG' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
#grep -l 'CAATCGGATCGTGTT.*GCTGTAGATGGATGCTTTGGTA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt ### ORIGINAL THAT ALWAYS WORKED BEFORE


#grep -l 'TCCATCTACAGC.*AACACGATC' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
#grep -l 'GATCGTGTT.*GCTGTAGATGGA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt

###I made this really short. I wanted to maximize the changes of finding reads containing the full length repeat. To ensure sensitivity, I made sure that if
### the forward primer signal was short, then the other had to be longer. I search for primer signatures in both the forward and reverse directions.


grep -l 'CCATCTACAGC.*AACAC' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
grep -l 'GTGTT.*GCTGTAGATGG' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt

grep -l 'TCTACAGC.*AACACGATC' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
grep -l 'GATCGTGTT.*GCTGTAGA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt

done

####because you search twice above now, you might want to remove duplicate rows.

sort -u ../$SPECIMEN_NAME.list_of_sequences_to_review.txt > ../$SPECIMEN_NAME.2.list_of_sequences_to_review.txt





###actually, you will try clustering these - not assembling them.
mkdir ../assemble_these

for READS in `cat ../$SPECIMEN_NAME.2.list_of_sequences_to_review.txt`
do
cat $READS >> $SPECIMEN_NAME.assemble_these_seqs
done
cp $SPECIMEN_NAME.assemble_these_seqs ../assemble_these
cd ../assemble_these


#### Need to optimize this 95 will merge different Junction types so be warned but reduces false positives. ##################################

#the line below is the original that uses CAP3
#$working_directory/CAP3/cap3.macosx.intel64/cap3 $SPECIMEN_NAME.assemble_these_seqs.fasta -o 90 -p 99


#length1=`cat *.contigs | wc -l`
#if [ $length1 -gt 1 ];
#then cat *.contigs >> ../$SPECIMEN_NAME.map_to_this.fasta

#fi

#length2=`cat *.singlets | wc -l`
#if [ $length2 -gt 1 ] && [ $length1 == 0 ];
#then cat *.singlets >> ../$SPECIMEN_NAME.map_to_this.fasta

#fi

##all above is good code that works, but it is slow. Will the code below produce the same result?



############################ Will this suffice as it is quicker than CAP3?

awk '/^>/{print ">Contig_" ++i; next}{print}' < $SPECIMEN_NAME.assemble_these_seqs > $SPECIMEN_NAME.cluster_these_seqs.fasta


##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $SPECIMEN_NAME.cluster_these_seqs.fasta \
-o $temp_folder/assemble_these/$SPECIMEN_NAME.junction_clusters.fasta -c 1 -g 1 -d 0 -T $number_of_threads 
##########################################################################################################################################################################################################################

$working_directory/CAP3/cap3.macosx.intel64/cap3 $temp_folder/assemble_these/$SPECIMEN_NAME.junction_clusters.fasta -o 95 -p 99.99



cp  $SPECIMEN_NAME.junction_clusters.fasta.cap.contigs ../$SPECIMEN_NAME.map_to_this.fasta
cat $SPECIMEN_NAME.junction_clusters.fasta.cap.singlets >> ../$SPECIMEN_NAME.map_to_this.fasta



###### below here is as it was.
cd ..


###TIDY UP

rm -rf mapped_only *mapped_only_single_line.fasta assemble_these *list_of_sequences_to_review


cat $SPECIMEN_NAME.map_to_this.fasta | awk '/^>/{print ">new_junction_" ++i "_sequence"; next}{print}' > $SPECIMEN_NAME.map_to_this_2.fasta

seqtk seq -F '#' $SPECIMEN_NAME.map_to_this_2.fasta > $SPECIMEN_NAME.map_to_this.fastq

cp $SPECIMEN_NAME.map_to_this.fastq map_to_this.fastq
cp $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq clean_merged.fastq

$working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/mira /$working_directory/MOD_Junction_manifest_file.txt

cp $temp_folder/novel_junction_finder_assembly/novel_junction_finder_d_results/novel_junction_finder_out_new_junc.unpadded.fas* .

$working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/miraconvert novel_junction_finder_out_new_junc.unpadded.fasta $SPECIMEN_NAME.newly_found_junction.fastq

## TIDY UP
rm -rf novel_junction_finder_assembly novel_junction_finder_out* *map_to_this.fastq


## now we need to find reads that have the primer between each end (or partial primer). Then we need to list the names of those reads.
                              #>junction_primer_forward
                              #TACCAAAGCATCCATCTACAGC
                              #>junction_primer_reverse-reverse_complement
                              #AACACGATCCGATTGCTTGGG --- this is the full length primer, but shorten to this: AACACGATCCGATTG so that we detect China.


                              #>reverse_sequence_of_forward_junction_primer
                              #GCTGTAGATGGATGCTTTGGTA
                              #>junction_reverse_primer
                              #CAATCGGATCGTGTT --- shortened so we detect china.


#### CUTADAPT TO TRIM OFF THE PRIMERS -- primer lengths were shortened to less than the full primer.. sometimes when there was an error in the priming site, cutadapt wouldnt cut, so you end up with this wierd long sequence
#### that had a piece of rubbish at the end, but an actual real junction type. If you continue to catch these wierd artifacts, you could try dropping a base or two from each primer -- but not too short! you will lose specificity.

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.newly_found_junction.fastq -g ACAGC --output $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq ### five prime adapter only -- was this TACCAAAGCATCCATCTACAGC --- changed to this: CCATCTACAGC
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACAC --output $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq ### THREE prime adapter only # was this AACACGATCCGATTG but changed to this: AACACGATCC
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g GTGTT --output $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq ### five prime adapter only   ----- IF CONTIG IS REVERSED --- was this CAATCGGATCGTGTT -- changed to this: GGATCGTGTT
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT --output $SPECIMEN_NAME.REV_three_prime_trim.fasta  --fasta ### five prime adapter only   ----- IF CONTIG IS REVERSED -- was this GCTGTAGATGGATGCTTTGGTA -- changed to this: GCTGTAGATGG

#cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $SPECIMEN_NAME.novel_junction_sequence_found.fasta







cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_0_junction_sequences.fasta ###### save these sequences and add them to the clusters at the end of the hap caller script just before blasting.

### Now we need to remove contigs that are not real junction sequences (artefacts), by searcher for various signatures using "or" loops.



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


for my_fastas in *
do

if 

#my_fastas=new_junction_1_sequence_bb.fasta

   
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



####FIRST PASS TO FIND REAL JUNCTION CONTIGS
##########################################################################################################################################
# MAP READS TO THE OUTPUT OF PREVIOUS CLUSTERING WITH BOWTIE - USE HIGHLY STRINGENT PARAMETERS
#Parameters: -N = number of mismatches -L = length of seed substring -i = specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length.
# -D = Up to <int> consecutive seed extension attempts can “fail” before Bowtie 2 moves on. -R the maximum number of times Bowtie 2 will “re-seed” reads with repetitive seeds.
$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam




##### SORT BAM BEFORE MAKING CONSENSUS
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam


# call variants
bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam | bcftools call -mv -Oz -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.consensus.fasta





##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta \
-o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads 
##########################################################################################################################################################################################################################




$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam





##### SORT BAM BEFORE MAKING CONSENSUS
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam


# call variants
bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.bam | bcftools call -mv -Oz -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.vcf.gz > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001.consensus.fasta


seqtk seq -F '#' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fastq




#### CUTADAPT TO TRIM OFF THE PRIMERS -- primer lengths were shortened to less than the full primer.. sometimes when there was an error in the priming site, cutadapt wouldnt cut, so you end up with this wierd long sequence
#### that had a piece of rubbish at the end, but an actual real junction type. If you continue to catch these wierd artifacts, you could try dropping a base or two from each primer -- but not too short! you will lose specificity.

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fastq -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq ### five prime adapter only -- was this TACCAAAGCATCCATCTACAGC --- changed to this: CCATCTACAGC
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACAC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq ### THREE prime adapter only # was this AACACGATCCGATTG but changed to this: AACACGATCC
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g GTGTT --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq ### five prime adapter only   ----- IF CONTIG IS REVERSED --- was this CAATCGGATCGTGTT -- changed to this: GGATCGTGTT
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim.fasta  --fasta ### five prime adapter only   ----- IF CONTIG IS REVERSED -- was this GCTGTAGATGGATGCTTTGGTA -- changed to this: GCTGTAGATGG

#cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $SPECIMEN_NAME.novel_junction_sequence_found.fasta


/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim.fasta -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim2.fasta

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim2.fasta -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim3.fasta

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim3.fasta -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim4.fasta

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim4.fasta -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim5.fasta

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim5.fasta -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim6.fasta

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim6.fasta -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta





$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT2.bam



##### SORT BAM BEFORE MAKING CONSENSUS
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT2.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted2.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted2.bam


# call variants
bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted2.bam | bcftools call -mv -Oz -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls2.vcf.gz

bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls2.vcf.gz

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls2.vcf.gz > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00002.consensus.fasta








###get your sequences on one line
/usr/local/bin/gsed -e 's/\(^>.*$\)/#\1#/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00002.consensus.fasta \
| tr -d "\r" | tr -d "\n" | /usr/local/bin/gsed -e 's/$/#/' | tr "#" "\n" | /usr/local/bin/gsed -e '/^$/d' > \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

### now blast the consensus.
### if the original had a hit in the database, remove the consensus



for CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00002.consensus.fasta | /usr/local/bin/gsed 's/>//g'`

do

echo $CHECK_CONSENSUS_DRAFT_1 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt


# echo ">new_junction_4_sequence_bb" $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta
grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta > $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta

rm $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta | tail -1` == `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta | tail -1` ];

then /usr/local/bin/gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

fi

done




#####check for contigs that have N bases in their sequence and other junk

for CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001.consensus.fasta | /usr/local/bin/gsed 's/>//g'`

do

echo $CHECK_CONSENSUS_DRAFT_1 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt


# echo ">new_junction_4_sequence_bb" $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_CONSENSUS.fasta
grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta > $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_BEFORE_consensus.fasta

rm $temp_folder/CONFIRMED_JUNCTIONS/reference.txt


if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_CONSENSUS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];

then /usr/local/bin/gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

#then echo "The consensus has an N base"

fi

if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1_BEFORE_consensus.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];

#then echo "The original has an N base"

then /usr/local/bin/gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

fi


if [ `/usr/local/bin/gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta | tail -1 | wc -l` -eq 0 ] && \
[ `/usr/local/bin/gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];

#then echo "The original has an N base"


then /usr/local/bin/gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta

fi

if [ `/usr/local/bin/gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ] && \
[ `/usr/local/bin/gsed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];

#then echo "The original has an N base"


then /usr/local/bin/gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta
/usr/local/bin/gsed -i "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

fi

done






for SCND_CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.00001_junction_sequences.fasta | /usr/local/bin/gsed 's/>//g'`

do

echo $SCND_CHECK_CONSENSUS_DRAFT_1 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt


# echo ">new_junction_4_sequence_bb" $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.fasta
grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta > $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.fasta

#rm $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb \
-in $temp_folder/JUNCTION_REFS.fasta \
-input_type fasta \
-dbtype nucl

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



if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` == `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

fi

if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -gt 0 ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

fi


if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -eq 0 ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta

fi


if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -eq 0 ] && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT_1_CONSENSUS.blast_result | wc -l` -gt 0 ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta

fi

done


cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.00002.consensus.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.01_junction_sequences.fasta
cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_HEURISTICS.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.01_junction_sequences.fasta


rm $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.calls.*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT.bam
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta_BT_*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta.fai
#rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.aln.sorted.ba*
#rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1_junction_sequences.fasta
rm -rf $temp_folder/CONFIRMED_JUNCTIONS/artefact_remove
#rm $temp_folder/CONFIRMED_JUNCTIONS/X_DRAFT*
#rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1.consensus.fasta
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.OUTPUT2.bam
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REV_three_prime_trim*
rm $temp_folder/CONFIRMED_JUNCTIONS/Y_DRAFT_1*



##### JOEL YOU MIGHT WANT TO CLUSTER file Draft 1.1 here in later versions ---- some of your contigs are identical
###This is because there is potential for 


############# NOW LETS CLUSTER ---new addition. test on MO023

##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.01_junction_sequences.fasta \
-o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads 
##########################################################################################################################################################################################################################




seqtk seq -F '#' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fastq



####THIS NEXT PART IS NEW TEST ---- NEED TO CONVERT TO A FASTQ OR IT WILL NOT WORK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


## CUTADAPT TO TRIM OFF THE PRIMERS -- primer lengths were shortened to less than the full primer.. sometimes when there was an error in the priming site, cutadapt wouldnt cut, so you end up with this wierd long sequence
#### that had a piece of rubbish at the end, but an actual real junction type. If you continue to catch these wierd artifacts, you could try dropping a base or two from each primer -- but not too short! you will lose specificity.

/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.02_junction_sequences.fastq -g ACAGC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq ### five prime adapter only -- was this TACCAAAGCATCCATCTACAGC --- changed to this: CCATCTACAGC
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACAC --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq ### THREE prime adapter only # was this AACACGATCCGATTG but changed to this: AACACGATCC
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g GTGTT --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq ### five prime adapter only   ----- IF CONTIG IS REVERSED --- was this CAATCGGATCGTGTT -- changed to this: GGATCGTGTT
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT --output $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.1_junction_sequences.fasta  --fasta ### five prime adapter only   ----- IF CONTIG IS REVERSED -- was this GCTGTAGATGGATGCTTTGGTA -- changed to this: GCTGTAGATGG


##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.1_junction_sequences.fasta \
-o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads 
##########################################################################################################################################################################################################################





##################### SECOND PASS FOR IDENTIFYING NEXT JUNCTION ##################################################

### this is where we generate draft 2 of the junction sequences

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.OUTPUT.bam


##### EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.aln.sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.aln.sorted.bam


#cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.consensus_after_mira.fasta | /usr/local/bin/gsed '/>/!d' > $temp_folder/GOOD_HAPS_AMONG_CLUSTERS.txt

samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.aln.sorted.bam > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.coverage # exports a file containing the coverage.


for JUNK_REMOVAL in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.coverage | cut -f 1 | uniq`

do

         if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$JUNK_REMOVAL\t/$JUNK_REMOVAL\t/g" | grep "$JUNK_REMOVAL" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt 10 ]

then echo $JUNK_REMOVAL >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.contigs_with_coverage_less_than.10.txt

fi

done


###get your sequences on one line
#/usr/local/bin/gsed -e 's/\(^>.*$\)/#\1#/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.1_junction_sequences.fasta \
#| tr -d "\r" | tr -d "\n" | /usr/local/bin/gsed -e 's/$/#/' | tr "#" "\n" | /usr/local/bin/gsed -e '/^$/d' > \
#$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta



cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1.11_junction_sequences.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta


/usr/local/bin/gsed -i 's/^/>/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.contigs_with_coverage_less_than.10.txt


## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST OF NEW ONES THAT POSSESS LOW COVERAGE
for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.contigs_with_coverage_less_than.10.txt`
do
/usr/local/bin/gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta #removes the line containing the fasta sequence name. The '+1' removes the line after as well.
done


rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_1*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_0*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CUTADAPT* 
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FILTER_ON_H*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.ONE_LINE.DRAFT_1*


#####DRAFT 2 RELAXED



$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_consensus.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_consensus.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED.OUTPUT.bam



##### SORT BAM BEFORE MAKING CONSENSUS
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_sorted.bam


# call variants
bcftools mpileup -Ou -f $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_sorted.bam | bcftools call -mv -Oz -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.calls3.vcf.gz

bcftools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.calls3.vcf.gz

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta | bcftools consensus $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.calls3.vcf.gz > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.RELAXED.consensus.fasta


### NOW WE DO SOME BLASTING TO CHECK FOR THE BUGS.





###get your sequences on one line
/usr/local/bin/gsed -e 's/\(^>.*$\)/#\1#/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.RELAXED.consensus.fasta \
| tr -d "\r" | tr -d "\n" | /usr/local/bin/gsed -e 's/$/#/' | tr "#" "\n" | /usr/local/bin/gsed -e '/^$/d' > \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta






for SCND_CHECK_CONSENSUS_DRAFT_2 in `grep -h '>' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta | /usr/local/bin/gsed 's/>//g'`

do

echo $SCND_CHECK_CONSENSUS_DRAFT_2 > $temp_folder/CONFIRMED_JUNCTIONS/reference.txt


# echo ">new_junction_4_sequence_bb" $temp_folder/CONFIRMED_JUNCTIONS/reference.txt

grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.fasta
grep -A 1 -wFf $temp_folder/CONFIRMED_JUNCTIONS/reference.txt $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.fasta

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

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta

fi

if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -gt 0 ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta

fi


if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -eq 0 ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta

fi


if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -eq 0 ] && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -gt 0 ];

then /usr/local/bin/gsed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta

fi

done







cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX.fasta
cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.consensus.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX.fasta

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta








################################################## THIRD PASS FOR IDENTIFYING THE NEXT JUNCTION

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.OUTPUT.bam


##### EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.aln.sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.aln.sorted.bam


samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.aln.sorted.bam > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.coverage # exports a file containing the coverage.






JUNCTION_rough_coverage_cutoff=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.010)/1" | bc` ################## test nepal make sure it only has 2 junctions.



for JUNK_REMOVAL in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.coverage | cut -f 1 | uniq`

do

         if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$JUNK_REMOVAL\t/$JUNK_REMOVAL\t/g" | grep "$JUNK_REMOVAL" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff ]

then echo $JUNK_REMOVAL >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt

fi

done


/usr/local/bin/gsed -i 's/^/>/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2.200.consensus.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.consensus.fasta

## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST OF NEW ONES THAT POSSESS LOW COVERAGE
for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt`
do
/usr/local/bin/gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.consensus.fasta #removes the line containing the fasta sequence name. The '+1' removes the line after as well.
done


rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_2*

############# NOW LETS CLUSTER

##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.consensus.fasta \
-o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta -c 1 -g 1 -d 0 -T $number_of_threads 
##########################################################################################################################################################################################################################



rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.AFTER_RELAX*



$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.OUTPUT.bam


##### EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.bam


samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.bam > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage # exports a file containing the coverage.



JUNCTION_rough_coverage_cutoff_2=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.015)/1" | bc` ################## DO NOT CHANGE THIS CUTOFF -- TRUST ME. JUNCTION I A DIFFICULT LOCUS TO DEAL WITH AND THIS CUTOFF IS OPTIMAL FOR THE DATA.



### find out which haplotypes contain bases with coverage less than X reads
for NEW_HAPS in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | cut -f 1 | uniq` 
      do
          if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $REQUIRED_DEPTH_NEW_HAP ] \
          && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff_2 ];


               then echo $NEW_HAPS >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.$REQUIRED_DEPTH_NEW_HAP.txt
          fi
      done


/usr/local/bin/gsed -i 's/^/>/' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.$REQUIRED_DEPTH_NEW_HAP.txt




## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST OF NEW ONES THAT POSSESS LOW COVERAGE
for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.$REQUIRED_DEPTH_NEW_HAP.txt`
do
/usr/local/bin/gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta #removes the line containing the fasta sequence name. The '+1' removes the line after as well.
done

cp $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3.CLUSTERS_consensus.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.validated_junction_sequences.fasta 


rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_3*



###### NOW WE WILL LOOK FOR ANY SHARP DROPS AND RISES IN THE COVERAGE OVER THE CONTIG (in DRADT_4). IF THESE EXIST, IT IS MORE THAN LIKELY A FALSE CONTIG


for JUNCTION_TESTER in `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | cut -f 1 | uniq` 
      do

#candidate_junctions=`cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | cut -f 1 | uniq`  ## list of the contigs that could be real.


coverage_at_each_base=`cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage | /usr/local/bin/gsed -n "/$JUNCTION_TESTER/p" | cut -f 3`
number_of_bases=`echo $coverage_at_each_base | wc -w`
sum_of_bases=`echo $coverage_at_each_base | xargs | sed -e 's/\ /+/g' | bc`
average_coverage=`echo "($sum_of_bases/$number_of_bases)" | bc`

standardDeviation=$(
    echo "$coverage_at_each_base" |
        awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}'
)

two_SDs=`echo "(2*$standardDeviation)" | bc`

### now if the average coverage is less than 2 standard deviations from the average coverage, then delete the contig.
### It means that there is a sharp drop in coverage somehwere, indicative of some erroneous bases.

rounded_2SDs=`echo ${two_SDs%%.*}`

if [ $rounded_2SDs -gt $average_coverage ];

##delete the contig in question
then /usr/local/bin/gsed -i "/$JUNCTION_TESTER/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.validated_junction_sequences.fasta

fi

done

## THE ABOVE WILL NOT WORK IF THE VARIABLES ARE DECIMALS


#echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.010)/1" | bc




####THESE LINES WHERE HERE BEFORE

rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.OUTPUT.bam
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.aln.sorted.ba*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.coverage
rm -rf $temp_folder/CONFIRMED_JUNCTIONS/artefact_remove



cd $temp_folder





##### you are now done confirming the sequence of the mt junctions that you have found

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4.validated_junction_sequences.fasta > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.FINAL.validated_junction_sequences.fasta


##########################################################################################################################################

#This will split the newly discovered junction haplotypes
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

rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.DRAFT_4*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.contigs_with_coverage*
rm $temp_folder/CONFIRMED_JUNCTIONS/Z_DRAFT_2*


for FOOBAR in `ls -1`
do 


#FOOBAR=`ls -1`

junction_count=`grep 'Junction_Hap' $temp_folder/JUNCTION_REFS.fasta | wc -l`
new_number=`echo "$junction_count + 1" | bc`

cat $FOOBAR | awk '/^>/{print ">Junction_Hap_'$new_number'"; next}{print}' > $SPECIMEN_NAME.junction_recount_$FOOBAR

# NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.
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
rep_length_2=`echo "$repeat_length + 21 + 22" | bc`


/usr/local/bin/gsed -i "1s/.*/>Mt_Cmt$rep_length_2.X_Junction_Hap_$new_number/" $SPECIMEN_NAME.junction_recount_$FOOBAR

cat $SPECIMEN_NAME.junction_recount_$FOOBAR >> $temp_folder/JUNCTION_REFS.fasta
cat $SPECIMEN_NAME.junction_recount_$FOOBAR >> $working_directory/REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.Mt_Cmt$rep_length_2.X_Junction_Hap_$new_number.fasta

elif [ $match_present -gt 0 ];
then echo "no new junction sequences were found"
fi

done

cd ..

rm -rf new_potential_junctions
rm -rf *fasta.amb *fasta.ann *fasta.bwt *fasta.nhr *fasta.nin *fasta.nsq *fasta.pac *fasta.sa *junction_recount_* *review.txt *RESULT* *novel* *mapped_only* *map_to_this*


done

rm $working_directory/MOD_Junction_manifest_file.txt

cd $temp_folder



