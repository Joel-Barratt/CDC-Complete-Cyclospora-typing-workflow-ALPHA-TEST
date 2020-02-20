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


grep -l 'TCCATCTACAGC.*AACACGATC' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
grep -l 'GATCGTGTT.*GCTGTAGATGGA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt



done


###actually, you will try clustering these - not assembling them.
mkdir ../assemble_these

for READS in `cat ../$SPECIMEN_NAME.list_of_sequences_to_review.txt`
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



#### CUTADAPT TO TRIM OFF THE PRIMERS
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.newly_found_junction.fastq -g TACCAAAGCATCCATCTACAGC --output $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq ### five prime adapter only
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_five_prime_trim.fastq -a AACACGATCCGATTG --output $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq ### THREE prime adapter only
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_three_prime_trim.fastq -g CAATCGGATCGTGTT --output $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq ### five prime adapter only   ----- IF CONTIG IS REVERSED
/Users/joelbarratt/Library/Python/3.7/bin/cutadapt $SPECIMEN_NAME.CUTADAPT_REV_five_prime_trim.fastq -a GCTGTAGATGGATGCTTTGGTA --output $SPECIMEN_NAME.REV_three_prime_trim.fasta  --fasta ### five prime adapter only   ----- IF CONTIG IS REVERSED

cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $SPECIMEN_NAME.novel_junction_sequence_found.fasta



cat $SPECIMEN_NAME.REV_three_prime_trim.fasta >> $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.validated_junction_sequences.fasta ###### save these sequences and add them to the clusters at the end of the hap caller script just before blasting.

####NOTE THAT CHINA IS NOT CONSIDERED IN THE CUTADAPT STEP YET. CONSIDER ITS ADDITION. BY CUTTING THE PRIMERS IN HALF TO THE POINT WHERE THE CHINESE PRIMER ENDS.



rm -rf $SPECIMEN_NAME.CUTADAPT_* $SPECIMEN_NAME.newly_found_junction.fastq $SPECIMEN_NAME.REV_three_prime_trim.fasta clean_merged*


#THE CODE BELOW IS EXPERIMENTAL.




##########################################################################################################################################################################################################################
# WE CLUSTER THE CLUSTERS PREVIOUSLY FOUND. THIS IS BECAUSE WE EXTRACTED REDUNDANT SEQUENCES IN THE PREVIOUS STEP 
$working_directory/CD-HIT/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.validated_junction_sequences.fasta \
-o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta -c 1 -g 1 -d 0 -T $number_of_threads
##########################################################################################################################################################################################################################



##########################################################################################################################################
# MAP READS TO THE OUTPUT OF PREVIOUS CLUSTERING WITH BOWTIE - USE HIGHLY STRINGENT PARAMETERS
#Parameters: -N = number of mismatches -L = length of seed substring -i = specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length.
# -D = Up to <int> consecutive seed extension attempts can “fail” before Bowtie 2 moves on. -R the maximum number of times Bowtie 2 will “re-seed” reads with repetitive seeds.
$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta \
$temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta_BT_INDEX

$working_directory/BOWTIE/bowtie2-2.3.5.1-macos-x86_64/bowtie2 -x $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta_BT_INDEX \
-U $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads $number_of_threads --local > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.OUTPUT.bam



##### EXTRACT COVERAGE INFORMATION AFTER BOWTIE MAPPING SO THAT HAPS WITH LOW COVERAGE CAN BE THROWN AWAY.
samtools sort -T -n $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.OUTPUT.bam -o $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.aln.sorted.bam
samtools index $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.aln.sorted.bam

cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta | /usr/local/bin/gsed '/>/!d' > $temp_folder/GOOD_HAPS_AMONG_CLUSTERS.txt
samtools depth -aa -d 0 $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.aln.sorted.bam > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.coverage # exports a file containing the coverage.




## tidy up un-needed files so you dont run out of space on your hard drive.

rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.aln.sorted.bam
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.aln.sorted.bam.bai
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.OUTPUT.bam
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta_BT_INDEX*
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta.clstr
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.coverage
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.possible_new_junctions.fasta
rm $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta


##########################################################################################################################################################################################################################



#$READS_MAPPING_TO_YOUR_JUNCTION
#$REQUIRED_DEPTH_NEW_HAP    this is the one that is set at the very begining of the script.


JUNCTION_rough_coverage_cutoff=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.20)/1" | bc`


#cd $temp_folder/CONFIRMED_JUNCTIONS/

### find out which haplotypes contain bases with coverage less than X reads
for NEW_HAPS in `cat GOOD_HAPS_AMONG_CLUSTERS.txt` 
      do
          if [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $REQUIRED_DEPTH_NEW_HAP ] \
          && [ `cat $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.coverage | /usr/local/bin/gsed 's/^/>/' | /usr/local/bin/gsed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff ];


               then echo $NEW_HAPS >> $temp_folder/CONFIRMED_JUNCTIONS/clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt
          fi
      done


#####  THIS IS WHERE YOU SET THE REQUIRED DEPTH ---- SEE ABOVE


##########################################################################################################################################################################################################################
cp $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.CLUSTERED_validated_junction_sequences.fasta $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.possible_new_junctions.fasta

## this will put underscores at the begining and end of fasta names -- if they are kept solely numerical, this can create problems
/usr/local/bin/gsed -i 's/>/>_/g' $temp_folder/CONFIRMED_JUNCTIONS/clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt
/usr/local/bin/gsed -i '/^>_/ s/$/_/' $temp_folder/CONFIRMED_JUNCTIONS/clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt
/usr/local/bin/gsed 's/>/>_/g' $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.possible_new_junctions.fasta | /usr/local/bin/gsed '/^>_/ s/$/_/' > $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta
##########################################################################################################################################################################################################################

## REMOVES HAPLOTYPES FROM YOUR POTENTIAL LIST OF NEW ONES THAT POSSESS LOW COVERAGE
for HAPS_TO_REMOVE in `cat $temp_folder/CONFIRMED_JUNCTIONS/clusters_with_coverage_less_than.$REQUIRED_DEPTH_NEW_HAP.txt`
do
/usr/local/bin/gsed -i "/$HAPS_TO_REMOVE/,+1 d" $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta #removes the line containing the fasta sequence name. The '+1' removes the line after as well.
done



#THIS IS WHAT YOU ORIGINALLY HAD


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
done < $temp_folder/CONFIRMED_JUNCTIONS/$SPECIMEN_NAME.REMOVE_BAD_SEQUENCES.fasta


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



