#!/bin/bash
#####################################################################################################################################################################
####     DISCOVER NEW MT JUNCTION TYPES    ##########################################################################################################################
#####################################################################################################################################################################

#cd TMP

#rm SPECIMENS_TO_SEARCH_JUNCTION

working_directory=`cat RUN_DIR` # /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER
complete_reference_database=`cat ALL_REF` # /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER/REF_SEQS/BLASTING/ORIGINAL_REFS/2019_ORIGINAL_REFERENCES_SHORTENED.fasta
temp_folder=`cat TMP_FOL` #/Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER/TMP



# Mt Junction sequences discovered prior to 2019 are found in this location.
# I have left the primer sequences on these references because the purpose of this file is to maximize read recovery -- reads will contain primer sequence
junction_with_primer=$working_directory/REF_SEQS/READ_RECOVERY/FOR_READ_RECOVERY_Junction_sequences_with_primers_2019.fasta



############################################################################################################## get all known junction sequences for read recovery
cat $junction_with_primer > $temp_folder/JUNCTION_REFS.fasta # pull only junction sequences from original reference file.

file_check=`echo "$working_directory/REF_SEQS/BLASTING/NEW_HAPS/*Junction*.fasta"`

if [ -f $file_check ]; 
then 
cat $working_directory/CYCLO_REF_SEQS/BLASTING/NEW_HAPS/*Junction*.fasta >> $temp_folder/JUNCTION_REFS.fasta

#else echo "NO NEW HAPLOTYPES FOUND PREVIOUSLY -- LETS CHECK THE LATEST SPECIMEN"
fi

#############################################################################################################




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


#SPECIMEN_NAME=`cat SPECIMENS_TO_SEARCH_JUNCTION`


               bwa index JUNCTION_REFS.fasta
               bwa mem JUNCTION_REFS.fasta $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq > $SPECIMEN_NAME.alignment.sam
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
grep -l 'TACCAAAGCATCCATCTACAGC.*AACACGATCCGATTG' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt
grep -l 'CAATCGGATCGTGTT.*GCTGTAGATGGATGCTTTGGTA' $my_fastas >> ../$SPECIMEN_NAME.list_of_sequences_to_review.txt

done



mkdir ../assemble_these

for READS in `cat ../$SPECIMEN_NAME.list_of_sequences_to_review.txt`
do
cat $READS >> $SPECIMEN_NAME.assemble_these_seqs
done
cp $SPECIMEN_NAME.assemble_these_seqs ../assemble_these
cd ../assemble_these

$working_directory/CAP3/cap3.macosx.intel64/cap3 $SPECIMEN_NAME.assemble_these_seqs -o 90 -p 95


length1=`cat *.contigs | wc -l`
if [ $length1 -gt 1 ];
then cat *.contigs >> ../$SPECIMEN_NAME.map_to_this.fasta

fi

length2=`cat *.singlets | wc -l`
if [ $length2 -gt 1 ] && [ $length1 == 0 ];
then cat *.singlets >> ../$SPECIMEN_NAME.map_to_this.fasta

fi

cd ..


###TIDY UP

rm -rf mapped_only *mapped_only_single_line.fasta assemble_these *list_of_sequences_to_review


cat $SPECIMEN_NAME.map_to_this.fasta | awk '/^>/{print ">new_junction_" ++i "_sequence"; next}{print}' > $SPECIMEN_NAME.map_to_this_2.fasta

seqtk seq -F '#' $SPECIMEN_NAME.map_to_this_2.fasta > $SPECIMEN_NAME.map_to_this.fastq

cp $SPECIMEN_NAME.map_to_this.fastq map_to_this.fastq
cp $temp_folder/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq clean_merged.fastq

$working_directory/MIRA/mira_4.0.2_darwin13.1.0_x86_64_static/bin/mira /$working_directory/Junction_manifest_file.txt

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
####NOTE THAT CHINA IS NOT CONSIDERED IN THE CUTADAPT STEP YET. CONSIDER ITS ADDITION. BY CUTTING THE PRIMERS IN HALF TO THE POINT WHERE THE CHINESE PRIMER ENDS.



rm -rf $SPECIMEN_NAME.CUTADAPT_* $SPECIMEN_NAME.newly_found_junction.fastq $SPECIMEN_NAME.REV_three_prime_trim.fasta clean_merged*






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
done < ../$SPECIMEN_NAME.novel_junction_sequence_found.fasta








for FOOBAR in `ls -1`
do 


#FOOBAR=`ls -1`

junction_count=`grep 'Junction_Hap' $temp_folder/JUNCTION_REFS.fasta | wc -l`
new_number=`echo "$junction_count + 1" | bc`

cat $FOOBAR | awk '/^>/{print ">Junction_Hap_'$new_number'"; next}{print}' > $SPECIMEN_NAME.junction_recount_$FOOBAR

# NOW IF THAT SEQUENCE IS BLASTED AND IT IS NOT 100 PERCENT IDENTICAL WITH AND NOT 100% COVERAGE WITH ANYTHING IN THE DATABASE, THEN WE KEEP IT.
$working_directory/BLAST/ncbi-blast-2.9.0+/bin/makeblastdb -in $temp_folder/JUNCTION_REFS.fasta \
-input_type fasta \
-dbtype nucl \
-title junction_check


$working_directory/BLAST/ncbi-blast-2.9.0+/bin/blastn -db ../JUNCTION_REFS.fasta \
-query $SPECIMEN_NAME.junction_recount_$FOOBAR \
-word_size 7 \
-evalue 0.001 \
-perc_identity 100 \
-qcov_hsp_perc 100 \
-out $temp_folder/new_potential_junctions/$SPECIMEN_NAME.RESULT_$FOOBAR.blast_result \
-max_target_seqs 1 \
-outfmt "6 qseqid pident mismatch gapopen gaps sseqid sseq evalue bitscore"

match_present=`cat $temp_folder/new_potential_junctions/$SPECIMEN_NAME.RESULT_$FOOBAR.blast_result | wc -l`

if [ $match_present == 0 ];
then echo "you have discovered a new junction type"

repeat_length=`cat $SPECIMEN_NAME.junction_recount_$FOOBAR | sed "1d" | tr -cd '[:alpha:]' | wc -m`
rep_length_2=`echo "$repeat_length + 21 + 22" | bc`


/usr/local/bin/gsed -i "1s/.*/>Cmt$rep_length_2.X_Junction_Hap_$new_number/" $SPECIMEN_NAME.junction_recount_$FOOBAR

cat $SPECIMEN_NAME.junction_recount_$FOOBAR >> $temp_folder/JUNCTION_REFS.fasta
cat $SPECIMEN_NAME.junction_recount_$FOOBAR >> $working_directory/REF_SEQS/BLASTING/NEW_HAPS/$SPECIMEN_NAME.Cmt$rep_length_2.X_Junction_Hap_$new_number.fasta

elif [ $match_present -gt 0 ];
then echo "no new junction sequences were found"
fi

done

cd ..
rm -rf new_potential_junctions


rm -rf *fasta.amb *fasta.ann *fasta.bwt *fasta.nhr *fasta.nin *fasta.nsq *fasta.pac *fasta.sa *junction_recount_* *review.txt *RESULT* *novel* *mapped_only* *map_to_this*
done

cd $temp_folder



