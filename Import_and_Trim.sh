#!/bin/bash
#####################################################################################################################################################################
#### IMPORT DATA AND TRIM RAW READS MODULE ##########################################################################################################################
#####################################################################################################################################################################

cd /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER
rm -rf TMP
#####################################################################################################################################################################
######                                                                         ######################################################################################
######  SPECIFY LOCATION OF REFERENCE SEQUENCES AND RAW READS FOR ANALYSIS     ######################################################################################
######                                                                         ######################################################################################
#####################################################################################################################################################################
HAPLOTYPE_CALLER=/HAPLOTYPE_CALLER  # do not modify this line                   ######################################################################################
#####################################################################################################################################################################   
																


                                            ##################################################### THE USER SHOULD MODIFY THESE LINES OF THE SCRIPT ONLY. PROVIDE THE 
                                            ################################################### PIPELINE WITH THE LOCATION OF YOUR RAW READS AND REFERENCE SEQUENCES
                                            ###
                                            ###
                                            ###
                                         #########
                                          #######
                                           #####
                                            ###
                                             #

#### IF YOU FIND A NEW HAPLOTYPE, WHAT IS THE DEPTH REQUIRED FOR IT TO BE CONSIDERED REAL?
REQUIRED_DEPTH_NEW_HAP=500


#### WHAT IS THE DEPTH REQUIRED FOR DETECTION OF A KNOWN HAPLOTYPE IN A SPECIMEN?
CALL_DEPTH=50



# ONLY MODIFY TEXT TO THE RIGHT OF THE 'EQUALS' SIGN "=" BELOW - DO NOT MODIFY TEXT TO THE LEFT THE 'EQUAL' SIGN

# IN WHAT DIRECTORY DID YOU PLACE THE "HAPLOTYPE_CALLER" FOLDER?

working_directory=/Users/joelbarratt/Documents/CYCLOSPORA



# WHAT IS THE LOCATION OF YOUR ILLUMINA ADAPTER SEQUENCES FOR ADAPTER TRIMMING?
# Modify this location if required, but I have already provided a list of defauly Illumina adapters that should work for most workflows.

illumina_adapaters=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/TRIMMING/Illumina_adapters.fasta                                           
                                                                                                                                               

                                                                                                                                             
# WHAT IS THE LOCATION OF THE RAW ILLUMINA READS THAT YOU WANT ANALYZED?                                                                                   
# R1 and R2 reads must be placed in the folder that you specify here. They will not be moved or modified.

input_reads=$working_directory$HAPLOTYPE_CALLER/INPUT_READS
                                                                  


# WHAT IS THE LOCATION OF YOUR REFERENCE SEQUENCES (YOUR KNOWN CYCLOSPORA HAPLOTYPES)?
# The location below should be fine for the automatic Cyclospora workflow - you probably shouldn't modify this.
                            
complete_reference_database=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/BLASTING/ORIGINAL_REFS/2019_ORIGINAL_REFERENCES_SHORTENED.fasta


######################################################################################################################################################################## DO NOT MODIFY SCRIPT BELOW THIS POINT.

###########################################################################################################################################################
######                                                  ###################################################################################################
######  WRITE DIRECTORIES FOR TMP FILES & VARIABLES     ###################################################################################################
######                                                  ###################################################################################################
###########################################################################################################################################################

LOC=$working_directory$HAPLOTYPE_CALLER
cd $LOC

mkdir TMP
cd TMP

echo $LOC > RUN_DIR
echo $complete_reference_database > ALL_REF
echo $LOC/TMP > TMP_FOL
echo $REQUIRED_DEPTH_NEW_HAP > DEPTH_NEW

mkdir R1_FILES
mkdir R2_FILES
mkdir $LOC/RESULTS

##########################################################################################################################################################
######                                       #############################################################################################################
######  MOVE AND RENAME THE R1 AND R2 READS  #############################################################################################################
######                                       #############################################################################################################
##########################################################################################################################################################
cd $input_reads

cp *_R1_* $LOC/TMP/R1_FILES
cp *_R2_* $LOC/TMP/R2_FILES


cd $LOC/TMP/R1_FILES
##RENAME THE FASTQ FILES IN R1 FOLDER
   for file in *
                do
                     mv "$file" `echo "$file" | sed -e 's/\(.\{11\}\).*.\(.\{3\}\)/\1.R1.fastq\2/'`
                done

##MAKE LIST OF SPECIMENS TO GENOTYPE
   for file in *
                 do
                      echo "$file" | sed -e 's/^\(.\{11\}\).*/\1/' >> $LOC/TMP/SPECIMENS_TO_TRIM   # MAKES A LIST OF SPECIMENS THAT NEED TRIMMING
                 done

cd $LOC/TMP/R2_FILES
##RENAME THE FASTQ FILES IN R2 FOLDER
   for file in *
                  do
                       mv "$file" `echo "$file" | sed -e 's/\(.\{11\}\).*.\(.\{3\}\)/\1.R2.fastq\2/'`

                  done

cd $LOC/TMP
mkdir PROCESSED_READS
##########################################################################################################################################################
######                               #####################################################################################################################
######    TRIM THE R1 AND R2 READS   #####################################################################################################################
######                               #####################################################################################################################
##########################################################################################################################################################



         for SPECIMEN_NAME in `cat SPECIMENS_TO_TRIM`
                                                            do

# TRIM USING BBDUK

bash $LOC/BBMAP/bbmap/bbduk.sh in1=$LOC/TMP/R1_FILES/$SPECIMEN_NAME.R1.fastq.gz in2=$LOC/TMP/R2_FILES/$SPECIMEN_NAME.R2.fastq.gz out1=$LOC/TMP/$SPECIMEN_NAME.clean1.fq out2=$LOC/TMP/$SPECIMEN_NAME.clean2.fq minlen=50 ref=$illumina_adapaters ktrim=r k=23 mink=11 hdist=1

# MERGE USING BBMERGE

bash $LOC/BBMAP/bbmap/bbmerge-auto.sh in1=$LOC/TMP/$SPECIMEN_NAME.clean1.fq in2=$LOC/TMP/$SPECIMEN_NAME.clean2.fq out=$LOC/TMP/$SPECIMEN_NAME.merged.fq outu=$LOC/TMP/$SPECIMEN_NAME.unmerged.fq rem k=62 extend2=50 ecct


 # NOTE FROM BBTOOLS MANUAL ADVISING AGAINST QULITY TRIMMING PRIOR TO MERGING READS:
 # Adapter-trimming reads is fine and is recommended prior to running BBDuk, particularly if kmers will be used. ***Quality-trimming is 
 # usually not recommended*****, unless the library has major quality issues at the right end of reads resulting in a low merge rate, in which 
 # case only weak trimming (such as qtrim=r trimq=8) should be used.




                                                # CONCATENATE MERGED AND UNMERGED READS AND PLACE IN A 'PROCESSED_READS' DIRECTORY FOR DOWNSTREAM USE

                                                
                                                mv $LOC/TMP/$SPECIMEN_NAME.merged.fq $LOC/TMP/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq
                                                cat $LOC/TMP/$SPECIMEN_NAME.unmerged.fq >> $LOC/TMP/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq
                                                rm -rf $SPECIMEN_NAME.clean1.fq $SPECIMEN_NAME.clean2.fq $SPECIMEN_NAME.unmerged.fq
 

                                                           done


cd $LOC/TMP/

bash ../Finding_New_Haps_SNP_Based_Module.sh
bash ../Finding_New_Junction_Types_Module.sh
bash ../Call_Haplotypes.sh



