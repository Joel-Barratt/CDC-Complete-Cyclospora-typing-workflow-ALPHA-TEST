#!/bin/bash
#####################################################################################################################################################################
#### IMPORT DATA AND TRIM RAW READS MODULE ##########################################################################################################################
#####################################################################################################################################################################

cd /Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER_CYCLO_V2
rm -rf TMP
#####################################################################################################################################################################
######                                                                         ######################################################################################
######  SPECIFY LOCATION OF REFERENCE SEQUENCES AND RAW READS FOR ANALYSIS     ######################################################################################
######                                                                         ######################################################################################
#####################################################################################################################################################################
HAPLOTYPE_CALLER=/HAPLOTYPE_CALLER_CYCLO_V2    # do not modify this line                   ##########################################################################
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

#### IF YOU FIND A NEW HAPLOTYPE, WHAT IS THE MINIMUM DEPTH REQUIRED FOR IT TO BE CONSIDERED REAL? THIS IS SPECIFICALLY FOR NEW HAPLOTYPE DISCOVERY.

##### PLEASE NOTE THAT THERE IS A DYNAMIC DEPTH DETERMINATION FEATURE OF THIS SOFTWARE FOR DETECTING NEW HAPLOTYPES, AND THIS WILL SLIDE DEPENDING ON THE COVERAGE OBTAINED FOR A PARTICULAR AMPLICON
##### The setting below only comes into affect if you have very low coverage for a particular amplicon. For example, if you have 20 fold coverage of your marker, would you 
##### feel comfortable calling two new haplotypes from this (one with coverage 12 and another with coverage 8 for example)? I don't, so I set this threshold to 100 per new haplotype as a default.
##### Alternatively, if the locus has 500 fold coverage, new haplotypes are called only if their coverage is equal to, or rises above 25% of the total population of reads that map to this locus;
##### that is 125 reads in this case. If the coverage of a locus is 1000, new haplotypes must have coverage of 250 reads at least. Note that this is based on average coverage across the locus.

REQUIRED_DEPTH_NEW_HAP=100  #### I set this to 100 as a default -- remember this only comes into affect when the coverage of a locus is very low. Otherwise a sliding cutoff is used.

###########



### THIS IS NOT FOR HAPLOTYPE DISCOVERY -- THIS IS THE DEPTH REQUIRED TO DETERMINE IF A SPECIMEN POSSESSES A HAPLOTYPE. THIS IS WHY THE DEPTH IS LESS.
#### THIS IS A LESS STRINGENT CUTOFF
##### PLEASE NOTE THAT THERE IS A DYNAMIC DEPTH DETERMINATION FEATURE OF THIS SOFTWARE FOR DETECTING NEW HAPLOTYPES, AND THIS WILL SLIDE DEPENDING ON THE COVERAGE OBTAINED FOR A PARTICULAR AMPLICON
##### The setting below only comes into affect if you have very low coverage for a particular amplicon. For example, if you have 20 fold coverage of your marker, would you 
##### feel comfortable calling two haplotypes from this (one with coverage 12 and another with coverage 8 for example)? I don't, so I set this threshold to 20 to determine if a haplotype exists in a specimen.
##### Alternatively, if the locus has 500 fold coverage, new haplotypes are called only if their coverage is equal to, or rises above 10% of the total population of reads that map to this locus;
##### that is 50 reads in this case. If the coverage of a locus is 1000, new haplotypes must have coverage of 250 reads at least. Note that this is based on average coverage across the locus.

MINIMUM_DEPTH_TO_ASSIGN_HAPLPLOTYPES_TO_SPECIMENS=20  #### I set this to 20 as a default -- remember this only comes into affect when the coverage of a locus is very low. Otherwise a sliding cutoff is used.



# ONLY MODIFY TEXT TO THE RIGHT OF THE 'EQUALS' SIGN "=" BELOW - DO NOT MODIFY TEXT TO THE LEFT THE 'EQUAL' SIGN
# IN WHAT DIRECTORY DID YOU PLACE THE "HAPLOTYPE_CALLER" FOLDER?
working_directory=/Users/joelbarratt/Documents/CYCLOSPORA

# WHAT IS THE LOCATION OF YOUR ILLUMINA ADAPTER SEQUENCES FOR ADAPTER TRIMMING?
# Modify this location if required, but I have already provided a list of defauly Illumina adapters that should work for most workflows.

illumina_adapaters=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/TRIMMING/Illumina_adapters.fasta                                           
                                                                                                                                               
#### PROVIDE THE LOCATION OF A FASTA FILE THAT CONTAINS YOUR REFERENCES ****BEFORE**** YOU SPLIT THE HAPLOTYPES
## THESE SEQUENCES WILL BE USED TO RECOVER READS THAT MAP TO YOUR REFERENCES

#### NON JUNCTION BELOW
read_recovery_refs=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/READ_RECOVERY/MAPPING_NON_JUNCTION_REFERENCES_WITH_PRIMER_FEB_2020.fasta


#### JUNCTION BELOW
read_recovery_refs_JUNCTION=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/READ_RECOVERY/MAPPING_JUNCTION_WITH_PRIMERS_FEB_2020.fasta




# WHAT IS THE LOCATION OF YOUR REFERENCE SEQUENCES (YOUR KNOWN CYCLOSPORA HAPLOTYPES -- THESE ARE THE ONES WHERE YOU IDENTIFY HAPLOTYPES AFTER BLASTING)?
# The location below should be fine for the automatic Cyclospora workflow - you probably shouldn't modify this.





#FOR CYCLOSPORA THIS MUST INCLUDE ALL HAPLOTYPES AFTER SEGMENTING, AS WELL AS ALL JUNCTION TYPES WITHOUT PRIMER SEQUENCE.
complete_reference_database=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/BLASTING/ORIGINAL_REFS/CYCLOSPORA_REFERENCES_SHORTENED_FEB_2020.fasta


# WHAT IS THE LOCATION OF THE RAW ILLUMINA READS THAT YOU WANT ANALYZED?                                                                                   
# R1 and R2 reads must be placed in the folder that you specify here. They will not be moved or modified.

input_reads=$working_directory$HAPLOTYPE_CALLER/INPUT_READS



###### ADD LOCATION OF THE BED FILE AS A VARIABLE HERE --- IF RUNNING CYCLOSPORA DO NOT INCLUDE THE JUNCTION SEQUNECES HERE.
bed_references=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/BED_FILE/CYCLOSPORA_FEB_11_2020.bed

####IF YOU MADE THIS BED FILE USING XLSX YOU WILL NEED TO RUN DOS2UNIX ON THE FILE TO CONVERT IT OR THE SCRIPT WONT RUN CORRECTLY!

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
echo $read_recovery_refs > READ_REC
echo $MINIMUM_DEPTH_TO_ASSIGN_HAPLPLOTYPES_TO_SPECIMENS > DEPTH_FOR_CALLING
echo $bed_references > $LOC/TMP/BED_LOCATION
echo $read_recovery_refs_JUNCTION > $LOC/TMP/READ_REC_JUNCTION

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
                     mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e 's/\(.\{11\}\).*.\(.\{3\}\)/\1.R1.fastq\2/'`  #### NOTE THAT THIS WILL TRUNCATE THE FILE NAME TO 11 CHARACTERS -- IS THIS OK?
                done

##MAKE LIST OF SPECIMENS TO GENOTYPE
   for file in *
                 do
                      echo "$file" | sed -e 's/-/_/g' | sed -e 's/^\(.\{11\}\).*/\1/' >> $LOC/TMP/SPECIMENS_TO_TRIM   # MAKES A LIST OF SPECIMENS THAT NEED TRIMMING -  NOTE THAT THIS WILL TRUNCATE THE FILE NAME TO 11 CHARACTERS -- IS THIS OK?
                 done

cd $LOC/TMP/R2_FILES
##RENAME THE FASTQ FILES IN R2 FOLDER
   for file in *
                  do
                       mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e 's/\(.\{11\}\).*.\(.\{3\}\)/\1.R2.fastq\2/'` #  NOTE THAT THIS WILL TRUNCATE THE FILE NAME TO 11 CHARACTERS -- IS THIS OK?

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

rm -rf R1_FILES R2_FILES

bash ../Finding_New_Haps_SNP_Based_Module.sh

bash ../Finding_New_Junction_Types_Module.sh

bash ../Call_Haplotypes.sh

cd $LOC/

Rscript START_haplotype_sheet_generator.R

