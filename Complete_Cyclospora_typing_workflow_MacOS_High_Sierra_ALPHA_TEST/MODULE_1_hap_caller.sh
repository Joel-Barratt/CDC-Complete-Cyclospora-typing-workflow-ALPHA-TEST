#!/bin/bash

#UPDATED August 5, 2021

package=`echo "MODULE_1_hap_caller.sh"`
alpha_version=`echo "Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_ALPHA_TEST"`
HAPLOTYPE_CALLER=`echo "HAPLOTYPE_CALLER_CYCLO_V2"`


															
while test $# -lt 1; do

echo "


 ERROR: You have not supplied any arguments.                                              

"
            echo " 

                                                                 ################################################################################
                                                                 #                                                                              #
                                                                 #             USAGE -- haplotype caller optimised for Cyclospora               #
                                                                 #                                                                              #
                                                                 ################################################################################

"
      echo "                    $package [ARGUMENTS]"
      echo " "
      echo "                    ARGUMENTS:"
      echo "                    -A                specify depth required to assign a haplotype to a specimen"
      echo "                    -C                specify the location of the typing_workflow_[version] directory"
      echo "                    -D                specify the directory containing all paired 'fastq.gz' files. PLEASE PROVIDE THE FULL PATH - eg. Users/your_name/Documents/etc/etc."
      echo "                    -h                show brief help"
      echo "                    -H                indicate if you would like to generate a haplotype sheet '-H Yes' or '-H Y' or '-H YES' or '-H yes' or '-H No' or '-H No' or '-H NO' or '-H no'."
      echo "                    -I                provide complete path to list of Illumina adapter sequences (fasta format only)"
      echo "                    -L                provide length of specimen names. This will truncate the name of your 'fastq.gz' files to length 'L'."
      echo "                    -R                amount of RAM (MB) you wish to allocate (only comes into effect during the clustering)."
      echo "                    -T                specify number of threads to use at various steps"
      echo "                    -N                specify depth required to add a new haplotype to your haplotype database"

      echo "

 "
      exit 0


done


# colon after a flag means it is a required flag.


Aflag=
Dflag=
Cflag=
Nflag=
Tflag=
Iflag=
Rflag=
Lflag=
Hflag=


while getopts ":N:T:A:D:C:I:R:L:H:h" opt; do
  case $opt in

    h)
            echo " 

                                                                 ################################################################################
                                                                 #                                                                              #
                                                                 #             USAGE -- haplotype caller optimised for Cyclospora               #
                                                                 #                                                                              #
                                                                 ################################################################################

"
      echo "                    $package [ARGUMENTS]"
      echo " "
      echo "                    ARGUMENTS:"
      echo "                    -A                specify depth required to assign a haplotype to a specimen"
      echo "                    -C                specify the location of the typing_workflow_[version] directory"
      echo "                    -D                specify the directory containing all paired 'fastq.gz' files. PLEASE PROVIDE THE FULL PATH - eg. Users/your_name/Documents/etc/etc."
      echo "                    -h                show brief help"
      echo "                    -H                indicate if you would like to generate a haplotype sheet '-H Yes' or '-H Y' or '-H YES' or '-H yes' or '-H No' or '-H No' or '-H NO' or '-H no'."
      echo "                    -I                provide complete path to list of Illumina adapter sequences (fasta format only)"
      echo "                    -L                provide length of specimen names. This will truncate the name of your 'fastq.gz' files to length 'L'."
      echo "                    -R                amount of RAM (MB) you wish to allocate (only comes into effect during the clustering)."
      echo "                    -T                specify number of threads to use at various steps"
      echo "                    -N                specify depth required to add a new haplotype to your haplotype database"

      echo "

 "
      exit 0

      ;;

      \?)
      echo "Invalid option: -$OPTARG." >&2
      exit 1
      ;;

    :)
      echo "ERROR: No argument was supplied for -$OPTARG." >&2
      exit 1
      ;;

    A)
      Aflag=$OPTARG
      ;;

    C)
      Cflag=$OPTARG
      ;;

    D)
      Dflag=$OPTARG
      ;;

    H)
      Hflag=$OPTARG
      ;;

    I)
      Iflag=$OPTARG
      ;;

    L)
      Lflag=$OPTARG
      ;;


    N)
      Nflag=$OPTARG
      ;;

    R)
      Rflag=$OPTARG
      ;;
  
    T)
      Tflag=$OPTARG
      ;;

  esac

done


echo ""
echo ""



#Check that there is a data location (D) flag

if [ -z "$Dflag" ]; then
    echo "ERROR: Option -D (location of paired 'fastq.gz' files) was not specified. Exiting."
exit
fi

if [ ! -d "$Dflag" ] 
then
    echo "ERROR: The directory $Dflag does not exist. Exiting." 
exit
fi




if [ -z "$Aflag" ]; then

      echo "No value specified for '-A'. A default value of 20 will be applied."
      Aflag=20


      elif [[ $Aflag = *"="* ]] || [[ $Aflag = *":"* ]] || [[ $Aflag = *";"* ]]|| [[ $Aflag = *","* ]] || [[ $Aflag = *"-"* ]] || [[ $Aflag = *"--"* ]]; then
      echo "ERROR: Please do not provide an '=' sign, colon ':' or any other separator in front of your arguments. Exiting."
      exit 1

      elif [[ $Aflag != *[[:digit:]]* ]]; then
      echo "ERROR: Please provide a numeric value for argument '-A' (min. depth required assign haplotypes to specimens). Exiting."
      exit 1

      fi




if [ -z "$Nflag" ]; then
    

      echo "No value specified for '-N'. A default value of 100 will be applied."
      Nflag=100

      elif [[ $Nflag = *"="* ]] || [[ $Nflag = *":"* ]] || [[ $Nflag = *";"* ]]|| [[ $Nflag = *","* ]] || [[ $Nflag = *"-"* ]] || [[ $Nflag = *"--"* ]]; then
      echo "ERROR: Please do not provide an '=' sign, colon ':' or any other separator in front of your arguments. Exiting."
      exit 1

      elif [[ $Nflag != *[[:digit:]]* ]]; then
      echo "ERROR: Please provide a numeric value for argument '-N' (min. depth required to write a newly discovered haplotype to the database). Exiting."
      exit 1

      fi



if [ -z "$Tflag" ]; then

echo "No value specified for '-T'. A default value of 10 threads will be applied."

Tflag=10

      elif [[ $Tflag = *"="* ]] || [[ $Tflag = *":"* ]] || [[ $Tflag = *";"* ]]|| [[ $Tflag = *","* ]] || [[ $Tflag = *"-"* ]] || [[ $Tflag = *"--"* ]]; then
      echo "ERROR: Please do not provide an '=' sign, colon ':' or any other separator in front of your arguments. Exiting."
      exit 1

      elif [[ $Tflag != *[[:digit:]]* ]]; then
      echo "ERROR: Please provide a numeric value for argument '-T' (number of threads). Exiting."
      exit 1

      fi




if [ -z "$Cflag" ]; then
    echo "ERROR: Option -C (location of software directory) was not specified. Exiting."
exit
fi

if [ ! -d "$Cflag/$alpha_version" ] 
then
    echo "ERROR: The directory $Cflag/$alpha_version does not exist. 
This is because you have either, (a) changed the original name of the software directory, or (b) have provided the wrong directory location.
Exiting." 
exit
fi



if [ -z "$Iflag" ]; then
    echo "No value specified for '-I'. A default set of Illumina adapters will be used for adapter trimming."

Iflag=$Cflag/$alpha_version/$HAPLOTYPE_CALLER/REF_SEQS/TRIMMING/Illumina_adapters.fasta

fi

if [ ! -f "$Iflag" ] 
then
    echo "ERROR: The file $Iflag does not exist. Please provide the complete path to your list of Illumina adapter sequences in fasta format. Exiting." 
exit
fi


if [ -z "$Rflag" ]; then


echo "No value specified for '-R'. A default value of 1000 MB will be allocated for clustering."
Rflag=1000


      elif [[ $Rflag != *[[:digit:]]* ]]; then
      echo "ERROR: Please provide a numeric value for argument '-R' (amount of RAM allocated for clustering). Exiting."
      exit 1

      fi



if [ -z "$Lflag" ]; then

echo "No value specified for '-L'. A default specimen name length of 11 characters will be applied."
Lflag=11


      elif [[ $Lflag != *[[:digit:]]* ]]; then
      echo "ERROR: Please provide a numeric value for argument '-L' (length of specimen names). Exiting."
      exit 1

      fi



if [ -z "$Hflag" ]; then

      echo "No value specified for '-H'. The default option is to NOT print a haplotype sheet. No haplotype sheet will be printed."
      Hflag=N


      elif [[ "$Hflag" == "N" ]] || [[ "$Hflag" == "No" ]] || [[ "$Hflag" == "NO" ]] || [[ "$Hflag" == "no" ]]; then

      echo "You have opted NOT to print a haplotype sheet at the completion of this workflow."

      elif [[ "$Hflag" == "Y" ]] || [[ "$Hflag" == "Yes" ]] || [[ "$Hflag" == "YES" ]] || [[ "$Hflag" == "yes" ]]; then
 
      echo "You HAVE opted to print a haplotype sheet at the completion of this workflow."

      else

      echo "ERROR: You have provided an invalid argument for '-H'. Valid options include: Yes, Y, YES, yes, No, N, NO, no. Exiting."

      exit 1

fi

echo ""
echo ""







#### IF YOU FIND A NEW HAPLOTYPE, WHAT IS THE MINIMUM DEPTH REQUIRED FOR IT TO BE CONSIDERED REAL? THIS IS SPECIFICALLY FOR *NEW* HAPLOTYPE DISCOVERY.

REQUIRED_DEPTH_NEW_HAP=$Nflag  ## a default of 100 is provided if no argument is supplied for this flag.



##### REQUIRED DEPTH OF SEQUENCING TO SUPPORT THE PRESENCE OF A HAPLOTYPE WITHIN A SPECIMEN

MINIMUM_DEPTH_TO_ASSIGN_HAPLPLOTYPES_TO_SPECIMENS=$Aflag  #### I set this to 20 as a default -- remember this is only the absolute minimum.




working_directory=$Cflag/$alpha_version/


illumina_adapaters=$Iflag                                          
                                                                                                                                               

input_reads=$Dflag


NUMBER_OF_THREADS=$Tflag


LENGTH_OF_SPECIMEN_NAMES=$Lflag


AMOUNT_OF_RAM=$Rflag



# REFERENCE FILE 1
read_recovery_refs=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/READ_RECOVERY/MAPPING_NON_JUNCTION_REFERENCES_WITH_PRIMER_FEB_2020.fasta


# REFERENCE FILE 2
read_recovery_refs_JUNCTION=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/READ_RECOVERY/MAPPING_JUNCTION_WITH_PRIMERS_FEB_2020.fasta


# REFERENCE FILE 3
complete_reference_database=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/BLASTING/ORIGINAL_REFS/CYCLOSPORA_REFERENCES_SHORTENED_FEB_2020.fasta


# REFERENCE FILE 4
bed_references=$working_directory$HAPLOTYPE_CALLER/REF_SEQS/BED_FILE/CYCLOSPORA_FEB_11_2020.bed

#IF YOU MADE THIS BED FILE USING XLSX YOU WILL NEED TO RUN DOS2UNIX ON THE FILE TO CONVERT IT OR THE SCRIPT WONT RUN CORRECTLY!




LOC=$working_directory$HAPLOTYPE_CALLER

if [ ! -d "$LOC" ] 
then
    echo "ERROR: The directory $LOC does not exist. Please check that the correct path was supplied for your software directory. Exiting." 
exit
fi



cd $LOC

rm -rf TMP

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
                    echo $NUMBER_OF_THREADS > $LOC/TMP/THREADS_TO_USE
                    echo $LENGTH_OF_SPECIMEN_NAMES > $LOC/TMP/NAME_LENGTH
                    echo $AMOUNT_OF_RAM > $LOC/TMP/RAM_ALLOCATION


mkdir R1_FILES
mkdir R2_FILES


if [ ! -d "$input_reads" ] 
then
    echo "ERROR: The directory $input_reads does not exist. This can happen if you did not supply the complete/correct path to your Illumina reads. Exiting." 
exit
fi

cd $input_reads


cp *_R1_* $LOC/TMP/R1_FILES
cp *_R2_* $LOC/TMP/R2_FILES

cd $LOC/TMP/


LOC=`cat RUN_DIR`
LENGTH_OF_SPECIMEN_NAMES=`cat NAME_LENGTH`



cd $LOC/TMP/R1_FILES
##RENAME THE FASTQ FILES IN R1 FOLDER
   for file in *
                do
                     mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e "s/\(.\{$LENGTH_OF_SPECIMEN_NAMES\}\).*.\(.\{3\}\)/\1.R1.fastq\2/"`
                done

##MAKE LIST OF SPECIMENS TO GENOTYPE
   for file in *
                 do
                      echo "$file" | sed -e 's/-/_/g' | sed -e "s/^\(.\{$LENGTH_OF_SPECIMEN_NAMES\}\).*/\1/" >> $LOC/TMP/SPECIMENS_TO_TRIM   

                 done



cd $LOC/TMP/R2_FILES
##RENAME THE FASTQ FILES IN R2 FOLDER
   for file in *
                  do
                       mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e "s/\(.\{$LENGTH_OF_SPECIMEN_NAMES\}\).*.\(.\{3\}\)/\1.R2.fastq\2/"` 

                  done

cd $LOC/TMP


mkdir PROCESSED_READS



######    TRIM THE R1 AND R2 READS   #####################################################################################################################

for SPECIMEN_NAME in `cat $LOC/TMP/SPECIMENS_TO_TRIM`

           do

                    ##check that the R1 and R2 read files exist and that they copied over correctly. If not, exit.

                    if [ ! -f "$LOC/TMP/R1_FILES/$SPECIMEN_NAME.R1.fastq.gz" ] 
                    then
                        echo "ERROR: Your 'R1' file for specimens $SPECIMEN_NAME... does not exist. Please check that the location $input_reads contains your a set of paired (R1 and R2) 'fastq.gz' files and check that all files are in 'fastq.gz' format. Exiting." 
                    exit
                    fi

                    if [ ! -f "$LOC/TMP/R2_FILES/$SPECIMEN_NAME.R2.fastq.gz" ] 
                    then
                        echo "ERROR: Your 'R2' file for specimens $SPECIMEN_NAME... does not exist. Please check that the location $input_reads contains your a set of paired (R1 and R2) 'fastq.gz' files. Exiting." 
                    exit
                    fi


                                        # TRIM USING BBDUK

                                        bash $LOC/BBMAP/bbmap/bbduk.sh \
                                        in1=$LOC/TMP/R1_FILES/$SPECIMEN_NAME.R1.fastq.gz \
                                        in2=$LOC/TMP/R2_FILES/$SPECIMEN_NAME.R2.fastq.gz \
                                        out1=$LOC/TMP/$SPECIMEN_NAME.clean1.fq \
                                        out2=$LOC/TMP/$SPECIMEN_NAME.clean2.fq \
                                        minlen=50 \
                                        ref=$illumina_adapaters \
                                        ktrim=r \
                                        k=23 \
                                        mink=11 \
                                        hdist=1
                                        

                                        # MERGE USING BBMERGE

                                        bash $LOC/BBMAP/bbmap/bbmerge-auto.sh \
                                        in1=$LOC/TMP/$SPECIMEN_NAME.clean1.fq \
                                        in2=$LOC/TMP/$SPECIMEN_NAME.clean2.fq \
                                        out=$LOC/TMP/$SPECIMEN_NAME.merged.fq \
                                        outu=$LOC/TMP/$SPECIMEN_NAME.unmerged.fq \
                                        rem k=62 \
                                        extend2=50 ecct
                    

                    mv $LOC/TMP/$SPECIMEN_NAME.merged.fq $LOC/TMP/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq
                    cat $LOC/TMP/$SPECIMEN_NAME.unmerged.fq >> $LOC/TMP/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq
                    rm -rf $LOC/TMP/$SPECIMEN_NAME.clean1.fq $LOC/TMP/$SPECIMEN_NAME.clean2.fq $LOC/TMP/$SPECIMEN_NAME.unmerged.fq
 

             done







cd $LOC/TMP/

rm -rf R1_FILES R2_FILES



bash ../Finding_New_Haps_SNP_Based_Module.sh
bash ../Finding_New_Junction_Types_Module-v4.sh
bash ../Assign_Haplotypes.sh



cd $LOC/
rm -rf TMP


               if [[ "$Hflag" == "Yes" ]] ||  [[ "$Hflag" == "Y" ]] ||  [[ "$Hflag" == "YES" ]] || [[ "$Hflag" == "yes" ]]; then

                    echo "

                    The option '-H Yes' was selected. A new haplotype sheet WILL be printed momentarily. Depding on the number of specimens in your dataset, this could take several minutes.

                    "

                    Rscript START_haplotype_sheet_generator.R

                    echo "Haplotype sheet printing complete."

               elif [[ "$Hflag" == "No" ]] ||  [[ "$Hflag" == "N" ]] ||  [[ "$Hflag" == "NO" ]] || [[ "$Hflag" == "no" ]]; then

                    echo "
The '-H NO' option was selected. A new haplotype sheet will NOT be printed.

                    "

               fi

cd $LOC/



rm -rf TMP_TYPES

echo "The Cyclospora haplotype calling workflow has run to completion.
"


