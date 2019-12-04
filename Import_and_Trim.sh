#!/bin/bash

###Haplotype_caller_pipeline - version - 1



running_directory=/Users/joelbarratt/Documents/CYCLOSPORA/HAPLOTYPE_CALLER
echo $running_directory > RUN_DIR

junction_with_primer=$running_directory/MOCK_FOR_READ_RECOVERY_junction_with_primers.fasta ### this is for detection of new junction types - need to include the primers in this one
echo $junction_with_primer > JUNC_REF

complete_reference_database=$running_directory/mock.list.fasta ###trim the primers off the junction sequences for this one.
echo $complete_reference_database > ALL_REF

illumina_adapaters=$running_directory/Illumina_adapters.fasta



##########################################################################################################################################################
######                                       #############################################################################################################
######  MOVE AND RENAME THE R1 AND R2 READS  #############################################################################################################
######                                       #############################################################################################################
##########################################################################################################################################################

cd $running_directory

mkdir R1_FILES
mkdir R2_FILES
mkdir FINAL_GENOTYPING_RESULTS

cd SPECIMENS_TO_GENOTYPE

cp *_R1_* ../R1_FILES
cp *_R2_* ../R2_FILES


cd ..

cd R1_FILES
##RENAME THE FASTQ FILES IN R1 FOLDER
	for file in *
		do
			mv "$file" `echo "$file" | sed -e 's/\(.\{11\}\).*.\(.\{3\}\)/\1.R1.fastq\2/'`
		done

##MAKE LIST OF SPECIMENS TO GENOTYPE
	for file in *
		do
			echo "$file" | sed -e 's/^\(.\{11\}\).*/\1/' >> ../SPECIMENS_TO_RUN
		done

cd ../R2_FILES
##RENAME THE FASTQ FILES IN R2 FOLDER
	for file in *
		do
			mv "$file" `echo "$file" | sed -e 's/\(.\{11\}\).*.\(.\{3\}\)/\1.R2.fastq\2/'`

		done

cd $running_directory

##########################################################################################################################################################
######                               #####################################################################################################################
######    TRIM THE R1 AND R2 READS   #####################################################################################################################
######                               #####################################################################################################################
##########################################################################################################################################################

		for SPECIMEN_NAME in `cat SPECIMENS_TO_RUN`
			do
			#TRIM USING BBDUK

bash $running_directory/BBMAP/bbmap/bbduk.sh in1=$running_directory/R1_FILES/$SPECIMEN_NAME.R1.fastq.gz in2=$running_directory/R2_FILES/$SPECIMEN_NAME.R2.fastq.gz out1=$running_directory/$SPECIMEN_NAME.clean1.fq out2=$running_directory/$SPECIMEN_NAME.clean2.fq minlen=50 ref=$illumina_adapaters ktrim=r k=23 mink=11 hdist=1

			#MERGE USING BBMERGE

bash $running_directory/BBMAP/bbmap/bbmerge-auto.sh in1=$running_directory/$SPECIMEN_NAME.clean1.fq in2=$running_directory/$SPECIMEN_NAME.clean2.fq out=$running_directory/$SPECIMEN_NAME.merged.fq outu=$running_directory/$SPECIMEN_NAME.unmerged.fq rem k=62 extend2=50 ecct

			#NOTE FROM BBTOOLS MANUAL ADVIZING AGAINST QULITY TRIMMING PRIOR TO MERGING READS:
			#Adapter-trimming reads is fine and is recommended prior to running BBDuk, particularly if kmers will be used. ***Quality-trimming is 
			#usually not recommended*****, unless the library has major quality issues at the right end of reads resulting in a low merge rate, in which 
			#case only weak trimming (such as qtrim=r trimq=8) should be used.




						#CONCATENATE MERGED AND UNMERGED READS AND PLACE IN A 'PROCESSED_READS' DIRECTORY FOR DOWNSTREAM USE

						mkdir $running_directory/PROCESSED_READS

						mv $running_directory/$SPECIMEN_NAME.merged.fq $running_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq

						cat $running_directory/$SPECIMEN_NAME.unmerged.fq >> $running_directory/PROCESSED_READS/$SPECIMEN_NAME.clean_merged.fastq
						rm -rf $SPECIMEN_NAME.clean1.fq $SPECIMEN_NAME.clean2.fq $SPECIMEN_NAME.unmerged.fq

			done


#TIDY UP BY REMOVING FILES NO LONGER NEEDED
rm -rf R2_FILES R1_FILES SPECIMENS_TO_RUN

bash Finding_New_Junction_Types_Module.sh




