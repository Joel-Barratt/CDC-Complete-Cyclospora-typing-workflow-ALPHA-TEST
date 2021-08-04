
#!/bin/bash

#UP November 19, 2020

###################################################################################################################################
# do not modify the next 2 lines.                                                           #######################################
CYCLONE=/Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_BETA_1.001/                  #######################################
CLUSTERING=$CYCLONE/CLUSTERING                                                              #######################################
###################################################################################################################################

###USER MUST MODIFY THE FOLLOWING LINES:



# TELL ME THE DIRECTORY WHERE YOU PASTED THE CYCLONE FOLDER. LITERALLY WHERE YOU UNZIPPED CYCLONE AND INTEND TO RUN/INSTALL IT.
cyclone_location=/Users/joelbarratt/Documents/CYCLOSPORA


# Stringency. This must be a number between 1 and 100
stringency=95


# What number of clusters do you want to start at?
cluster_min=1


# What number of clusters do you want to finish at?
cluster_max=50


# Tell me the name of the reference list file - be sure to past this document in the REFERENCE_CLUSTER_LIST directory!
your_list_of_reference_clusters=2018_gold_clusters.txt


# Tell me the number of threads you would like to use
number_of_threads=11





















#########################################################   DO NOT MODIFY BELOW THIS POINT


gold_standard_clusters=$cyclone_location/$CYCLONE/REFERENCE_CLUSTER_LIST/$your_list_of_reference_clusters


LOC=$cyclone_location$CLUSTERING

matrix_folder=$cyclone_location/$CYCLONE/ensemble_matrices


###########################################################################################################################################################
######                                                  ###################################################################################################
######  WRITE DIRECTORIES FOR TMP FILES & VARIABLES     ###################################################################################################
######                                                  ###################################################################################################
###########################################################################################################################################################


cd $LOC

rm -rf TMP_REP

mkdir TMP_REP
cd TMP_REP

echo $most_recent_ensemble_matrix > MATRIX
echo $gold_standard_clusters > GOLD_CLUSTERS
echo $stringency > STRINGENCY
echo $cluster_min > CLUSTER_MIN
echo $cluster_max > CLUSTER_MAX
echo $matrix_folder > MATRIX_LOCATION
echo $number_of_threads > THREADS




Rscript $LOC/CLUSTER_FINDER.R

cd $LOC

rm -rf TMP_REP

cd $cyclone_location/$CYCLONE

echo "CLUSTERING COMPLETE!"




