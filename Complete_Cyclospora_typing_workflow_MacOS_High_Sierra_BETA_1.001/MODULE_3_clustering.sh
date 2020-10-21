
#!/bin/bash

###################################################################################################################################
# do not modify the NEXT 2 LINES lines                                                      #######################################
CYCLONE=/CYCLONE_MacOS_High_Sierra_BETA_Cyclospora_1.001/                                   #######################################
CLUSTERING=$CYCLONE/CLUSTERING                                                              #######################################
###################################################################################################################################

###USER MUST MODIFY THE FOLLOWING LINES:



### TELL ME THE DIRECTORY WHERE YOU PASTED THE CYCLONE FOLDER. LITERALLY WHERE YOU UNZIPPED CYCLONE AND INTEND TO RUN/INSTALL IT.
cyclone_location=/Users/joelbarratt/Documents/CYCLOSPORA


# stringency this must be a number between 1 and 100
stringency=95



### what number of clusters do you want to start at?
cluster_min=5


#what number of clusters do you want to finish at?
cluster_max=50


#tell me the name of the reference list file
your_list_of_reference_clusters=2018_gold_clusters.txt




















#### STOP - DO NOT MODIFY ANYTHING BELOW THIS POINT.

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




Rscript $LOC/CLUSTER_FINDER.R  ##thi script is good, but it finishes the workflow just before generating the state process reports and cluster reports, and the cluster log. just need to make these and you are sweet!

## go to the end of the cluster finder sheet to find where the state reports script kicks off. Then you will want to modify the state process reports.

cd $LOC

rm -rf TMP_REP

cd $cyclone_location/$CYCLONE

echo "CLUSTERING COMPLETE!"




