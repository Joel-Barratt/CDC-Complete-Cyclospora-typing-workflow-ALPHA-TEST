#!/bin/bash

#UP November 19, 2020

#####################################################################################################################################################################
CYCLONE_DIRECTORY=/Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_BETA_1.001  # do not modify this line   ####################################################
#####################################################################################################################################################################   


#Tell me the folder where you pasted and extracted the Complete genotyping workflow zip file.
working_directory=/Users/joelbarratt/Documents/CYCLOSPORA


## How many threads would you like to use.
number_of_threads=11


### Epsilon. This must be any number between 0 and 1.
epsilon=0.3072













# DO NOT MODIFY BELOW THIS POINT
#####################################################################################################################################################################   


cd $working_directory/$CYCLONE_DIRECTORY/EUKARYOTYPING

echo $number_of_threads > THREADS
echo $epsilon > EPSILON

Rscript run.r

echo "EUKARYOTYPING COMPLETE!"

rm THREADS
rm EPSILON

cd $working_directory/$CYCLONE_DIRECTORY/

