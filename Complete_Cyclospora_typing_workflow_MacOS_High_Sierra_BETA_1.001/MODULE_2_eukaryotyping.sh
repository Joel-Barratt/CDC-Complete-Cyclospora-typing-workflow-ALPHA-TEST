#!/bin/bash


#####################################################################################################################################################################
CYCLONE=/CYCLONE_MacOS_High_Sierra_BETA_Cyclospora_1.001  # do not modify this line        ##########################################################################
#####################################################################################################################################################################   


#Tell me the folder where you pasted and extracted the CYCLONE zip file.
working_directory=/Users/joelbarratt/Documents/CYCLOSPORA


## How many threads would you like to use.
number_of_threads=10


### epsilon must be any number between 0 and 1.
epsilon=0.3072









#DO NOT MODIFY BELOW THIS LINE.
#####################################################################################################################################################################   


cd $working_directory/$CYCLONE_DIRECTORY/EUKARYOTYPING

echo $number_of_threads > THREADS
echo $epsilon > EPSILON

Rscript run.r

echo "EUKARYOTYPING COMPLETE!"

rm THREADS
rm EPSILON

cd $working_directory/$CYCLONE_DIRECTORY/

