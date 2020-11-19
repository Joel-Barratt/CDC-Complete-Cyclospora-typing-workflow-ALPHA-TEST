library(cluster)
library(tidyverse)
library(quantmod)
library(dplyr)
library(stringr)
library(filesstrings)
library(reshape)
library(dbscan)
library(spaa)
require(parallel)
require(stringr)

#UP November 19, 2020

date_today <- Sys.Date()

temporary_directory <- getwd()




#### First you need to import the most recent ensemble matrix as a variable in R called "matrix" in R

matrix_location <- readLines("MATRIX_LOCATION")

setwd(matrix_location)

myFiles <- list.files(pattern="*csv", all.files = FALSE)
myMatrix <- (sort(myFiles, decreasing = TRUE)[1])

matrix <- as.matrix(read.table(myMatrix, sep = ",", row.names=1, header=TRUE))


setwd(temporary_directory)

##### Now you need to calculate the number of genetic clusters using the "gold standard" list. Use this to calculate the number of genetic clusters. "number_of_genetic_clusters" needs to be the name of this variable.

import_gold_standards <- readLines("GOLD_CLUSTERS")
cluster_standards <- read.table(import_gold_standards, head=TRUE)  #THIS IS WHERE YOU DEFINE WHICH Seq_IDs are your gold standard "within cluster" distances

###set range of cluster numbers to test
min_cluster_number = as.numeric(readLines("CLUSTER_MIN"))
max_cluster_number = as.numeric(readLines("CLUSTER_MAX"))



source("../cluster_counter.R")



Ensemble <- matrix
Ensemble_x <- as.hclust(agnes(x=Ensemble, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))
cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=correct_number_of_clusters))))
colnames(cluster) <- NULL
cluster <- cbind(rownames(cluster), cluster)
colnames(cluster) <- c("Seq_ID", "Cluster")

setwd("..")

write.table(cluster, paste("../clusters_detected/",date_today,"_RESULTING_CLUSTERS_",f,".txt", sep=""), col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, na="")



