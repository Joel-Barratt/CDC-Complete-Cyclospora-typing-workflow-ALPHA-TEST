#UPDATED August 5, 2021

stringency_level = as.numeric(readLines("STRINGENCY"))
number_of_threads = as.numeric(readLines("THREADS"))



gold_standard_clusters_A <- cluster_standards
colnames(gold_standard_clusters_A) <- c("Seq_ID", "Standard_number_A")
gold_standard_clusters_B <- gold_standard_clusters_A
colnames(gold_standard_clusters_B) <- c("row", "Standard_number_B")
list <- dist2list(as.dist(matrix))
attach(list)
sorted_list <- list[order(value),] 
sorted_list <- as.data.frame(sorted_list)


colnames(sorted_list) = NULL
colnames(sorted_list) <- c("Seq_ID", "row", "distance")
total <-merge(sorted_list, gold_standard_clusters_A, by=c("Seq_ID"), all.x=TRUE)
total_final <- merge(total, gold_standard_clusters_B, by=c("row"), all.x=TRUE)
last_column <- paste0(total_final$Standard_number_A,"plus",total_final$Standard_number_B,"_")
set_gold_standard <- cbind(total_final, last_column)




cluster_to_cluster_A <- set_gold_standard %>% filter(str_detect(Standard_number_A, "NA") == FALSE)
cluster_to_cluster_A_and_B <- cluster_to_cluster_A %>% filter(str_detect(Standard_number_B, "NA") == FALSE)


within_cluster_standard_distance <- NULL

for(n in unique(cluster_standards[,2])){
	

current_standard <- paste(n,"plus",n,"_",sep="")

standard <- cluster_to_cluster_A_and_B %>% filter(str_detect(last_column, current_standard) == TRUE)

within_cluster_standard_distance <- rbind(within_cluster_standard_distance, standard)
	
	
}



golden_cluster_distance <- mean(within_cluster_standard_distance$distance) + 3*(sd(within_cluster_standard_distance$distance))
date_today <- Sys.Date()


Ensemble <- as.matrix(matrix)
Ensemble_x <- as.hclust(agnes(x=Ensemble, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))


check_all_cluster_numbers = mclapply(min_cluster_number:max_cluster_number, function (f) {

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=f))))
colnames(cluster) <- NULL
cluster <- cbind(rownames(cluster), cluster)
colnames(cluster) <- c("Seq_ID", "clus_A")


cluster_first_col <-merge(sorted_list, cluster, by=c("Seq_ID"), all.x=TRUE)
colnames(cluster) <- NULL
colnames(cluster) <- c("row", "clus_B")
cluster_second_col <-merge(cluster_first_col, cluster, by=c("row"), all.x=TRUE)
cluster_calc <- subset(cluster_second_col, clus_A == clus_B)

distances_below_gold_cluster_distance <- cluster_calc %>% filter(distance < golden_cluster_distance)

percent_within_cluster_distances_meeting_cutoff <- ((length(distances_below_gold_cluster_distance$distance))/length(cluster_calc$distance))*100 # this will be assigned as the output of the function.

}, mc.cores= number_of_threads)


names(check_all_cluster_numbers) <- c(paste0("clus_",min_cluster_number:max_cluster_number))


for (f in names(check_all_cluster_numbers)){


					if (stringency_level > check_all_cluster_numbers[[f]]) {
						
                                fix <- as.numeric(str_remove(f, "clus_"))
								print(paste("Searching for the most appropriate cluster number -",fix,"clusters is too small."))
	
												} else {
													
								fix <- as.numeric(str_remove(f, "clus_"))
								name_clus <- paste(date_today,"_there_are_",fix,"_genetic_clusters.txt", sep="")
								cluster2 <- as.data.frame(melt(factor(cutree(Ensemble_x, k=fix))))
								colnames(cluster2) <- NULL
								cluster2 <- cbind(rownames(cluster2), cluster2)
								colnames(cluster2) <- c("Seq_ID", "Assigned_cluster")


							write.table(cluster2, name_clus, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
							
						print(paste("The most likely number of clusters is:", fix))
						
										if (check_all_cluster_numbers[[f]] >= stringency_level) {
											
											break
										}
					    }
	}


correct_number_of_clusters <- fix


