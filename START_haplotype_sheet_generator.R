library(filesstrings)
library(dplyr)


setwd('./SPECIMEN_GENOTYPES')
Seq_ID <- list.files(all.files = FALSE)

new_db <- ""


### Remove files with size of zero from list
for (combine in Seq_ID){
	
if (file.info(combine)$size == 0){
	
Seq_ID <- Seq_ID[Seq_ID != combine]
	
}
	
	
}

#length(Seq_ID)


for (combine in Seq_ID){
	
blast_result <- read.table(combine, sep="\t",stringsAsFactors = FALSE)	

new_db <- rbind(new_db, blast_result)	
	
}

new_db <- as.data.frame(tail(new_db, - 1))
#head(new_db, 5)

markers <- sort(unique(new_db$V1))

genotype_sheet <- data.frame(matrix(ncol = (length(markers)), nrow = (length(Seq_ID))))



genotype_sheet <- cbind(Seq_ID, genotype_sheet)

rownames(genotype_sheet) <- c(Seq_ID)

colnames(genotype_sheet) <- c("Seq_ID", markers)


for(add_x in Seq_ID){
	
for(current_marker in colnames(genotype_sheet))	{
	
blast_result <- read.table(add_x, sep="\t",stringsAsFactors = FALSE)	

if (length(grep(current_marker, blast_result)) > 0) {

genotype_sheet[add_x, current_marker] <- "X"

         }
	
	}
	
}

#genotype_sheet <- as.matrix(genotype_sheet)





### Now add those with zero markers to the bottom. This is important for report generation.

Seq_ID <- list.files(all.files = FALSE)


for (combine in Seq_ID){
	
if (file.info(combine)$size == 0){

empty_rows <- c(combine, 1:(length(colnames(genotype_sheet))-1))

empty_rows[2:length(empty_rows)] <- NA

empty_rows <- t(as.data.frame(empty_rows))

colnames(empty_rows) <- colnames(genotype_sheet)
rownames(empty_rows) <- NULL
	
genotype_sheet <- rbind(genotype_sheet, empty_rows)

rownames(genotype_sheet) <- NULL

	
  }

}


genotype_sheet[is.na(genotype_sheet)] <- ""


setwd("..")

file_name_new <- paste("../haplotype_sheets/",Sys.Date(),"_Cyclospora_haplotype_data_sheet",".txt", sep="")

write.table(genotype_sheet, file_name_new, sep="\t", quote=FALSE, row.names=FALSE)


