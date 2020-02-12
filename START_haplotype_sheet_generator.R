library(filesstrings)


setwd('./RESULTS')
myFiles <- list.files(all.files = FALSE)

new_db <- ""

for (combine in myFiles){
	
blast_result <- read.table(combine, sep="\t",stringsAsFactors = FALSE)	

new_db <- rbind(new_db, blast_result)	
	
}

new_db <- as.data.frame(tail(new_db, - 1))
#head(new_db, 5)

markers <- sort(unique(new_db$V1))

genotype_sheet <- data.frame(matrix(ncol = (length(markers)), nrow = (length(myFiles))))


rownames(genotype_sheet) <- c(myFiles)
colnames(genotype_sheet) <- c(markers)


for(add_x in myFiles){
	
for(current_marker in colnames(genotype_sheet))	{
	
blast_result <- read.table(add_x, sep="\t",stringsAsFactors = FALSE)	

if (length(grep(current_marker, blast_result)) > 0) {

genotype_sheet[add_x, current_marker] <- "X"

         }
	
	}
	
}

genotype_sheet <- as.matrix(genotype_sheet)

genotype_sheet[is.na(genotype_sheet)] <- ""

file_name_new <- paste("../HAPLOTYPE_SHEETS/",Sys.Date(),"_Cyclospora_haplotype_data_sheet",".txt", sep="")



write.table(genotype_sheet, file_name_new, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)


