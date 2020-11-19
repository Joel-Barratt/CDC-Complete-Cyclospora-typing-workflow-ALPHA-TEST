#UP November 19, 2020

setwd("../haplotype_sheets")
## now you need to list the files in order, take the top line and then import that below.

myFiles <- list.files(pattern="*txt", all.files = FALSE)

myFiles <- (sort(myFiles, decreasing = TRUE)[1])

data = read.csv(myFiles, skip=0, stringsAsFactors = FALSE, sep = "\t")

locinames_and_haps <- colnames(data)
locinames_and_haps <- locinames_and_haps[2:length(locinames_and_haps)]

locinames <- unique(str_replace(str_remove(locinames_and_haps, "_Hap_.*"), "Mt_*.*_Junction", "Mt_Cmt"))

almost_base <- str_remove(str_remove(locinames, "_Hap_.*"), "_PART_.*")
locinames_base <- unique(str_replace(almost_base, "Mt_*.*_Junction", "Mt_Cmt"))


## below will remove columns where only one type exists (where distance will always be zero. This is because these sites are not informative)
## data transformation wont work if you dont have more than one type for a given locus.

for(c in locinames){
	
if (length(grep(c, locinames_and_haps))	== 1){
	
	
data <- data %>% select(-contains(c))

}

}


##second round after the non-diverse regions are deleted

locinames_and_haps <- colnames(data)
locinames_and_haps <- locinames_and_haps[2:length(locinames_and_haps)]
locinames <- unique(str_replace(str_remove(locinames_and_haps, "_Hap_.*"), "Mt_*.*_Junction", "Mt_Cmt"))
almost_base <- str_remove(str_remove(locinames, "_Hap_.*"), "_PART_.*")
locinames_base <- unique(str_replace(almost_base, "Mt_*.*_Junction", "Mt_Cmt"))


#### find ploidy of markers
organelle <- as.data.frame(cbind(locinames , str_remove(locinames, "_.*")))
organelle$V2 <- str_replace(organelle$V2, "Mt", "1")
organelle$V2 <- str_replace(organelle$V2, "Nu", "2")

ploidy = c(as.numeric(organelle$V2))


###########The code above is not in the eukaryotyping scripts online.

data = data[!is.na(data$Seq_ID) & data$Seq_ID != "",]
ids = data$Seq_ID
nids = length(ids)
nloci = length(locinames)

### transform data 

newdata = c()
datacompleteness = c()
# MaxMOI for each locus 
for (j in 1:nloci) {
	locicolumns = grepl(paste(locinames[j],"",sep=""),colnames(data))
	raw_alleles = data[,locicolumns]
	maxMOI = max(c(2,rowSums(raw_alleles == "X",na.rm=TRUE)))
	MOI = rowSums(raw_alleles == "X",na.rm=TRUE)
	nalleles = sum(locicolumns,na.rm=TRUE)
	newdatatemp = rbind(matrix(NA,nids,maxMOI))
	sapply(1:nids,function (x) if (MOI[x] > 0) { newdatatemp[x,1:MOI[x]] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")})
	if (ploidy[j] > 1) {
		sapply(1:nids,function (x) if (MOI[x] == 1) { newdatatemp[x,1:2] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")})
	}
	colnames(newdatatemp ) = paste(1:maxMOI,"_Hap_",locinames[j],sep="")
	newdata = cbind(newdata,newdatatemp )
	datacompleteness  = cbind(datacompleteness,MOI)
}

data = cbind(ids,newdata)
data = data.frame(data)

# calculate locus-specific amplification success
datacompleteness_bylocus = sapply(locinames_base,function (x) rowSums(cbind(datacompleteness[,grepl(x,locinames)]))>0)

#print("Number of loci sequenced")
#table(rowSums(datacompleteness_bylocus))



#### USE THE NEXT FEW LINES TO DETERMINE THE MINIMUM CRITERIA FOR RETAINING THE SPECIMENS FOR ANALYSIS #################################################################

#"Mt_Cmt"   "Mt_MSR"   "Nu_360i2" "Nu_378"   "Nu_CDS1"  "Nu_CDS2"  "Nu_CDS3" "Nu_CDS4"

cleandata = data[(rowSums(datacompleteness_bylocus[,c("Mt_Cmt","Mt_MSR","Nu_360i2")]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | rowSums(datacompleteness_bylocus) >= 5 | (rowSums(datacompleteness_bylocus[,c("Mt_Cmt","Mt_MSR","Nu_378")]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c("Mt_MSR","Nu_360i2","Nu_378")]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c("Mt_Cmt","Nu_360i2","Nu_378")]) == 3 & rowSums(datacompleteness_bylocus) >= 4),] 

#########################################################################################################################################################################




data = cleandata
ids = data$ids
nids = length(ids)

#write.csv(data,"transformed_data.csv")

setwd("../EUKARYOTYPING/")