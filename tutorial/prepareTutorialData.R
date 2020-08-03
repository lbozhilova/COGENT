################################
##### COGENT tutorial data #####
################################

# First make sure you have downloaded and unzipped the data. You can find it here:
# https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5174/Downloads

# Parse the data
setwd("tutorial/")
yeastData <- read.table("E-MTAB-5174-query-results.tsv", sep="\t", header=TRUE)

# Remove unnecessary metadata columns
yeastData <- yeastData[,-1]

# Remove p-values
yeastData <- yeastData[,-seq(3, 53, 2)]

# Rename the "Gene.ID" column to "Name"
colnames(yeastData)[1] <- "Name"

# Remove gene symbols which occur on more than one row (none in this case)
badGenes <- names(table(yeastData$Name)[table(yeastData$Name)>1])
yeastData <- subset(yeastData, !(Name %in% badGenes))

# Remove genes with more than 25% of data points missing
fracAllowed <- 0.25
badGenes <- which(rowSums(is.na(yeastData[,-1]))>fracAllowed*(ncol(yeastData)-1)) # 1183 genes
yeastData <- yeastData[-badGenes,] # 1920 left

# Reduce to the 200 most variable genes, as measured by their median absolute deviation
geneMADs <- apply(yeastData[,-1], 1, mad, na.rm=TRUE)
names(geneMADs) <- yeastData$Name
geneMADs <- sort(geneMADs, decreasing=TRUE)
goodGenes <- names(geneMADs[1:200])
yeastData <- subset(yeastData, Name %in% goodGenes)

# Remove sample names
colnames(yeastData)[-1] <- paste0("sample", 1:(ncol(yeastData)-1))

# At this stage we won't do any data normalisation or apply additional filters
# Save results
save(yeastData, file="tutorialData.RData")
