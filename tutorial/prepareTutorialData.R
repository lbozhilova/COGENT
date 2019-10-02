################################
##### COGENT tutorial data #####
################################

# First make sure you have installed and unzipped the data. You can find it here:
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

# Parse the data
setwd("tutorial/")
liverData <- as.data.frame(read.table("LiverFemale3600.csv", sep=",", header=TRUE))

# Remove unnecessary metadata columns
liverData <- liverData[,-c(1, 3:8)]

# Rename the "gene_symbol" column to "Name"
colnames(liverData)[1] <- "Name"

# Remove gene symbols which occur on more than one row
badGenes <- names(table(liverData$Name)[table(liverData$Name)>1])
liverData <- subset(liverData, !(Name %in% badGenes))

# Reduce to the 200 most variable genes, as measured by their median absolute deviation
geneMADs <- apply(liverData[,-1], 1, mad)
names(geneMADs) <- liverData$Name
geneMADs <- sort(geneMADs, decreasing=TRUE)
goodGenes <- names(geneMADs[1:200])
liverData <- subset(liverData, Name %in% goodGenes)

# At this stage we won't do any data normalisation or apply additional filters
# Save results
save(liverData, file="tutorialData.RData")
