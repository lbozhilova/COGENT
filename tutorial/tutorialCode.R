#===== COGENT tutorial code =====#

## Introduction

# require("devtools")
# install_github("lbozhilova/COGENT")

#-----
require("COGENT")
require("ggplot2")
require("ggthemes")
require("parallel")
set.seed(101019)

## Gene expression data

load("tutorialData.RData")
head(yeastData[,1:6])

#-----
checkExpressionDF(yeastData)

## Network construction functions and network comparison

#-----
buildPearson <- function(df, quant=0.90){
   check <- checkExpressionDF(df)
   A <- cor(t(df[,colnames(df)!="Name"]), use="pairwise.complete.obs", method="pearson")
   threshold <- quantile(A, quant, na.rm=TRUE)
   A <- 1*(A>=threshold); diag(A) <- 0
   colnames(A) <- rownames(A) <- df$Name
   return(A)
}
pearsonA <- buildPearson(yeastData)

#-----
buildKendall<- function(df, quant=0.90){
   check <- checkExpressionDF(df)
   A <- cor(t(df[,colnames(df)!="Name"]), use="pairwise.complete.obs", method="kendall")
   threshold <- quantile(A, quant, na.rm=TRUE)
   A <- 1*(A>=threshold); diag(A) <- 0
   colnames(A) <- rownames(A) <- df$Name
   return(A)
}
kendallA <- buildKendall(yeastData)

#-----
PKcomparison <- getEdgeSimilarity(list(pearsonA, kendallA), align=FALSE)
PKcomparison$nodeCount

#-----
PKcomparison$globalSimilarity

#-----
hist(PKcomparison$localSimilarity,
     main="Local similarity between Pearson and Kendall networks",
     xlab="Similarity",
     breaks=seq(0, 1, .05), col="cornflowerblue")

## COGENT analysis

calculateDegrees <- function(A){
   deg <- rowSums(A)
   names(deg) <- colnames(A)
   return(deg)
}

#-----
x <- cogentSingle(df=yeastData, netwkFun=buildPearson, propShared=0.50)
x$globalSimilarity

#-----
x <- cogentSingle(df=yeastData, netwkFun=buildKendall, propShared=0.50,
                  nodeFun=calculateDegrees, nodeModes="all",
                  use="pairwise.complete.obs", method="pearson",
                  k.or.p=0.10,
                  scale=TRUE)
c("globalSimilarity"=x$globalSimilarity, "corDegrees"=x$corSimilarity)

#-----
stabilityPearson <- cogentLinear(df=yeastData, netwkFun=buildPearson, propShared=0.50,
                                 repCount=100, nodeFun=calculateDegrees, nodeModes="all",
                                 use="pairwise.complete.obs", method="pearson", 
                                 k.or.p=0.10, scale=TRUE)
head(stabilityPearson)

#-----
stabilityKendall <- cogentParallel(df=yeastData, netwkFun=buildKendall, propShared=0.50, 
                                   repCount=100, threadCount=4, nodeFun=calculateDegrees, 
                                   nodeModes="all", use="pairwise.complete.obs", 
                                   method="pearson", k.or.p=0.10, scale=TRUE)
head(stabilityKendall)

#-----
stabilityBoth <- rbind(
  cbind(stabilityPearson, method="Pearson"),
  cbind(stabilityKendall, method="Kendall")
)
ggplot(stabilityBoth, aes(x=method, y=globalSimilarity)) +
  theme_economist_white() +
  geom_boxplot() +
  ggtitle("Stability of Pearson and Kendall\ncorrelation networks") +
  scale_x_discrete("Method") +
  scale_y_continuous("Global similarity")

#-----
ggplot(stabilityBoth, aes(x=method, y=corSimilarity)) +
  theme_economist_white() +
  geom_boxplot() +
  ggtitle("Stability of Pearson and Kendall\ncorrelation networks") +
  scale_x_discrete("Method") +
  scale_y_continuous("Degree consistency via correlations")

## Choosing a threshold

getThresholdStability <- function(th){
   dfSplit <- splitExpressionData(yeastData, propShared=0)
   A <- lapply(dfSplit, function(df) buildPearson(df, th))
   return(getEdgeSimilarityCorrected(A, type="expected")) 
}

#-----
aggregateThresholdStability <- function(th, repCount=100){
   thS <- replicate(repCount, getThresholdStability(th), simplify=FALSE)
   thS <- do.call("rbind", thS); thS <- apply(thS, 2, unlist)
   return(as.data.frame(cbind(thS, threshold=th)))
}
thresholds <- seq(0.5, 0.99, 0.01)
thresholdComparisonDF <- mclapply(thresholds, aggregateThresholdStability, mc.cores=6)
thresholdComparisonDF <- do.call("rbind", thresholdComparisonDF)

#-----
thresholdComparisonDF <- subset(thresholdComparisonDF, 
                                !is.na(thresholdComparisonDF$correctedSimilarity))
ggplot(thresholdComparisonDF, aes(x=1-threshold, y=correctedSimilarity)) +
   geom_smooth() +
   theme_economist_white() +
   ggtitle("Threshold choice for Pearson correlation networks") +
   scale_y_continuous("Density-adjusted consistency") +
   scale_x_continuous("Network density", breaks=seq(0, 0.5, .05))

