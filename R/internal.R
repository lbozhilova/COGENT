############################
##### COGENT: INTERNAL ######
############################

#----- Alignment -----#

#' Align adjacency matrices
#'
#' \code{alignMatrices()} reduces two adjacency matrices to their common nodes
#' and ordering them to match. The matrices can then be compared using different
#' similarity measures. Note alignment is not always necessary, since many
#' network construction algorithms will preserve the same node order across
#' different runs.
#'
#' @param A A list of two square (weighted) adjacency matrices.
#'
#' @return The output of \code{alignMatrices()} follows the same format as its
#'   input. It produces two adjacency matrices in a list. Each matrix
#'   corresponds to the respective input entry. The input matrices have been
#'   reduced to their common row and column labels and reordered so their
#'   indeces match.
#'
#' @export
#'
#' @examples
#' # Create two symmetric matrices of different sizes with overlapping labels.
#' A1 <- matrix(runif(16), nrow=4, ncol=4); A1 <- A1 + t(A1)
#' rownames(A1) <- colnames(A1) <- LETTERS[1:4]
#' A2 <- matrix(runif(25), nrow=5, ncol=5); A2 <- A2 + t(A2)
#' rownames(A2) <- colnames(A2) <- LETTERS[2:6]
#' # Align them. This will result in two 3x3 matrices with row and
#' # column labels B, C, D.
#' alignMatrices(list(A1, A2))

alignMatrices <- function(A){
  check <- checkMatrixList(A)
  A1 <- A[[1]]; A2 <- A[[2]]; rm(A)
  if (is.null(colnames(A1)))
    commonNodes <- 1:min(ncol(A1), ncol(A2)) else
      commonNodes <- intersect(colnames(A1), colnames(A2))
  if (length(commonNodes)<2){
    warning('Fewer than two nodes in common.')
    return(NULL)
  }
  revisedA <- list(A1[commonNodes, commonNodes], A2[commonNodes, commonNodes])
  return(revisedA)
}

#' Align node metric arrays
#'
#' \code{alignArrays()} reduces two adjacency node metric arrays to their common
#' node labels and ordering them to match. The arrays can then be compared using
#' different similarity measures. Note alignment is not always necessary, since
#' many node metric functions will preserve the same node order across different
#' runs.
#'
#' @param D A list of two (named) numeric arrays.
#'
#' @return The output of \code{alignArrays()} follows the same format as its
#'   input. It produces two arrays in a list. Each array corresponds to the
#'   respective input entry. The input arrays have been reduced to their common
#'   labels and reordered so their indeces match.
#'
#' @export
#'
#' @examples
#' # Create two arrays of different sizes with overlapping labels
#' D1 <- rnorm(4)
#' names(D1) <- LETTERS[1:4]
#' D2 <- rnorm(5)
#' names(D2) <- LETTERS[2:6]
#' # Align them. This will result in two arrays of length 3 with entry
#' # labels B, C, D.
#' alignArrays(list(D1, D2))

alignArrays <- function(D){
  check <- checkArrayList(D)
  D1 <- D[[1]]; D2 <- D[[2]]; rm(D)
  if (is.null(names(D1)))
    commonNodes <- 1:min(length(D1), length(D2)) else
      commonNodes <- intersect(names(D1), names(D2))
  if (length(commonNodes)<2){
    warning('Fewer than two nodes in common.')
    return(NULL)
  }
  revisedD <- list(D1[commonNodes], D2[commonNodes])
  return(revisedD)
}

#----- Similarity and distance -----#

#' Rank k-similarity
#'
#' Rank k-similarity is the overlap between the highest ranking entries of two
#' numeric arrays. It is useful in network analysis, where two aligned networks
#' migh be considered similar if they share key nodes.
#'
#' @param D1,D2 Numeric arrays of equal length and labels, such as node metric
#'   values (e.g. two alternative degree sequences).
#' @param k.or.p Either an integer (k) or a fractiion(p). The number of highest
#'   ranking nodes to compare
#'
#' @return The rank k-similarity of two aligned vectors is the overlap of the
#'   top k highest ranking values across them. Rank k-similarity is given as a
#'   proportion.
#'
#'   For example, if the three highest ranking nodes in one network are A, B,
#'   and C, and in another they are B, C, and D, the overlap between the top
#'   three highest ranking nodes is of size 2. The rank k-similarity between the
#'   degree sequences for k=3 is therefore 2/3. If a proportion is given, e.g.
#'   \code{k.or.p=0.1}, the top \code{k.or.p*length(D1)} nodes are considered
#'   for the overlap.
#'
#' @references Bozhilova, Lyuba V., et al. "Measuring rank robustness in scored
#'   protein interaction networks." BMC Bioinformatics 20.1 (2019): 1-14.
#'
#'   Trajanovski, Stojan, et al. "Robustness envelopes of networks." Journal of
#'   Complex Networks 1.1 (2013): 44-62.
#'
#' @examples
#' # Generate two node metric vectors
#' D1 <- rpois(100, 4); D2 <- rpois(100, 3)
#' # Calculate their rank k-similarity for the top k=10 nodes.
#' calculateKSimilarity(D1, D2, k.or.p=10)
#' # Calculate their rank k-similarity for the top p=15% of nodes.
#' calculateKSimilarity(D1, D2, k.or.p=0.15)
#'
#' @export
calculateKSimilarity <- function(D1, D2, k.or.p=0.1){
  D1rank <- rank(D1, na.last=FALSE, ties.method="random")
  D2rank <- rank(D2, na.last=FALSE, ties.method ="random")
  nodeCount <- length(D1)
  if (k.or.p>0 & k.or.p<1)
    k <- round(k.or.p*nodeCount) else
      if (k.or.p>1)
        k <- min(round(k.or.p), nodeCount) else
          stop('Invalid k.')
  if (k==0)
    k==1
  ksim <- sum((D1rank>(nodeCount-k)) * (D2rank>nodeCount-k))/k
  return(ksim)
}

#' Euclidean distance
#'
#' This function computes the Euclidean distance between two numeric arrays. The
#' arrays can be rescaled prior to calculation.
#'
#' @param D1,D2 Numeric arrays of equal length and labels, such as node metric
#'   values (e.g. two alternative degree sequences).
#' @param scale Logical; Whether to scale \code{D1} and \code{D2} to [0; 1]
#'   prior to calculating the distance between them. Defaults to \code{FALSE}.
#'
#' @return The Euclidean distance between the two arrays. If \code{scale=TRUE},
#'   the arrays are first normalised using \deqn{D_scaled =
#'   (D-min(D))/(max(D)-min(D)).} The Euclidean distance is then calculated and
#'   is further divided by \code{sqrt(length(D1))} so it maps to [0; 1].
#'
#' @examples
#' # Generate two node metric vectors
#' D1 <- rpois(100, 4); D2 <- rpois(100, 3)
#' # Calculate their Euclidean distance without scaling
#' calculateEuclideanDistance(D1, D2, scale=FALSE)
#' # Calculate their Euclidean distance with scaling
#' calculateEuclideanDistance(D1, D2, scale=TRUE)
#'
#' @export
calculateEuclideanDistance <- function(D1, D2, scale=FALSE){
  nodeCount <- length(D1)
  if (scale){
    D1 <- (D1-min(D1))/(max(D1)-min(D1))
    D2 <- (D2-min(D2))/(max(D2)-min(D2))
  }
  L2dist <- sqrt(sum((D1-D2)^2))
  if (scale)
    L2dist <- L2dist/sqrt(nodeCount)
  return(L2dist)
}

#----- Data frame split -----#

#' Randomly subset expression samples
#'
#' This function takes a COGENT-compatible expression data frame and divides its
#' samples in two equally sized groups. Sample overlap between the groups is
#' allowed.
#'
#' @param df A COGENT-compatible data frame with rows corresponding to genes and
#'   columns corresponding to samples. A column called Name containing gene
#'   names is expected.
#' @param propShared The proportion of the samples in \code{df} shared across
#'   the two subsets.
#'
#' @return The function returns a list containing two COGENT-compatible data
#'   frames. Every sample of the original data frame belongs to at least one of
#'   the two data frames in the outputs. A total of propShared*numSamples are
#'   shared across the two. If the total number of samples is odd, a single
#'   sample will be omitted at random to ensure the equal sizing of the split.
#'
#' @examples
#' # Generate a COGENT-compatible data frame
#' df <- as.data.frame(matrix(runif(500), nrow=10, ncol=50))
#' df <- cbind(Name=LETTERS[1:10], df)
#' # Split with no sample overlap
#' splitExpressionData(df, propShared=0)
#' # Split so 50% of the samples are shared, and 25% each are unique to
#' # one of the two subsets.
#' splitExpressionData(df, propShared=0.50)
#'
#' @export
splitExpressionData <- function(df, propShared=0){
  check <- checkExpressionDF(df)
  sampleCount <- ncol(df)-1
  colIdx <- c(1, 1+sample.int(sampleCount))
  sharedSampleCount <- round(propShared*sampleCount)
  sharedColsIdx <- 1:(1+sharedSampleCount)
  specificSampleCount <- floor(0.5*(sampleCount-sharedSampleCount))
  df1 <- df[,colIdx[c(sharedColsIdx, (sharedSampleCount+2):(sharedSampleCount+1+specificSampleCount))]]
  df2 <- df[,colIdx[c(sharedColsIdx, (sharedSampleCount+2+specificSampleCount):(sharedSampleCount+1+2*specificSampleCount))]]
  dfList <- list(df1, df2)
  return(dfList)
}
