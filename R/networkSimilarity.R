######################################
##### COGENT: NETWORK SIMILARITY #####
######################################

#' Get the edge similarity for two networks
#'
#' \code{getEdgeSimilarity()} calculates the weighted or unweighted Jaccard
#' index between the edge sets of two networks. The calculation is performed
#' both for the full edge set and for each node neighbourhood.
#'
#' @param A A list of two square (weighted) adjacency matrices.
#' @param align Logical; Whether to align the two adjacency matrices. Only set
#'   to TRUE if this is not done automatically.
#'
#' @return The result is a list of similarity measures of the two networks with
#'   adjacency matrices in \code{A}. This includes:
#'   \itemize{
#'   \item \code{nodeCount} The number of non-isolated genes across the two
#'   networks.
#'   \item \code{globalSimilarity} The weighted Jaccard index of the edge sets
#'   of the two networks.
#'   \item \code{localSimilarity} The weighted Jaccard index for
#'   each gene neighbourhood.
#'   }
#'
#' @examples
#' # Generate two adjacency matrices
#' A1 <- matrix(0, ncol=10, nrow=10); A2 <- matrix(0, ncol=10, nrow=10)
#' A1[upper.tri(A1)] <- rbinom(45, 1, .2); A1 <- A1+t(A1)
#' A2[upper.tri(A1)] <- rbinom(45, 1, .4); A2 <- A2+t(A2)
#' colnames(A1) <- rownames(A1) <- LETTERS[1:10]
#' colnames(A2) <- rownames(A2) <- LETTERS[6:15]
#' # Calculate similarity
#' getEdgeSimilarity(list(A1, A2), align=TRUE)
#'
#' @export
getEdgeSimilarity <- function(A, align=FALSE){
  check <- checkMatrixList(A)
  if (align){
    A <- alignMatrices(A)
    if (is.null(A)){
      warning('Alignment failed.')
      return(NULL)
    }
  }
  A1 <- A[[1]]; A2 <- A[[2]]; rm(A)
  if (ncol(A1)!=ncol(A2))
    stop('Different network sizes. Consider setting align=TRUE.')
  intersectionMatrix <- pmin(A1, A2)
  intersectionDegrees <- rowSums(intersectionMatrix)
  unionMatrix <- pmax(A1, A2)
  unionDegrees <- rowSums(unionMatrix)
  isolatedNodesIdx <- which(unionDegrees==0)
  nodeLabels <- colnames(A1)
  nodeCount <- ncol(A1)
  if (length(isolatedNodesIdx)!=0){
    intersectionMatrix <- intersectionMatrix[-isolatedNodesIdx, -isolatedNodesIdx]
    intersectionDegrees <- intersectionDegrees[-isolatedNodesIdx]
    unionMatrix <- unionMatrix[-isolatedNodesIdx, -isolatedNodesIdx]
    unionDegrees <- unionDegrees[-isolatedNodesIdx]
    nodeLabels <- nodeLabels[-isolatedNodesIdx]
    nodeCount <- nodeCount-length(isolatedNodesIdx)
  }
  localSimilarity <- intersectionDegrees/unionDegrees
  names(localSimilarity) <- nodeLabels
  globalSimilarity <- sum(intersectionDegrees)/sum(unionDegrees)
  return(list(
    "nodeCount"=nodeCount,
    "globalSimilarity"=globalSimilarity,
    "localSimilarity"=localSimilarity
  ))
}

#' Get the node similarity of two networks
#'
#' \code{getNodeSimilarity()} compares two node metric vectors in one of three
#' ways: via a correlation coefficient, rank k-similarity or Euclidean distance.
#'
#' @param D  A list of two (named) numeric arrays.
#' @param mode Which measure to use for the comparison; one of "cor", "ksim", or
#'   "L2". For details, see \code{\link[stats]{cor}},
#'   \code{\link{calculateKSimilarity}} or
#'   \code{\link{calculateEuclideanDistance}} respectively.
#' @param align Logical; Whether to align the two node metric arrays. Only set
#'   to TRUE if this is not done automatically.
#' @param ... Parameters passed to the mode functions. See
#'   \code{\link[stats]{cor}}, \code{\link{calculateKSimilarity}} and
#'   \code{\link{calculateEuclideanDistance}}.
#'
#' @return The node metric similarity for the two arrays. In the case of
#'   \code{mode="cor"} this is a correlation coefficient between -1 and 1. For
#'   rank k-similarity (\code{mode="ksim"}) this is a similarity measure between
#'   0 and 1. Finally, Euclidean distance (\code{mode="L2"}) is non-negative and
#'   describes \emph{distance} rather than similarity. For correlations and rank
#'   k-similarity, higher values denote better overlap. For Euclidean distance
#'   higher values denote dissimilarity.
#'
#' @examples
#' # Generate two node metric vectors
#' D1 <- rpois(10, 4); D2 <- rpois(10, 3)
#' D <- list(D1, D2)
#' # Calculate their Spearman correlation coefficient
#' getNodeSimilarity(D, mode="cor", method="spearman")
#' # Calculate their rank k-similarity for k=3.
#' getNodeSimilarity(D, mode="ksim", k.or.p=3)
#' # Calculate their Euclidean distance without rescaling
#' getNodeSimilarity(D, mode="L2", scale=FALSE)
#'
#' @export
getNodeSimilarity <- function(D, mode=c("cor", "ksim", "L2"), align=FALSE, ...){
  check <- checkArrayList(D)
  if (align){
    D <- alignMatrices(D)
    if (is.null(D)){
      warning('Alignment failed.')
      return(NULL)
    }
  }
  D1 <- D[[1]]; D2 <- D[[2]]; rm(D)
  if (length(D1)!=length(D2))
    stop('Different network sizes. Consider setting align=TRUE.')
  mode <- mode[1]
  if (mode=="cor")
    nodeSimilarity <- cor(D1, D2, ...) else
      if (mode=="ksim")
        nodeSimilarity <- calculateKSimilarity(D1, D2, ...) else
          if (mode=="L2")
            nodeSimilarity <- calculateEuclideanDistance(D1, D2, ...) else
              stop("Mode must be one of cor, ksim or L2.")
  return(nodeSimilarity)
}

#' Get the edge similarity of two networks with configuration model correction
#'
#' Edge similarity as defined in \link{getEdgeSimilarity} is not corrected for
#' the density of the two networks. This means denser networks are more likely
#' to exhibit higher similarity than sparser ones, even when their overlap is
#' essentially random. The \code{getEdgeSimilarityCorrected()} function
#' introduces a density correction by considering the random overlap of networks
#' generated by the confiduration model with degree sequences matching the
#' input.
#'
#' @inheritParams getEdgeSimilarity
#' @param A  A list of two square \strong{unweighted} adjacency matrices.
#' @param type How to do the correction: either by generating random networks
#'   from the configuartion model or by using expectations. Must be one of
#'   "random" or "expected".
#'
#' @return The result is a list of two elements:
#'   \itemize{
#'   \item \code{nodeCount} The number of non-isolated genes across the two
#'   networks.
#'   \item \code{correctedSimilarity} The weighted Jaccard index of the edge
#'   sets of the two networks, corrected using the generated configuration model
#'   networks. } The corrected similarity corresponds to the global similarity
#'   in \link{getEdgeSimilarity}. The correction is done by generating two
#'   networks using the permuted degree sequence of the networks in \code{A}.
#'   These random networks are then overlapped to give a measure of the possible
#'   random overlap. The global similarity is corrected by removing the
#'   random overlap from the numerator and addding it to the denominator of the
#'   corresponging Jaccard index.
#'
#' @examples
#' # Generate two adjacency matrices
#' A1 <- matrix(0, ncol=10, nrow=10); A2 <- matrix(0, ncol=10, nrow=10)
#' A1[upper.tri(A1)] <- rbinom(45, 1, .2); A1 <- A1+t(A1)
#' A2[upper.tri(A1)] <- rbinom(45, 1, .4); A2 <- A2+t(A2)
#' colnames(A1) <- rownames(A1) <- LETTERS[1:10]
#' colnames(A2) <- rownames(A2) <- LETTERS[6:15]
#' # Calculate their corrected global similarity using random networks
#' getEdgeSimilarityCorrected(list(A1, A2), align=TRUE, type="random")
#' # Calculate their corrected global similarity using expectations
#' getEdgeSimilarityCorrected(list(A1, A2), align=TRUE, type="expected")
#'
#' @export
getEdgeSimilarityCorrected <- function(A, align=FALSE, type=c("random", "expected")){
  check <- checkMatrixList(A)
  if (align){
    A <- alignMatrices(A)
    if (is.null(A)){
      warning('Alignment failed.')
      return(NULL)
    }
  }
  A1 <- A[[1]]; A2 <- A[[2]]; rm(A)
  if (ncol(A1)!=ncol(A2))
    stop('Different network sizes. Consider setting align=TRUE.')
  intersectionMatrix <- A1*A2
  intersectionDegrees <- rowSums(intersectionMatrix)
  unionMatrix <- 1*((A1+A2)>0)
  unionDegrees <- rowSums(unionMatrix)
  isolatedNodesIdx <- which(unionDegrees==0)
  nodeLabels <- colnames(A1)
  nodeCount <- ncol(A1)
  if (length(isolatedNodesIdx)!=0){
    intersectionMatrix <- intersectionMatrix[-isolatedNodesIdx, -isolatedNodesIdx]
    intersectionDegrees <- intersectionDegrees[-isolatedNodesIdx]
    unionMatrix <- unionMatrix[-isolatedNodesIdx, -isolatedNodesIdx]
    unionDegrees <- unionDegrees[-isolatedNodesIdx]
    nodeLabels <- nodeLabels[-isolatedNodesIdx]
    nodeCount <- nodeCount-length(isolatedNodesIdx)
    A1 <- A1[-isolatedNodesIdx, -isolatedNodesIdx]
    A2 <- A2[-isolatedNodesIdx, -isolatedNodesIdx]
  }
  if (type=="random"){
    R1 <- as_adjacency_matrix(sample_degseq(sample(rowSums(A1)), method="simple.no.multiple"))
    R2 <- as_adjacency_matrix(sample_degseq(sample(rowSums(A2)), method="simple.no.multiple"))
    correctionTerm <- sum(A1*R2 + A2*R1)*0.5
  } else if (type=="expected"){
    if (all(rowSums(A1)==0) | all(rowSums(A2)==0)){
      correctionTerm <- 0
    } else {
      eR1 <- sample(rowSums(A1)); eR2 <- sample(rowSums(A2))
      eR1 <- eR1 %*% t(eR1); diag(eR1) <- 0; eR1 <- eR1*sum(A1)/sum(eR1)
      eR2 <- eR2 %*% t(eR2); diag(eR2) <- 0; eR2 <- eR2*sum(A2)/sum(eR2)
      correctionTerm <- sum(A1*eR2 + A2*eR1)*0.5
    }
  } else stop("Type should be one of \"random\" or \"expected\"")
  correctedSimilarity <- (sum(intersectionDegrees)-correctionTerm)/(sum(unionDegrees)+correctionTerm)
  return(list(
    "nodeCount"=nodeCount,
    "correctedSimilarity"=correctedSimilarity
  ))
}
