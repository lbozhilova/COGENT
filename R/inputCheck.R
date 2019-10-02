###############################
##### COGENT: INPUT CHECK ######
###############################

# Check correct input for the main parameters of COGENT functions.

#' Check adjacency matrix list
#'
#' \code{checkMatrixList()} checks whether its input is a list of matrices of
#' length two. It does \strong{not} check whether the matrices are symmetric.
#'
#' @param A A list of two square (weighted) adjacency matrices.
#'
#' @return If the input passes the check, the output is \code{TRUE}. Otherwise,
#'   one of three errors is thrown, corresponding to each of the three checks
#'   performed by \code{checkMatrixList()}. These checks are:
#'   \itemize{
#'   \item Is \code{A} a list of length 2?
#'   \item Are the matrices it contains square?
#'   \item Is there a label mismatch? Column names should be either present in
#'   both or absent in both matrices.
#'   }
#'
#' @examples
#' expressionData <- matrix(runif(25*20), nrow=25, ncol=20)
#' A1 <- cor(expressionData[,1:10]); A2 <- cor(expressionData[,11:20])
#' checkMatrixList(list(A1, A2))
#'
#' @export
checkMatrixList <- function(A){
  if (class(A)!="list" | length(A)!=2)
    stop("Input should be a list of length 2.")
  A1 <- A[[1]]; A2 <- A[[2]]; rm(A)
  if (nrow(A1)!=ncol(A1)| nrow(A2)!=ncol(A2))
    stop('Adjacency matrix not square.')
  if(is.null(colnames(A1))!=is.null(colnames(A2)))
    stop('Label mismatch: one matrix missing node labels.')
  return(TRUE)
}

#' Check node metric list
#'
#' \code{checkArrayList()} checks whether its input is a list of arrays.
#'
#' @param D A list of two numeric arrays, corresponding to node metrics.
#'
#' @return If the input passes the check, the output is \code{TRUE}. Otherwise,
#'   one of two errors is thrown, corresponding to each of the two checks
#'   performed by \code{checkArrayList()}. These checks are:
#'   \itemize{
#'   \item Is \code{D} a list of length 2?
#'   \item Is there a label mismatch? Either both
#'   arrays or neither array should be named.
#'   }
#'
#' @examples
#' D1 <- runif(10); names(D1) <- LETTERS[1:10]
#' D2 <- runif(10); names(D2) <- LETTERS[5:14]
#' checkArrayList(list(D1, D2))
#'
#' @export
checkArrayList <- function(D){
  if (class(D)!="list" | length(D)!=2)
    stop("Input should be a list of length 2.")
  D1 <- D[[1]]; D2 <- D[[2]]; rm(D)
  if(is.null(names(D1))!=is.null(names(D2)))
    stop('Label mismatch: one array missing node labels.')
  return(TRUE)
}

#' Check a gene expression data frame
#'
#' \code{checkExpressionDF()} checks whether its input is a COGENT-compatible
#' gene expressiom data frame.
#'
#' @param df A COGENT-compatible gene expression data frame.
#'
#' @return If the input passes the check, the output is \code{TRUE}. Otherwise,
#'   one of three errors is thrown, corresponding to each of the three checks
#'   performed by \code{checkExpressionDF()}. These checks are:
#'   \itemize{
#'   \item Is \code{df} a data frame?
#'   \item Does it have a column \code{Name} where gene names are recorded?
#'   \item Are all other columns of type numeric?
#'   }
#'   Note the \code{Name} column is not required to be of type \code{character}.
#'   Numeric gene indeces, \code{NA}s and \code{NULL}s are allowed.
#'
#' @examples
#' df <- as.data.frame(matrix(runif(500)), nrow=10, ncol=50)
#' df <- cbind(Name=LETTERS[1:10], df)
#' checkExpressionDF(df)
#'
#' @export
checkExpressionDF <- function(df){
  if (class(df)!="data.frame")
    stop("Input should be a data frame.")
  if (!("Name" %in% colnames(df)))
    stop("Gene names should be recorded in a Name column.")
  if (any(sapply(df[,colnames(df)!="Name"], class)!="numeric"))
    stop("Gene expression should be numeric.")
  return(TRUE)
}
#' Check function
#'
#' \code{checkFun()} checks whether its input is a function or a function name.
#' It is part of the input check for main COGENT functions.
#'
#' @param fun A function or a function name.
#'
#' @return If the input passes the check, the output is \code{TRUE}. Otherwise,
#'   an error is thrown.
#'
#' @examples
#' foo <- function(x) x+1
#' checkFun(foo)
#' checkFun("foo")
#'
#' @export
checkFun <- function(fun){
  if (class(fun) %in% c("function", "character", "string"))
    return(TRUE)
  stop("Input should be a function or a function name.")
}
