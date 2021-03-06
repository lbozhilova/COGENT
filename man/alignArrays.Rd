% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{alignArrays}
\alias{alignArrays}
\title{Align node metric arrays}
\usage{
alignArrays(D)
}
\arguments{
\item{D}{A list of two (named) numeric arrays.}
}
\value{
The output of \code{alignArrays()} follows the same format as its
  input. It produces two arrays in a list. Each array corresponds to the
  respective input entry. The input arrays have been reduced to their common
  labels and reordered so their indeces match.
}
\description{
\code{alignArrays()} reduces two adjacency node metric arrays to their common
node labels and ordering them to match. The arrays can then be compared using
different similarity measures. Note alignment is not always necessary, since
many node metric functions will preserve the same node order across different
runs.
}
\examples{
# Create two arrays of different sizes with overlapping labels
D1 <- rnorm(4)
names(D1) <- LETTERS[1:4]
D2 <- rnorm(5)
names(D2) <- LETTERS[2:6]
# Align them. This will result in two arrays of length 3 with entry
# labels B, C, D.
alignArrays(list(D1, D2))
}
