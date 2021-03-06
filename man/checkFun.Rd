% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inputCheck.R
\name{checkFun}
\alias{checkFun}
\title{Check function}
\usage{
checkFun(fun)
}
\arguments{
\item{fun}{A function or a function name.}
}
\value{
If the input passes the check, the output is \code{TRUE}. Otherwise,
  an error is thrown.
}
\description{
\code{checkFun()} checks whether its input is a function or a function name.
It is part of the input check for main COGENT functions.
}
\examples{
foo <- function(x) x+1
checkFun(foo)
checkFun("foo")

}
