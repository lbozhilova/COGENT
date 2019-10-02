test_that("cogentSingle works", {
  df <- as.data.frame(matrix(runif(500), nrow=10, ncol=50))
  df <- cbind(Name=LETTERS[1:10], df)
  foo <- function(df){
    A <- cor(t(df[,colnames(df)!="Name"]))
    A <- 1*(A>0.20)
    return(A)
  }
  x <- cogentSingle(df, foo)
  expect_equal(class(x), "list")
  expect_equal(length(x), 3)
  fooNode <- function(x) rowSums(x)
  x <- cogentSingle(df, foo, fooNode, nodeModes=c("cor", "ksim"), method="spearman")
  expect_equal(length(x), 5)
  x <- cogentSingle(df, foo, fooNode, nodeModes=c("all"), method="spearman")
  expect_equal(length(x), 6)
})

test_that("cogentLinear works", {
  df <- as.data.frame(matrix(runif(500), nrow=10, ncol=50))
  df <- cbind(Name=LETTERS[1:10], df)
  foo <- function(df){
    A <- cor(t(df[,colnames(df)!="Name"]))
    A <- 1*(A>0.20)
    return(A)
  }
  fooNode <- function(x) rowSums(x)
  x <- cogentLinear(df, foo, repCount=10)
  expect_equal(ncol(x)==3 & class(x)=="data.frame", TRUE)
  x <- cogentLinear(df, foo, fooNode, repCount=10)
  expect_equal(ncol(x)==6 & class(x)=="data.frame", TRUE)
})

test_that("cogentParallel works", {
  df <- as.data.frame(matrix(runif(500), nrow=10, ncol=50))
  df <- cbind(Name=LETTERS[1:10], df)
  foo <- function(df){
    A <- cor(t(df[,colnames(df)!="Name"]))
    A <- 1*(A>0.20)
    return(A)
  }
  fooNode <- function(x) rowSums(x)
  x <- cogentParallel(df, foo, repCount=10, threadCount=2)
  expect_equal(ncol(x)==3 & class(x)=="data.frame", TRUE)
  x <- cogentParallel(df, foo, fooNode, repCount=10, threadCount=2)
  expect_equal(ncol(x)==6 & class(x)=="data.frame", TRUE)
})
