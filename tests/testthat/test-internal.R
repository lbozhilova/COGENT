test_that("alignMatrices works", {
  A1 <- matrix(runif(16), nrow=4, ncol=4); A1 <- A1 + t(A1)
  rownames(A1) <- colnames(A1) <- LETTERS[1:4]
  A2 <- matrix(runif(25), nrow=5, ncol=5); A2 <- A2 + t(A2)
  rownames(A2) <- colnames(A2) <- LETTERS[2:6]
  A <- alignMatrices(list(A1, A2))
  expect_equal(checkMatrixList(A), TRUE)
  expect_equal(colnames(A[[1]]), colnames(A[[2]]))
  expect_equal(dim(A[[1]]), c(3, 3))
})

test_that("alignArrays works", {
  D1 <- rnorm(4)
  names(D1) <- LETTERS[1:4]
  D2 <- rnorm(5)
  names(D2) <- LETTERS[2:6]
  D <- alignArrays(list(D1, D2))
  expect_equal(checkArrayList(D), TRUE)
  expect_equal(names(D[[1]]), names(D[[2]]))
})

test_that("calculateKSimilarity works", {
  D1 <- rpois(100, 4); D2 <- rpois(100, 3)
  ksim <- calculateKSimilarity(D1, D2, k.or.p=10)
  expect_equal((ksim>=0)&(ksim<=1), TRUE)
  ksim <- calculateKSimilarity(D1, D2, k.or.p=0.15)
  expect_equal((ksim>=0)&(ksim<=1), TRUE)
})

test_that("calculateEuclideanDistance", {
  D1 <- rpois(100, 4); D2 <- rpois(100, 3)
  lds <- calculateEuclideanDistance(D1, D2, scale=FALSE)
  expect_equal(lds>=0, TRUE)
  lds <- calculateEuclideanDistance(D1, D2, scale=TRUE)
  expect_equal((lds>=0)&(lds<=1), TRUE)
})

test_that("splitExpressionData works", {
  df <- as.data.frame(matrix(runif(500), nrow=10, ncol=50))
  df <- cbind(Name=LETTERS[1:10], df)
  dfL <- splitExpressionData(df, propShared=0)
  expect_equal(ncol(dfL[[1]]), ncol(dfL[[2]]))
  expect_equal(nrow(dfL[[1]]), nrow(dfL[[2]]))
  expect_equal(ncol(dfL[[1]]), 26)
  expect_equal(nrow(dfL[[1]]), 10)
  dfL <- splitExpressionData(df, propShared=0.50)
  expect_equal(ncol(dfL[[1]]), ncol(dfL[[2]]))
  expect_equal(nrow(dfL[[1]]), nrow(dfL[[2]]))
  expect_equal(ncol(dfL[[1]]), 38)
  expect_equal(nrow(dfL[[1]]), 10)
})
