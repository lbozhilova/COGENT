test_that("checkMatrixList works", {
  expressionData <- matrix(runif(25*20), nrow=25, ncol=20)
  A1 <- cor(expressionData[,1:10]); A2 <- cor(expressionData[,11:20])
  expect_equal(checkMatrixList(list(A1, A2)), TRUE)
})

test_that("checkArrayList works", {
  D1 <- runif(10); names(D1) <- LETTERS[1:10]
  D2 <- runif(10); names(D2) <- LETTERS[5:14]
  expect_equal(checkArrayList(list(D1, D2)), TRUE)
})

test_that("checkExpressionDF works", {
  df <- as.data.frame(matrix(runif(500)), nrow=10, ncol=50)
  df <- cbind(Name=LETTERS[1:10], df)
  expect_equal(checkExpressionDF(df), TRUE)
})

test_that("checkFun works", {
  foo <- function(x) x+1
  expect_equal(checkFun(foo), TRUE)
  expect_equal(checkFun("foo"), TRUE)
})
