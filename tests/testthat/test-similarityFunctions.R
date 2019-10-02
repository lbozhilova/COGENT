test_that("getEdgeSimilarity works", {
  A1 <- matrix(0, ncol=10, nrow=10); A2 <- matrix(0, ncol=10, nrow=10)
  A1[upper.tri(A1)] <- rbinom(45, 1, .2); A1 <- A1+t(A1)
  A2[upper.tri(A1)] <- rbinom(45, 1, .4); A2 <- A2+t(A2)
  colnames(A1) <- rownames(A1) <- LETTERS[1:10]
  colnames(A2) <- rownames(A2) <- LETTERS[6:15]
  es <- getEdgeSimilarity(list(A1, A2), align=TRUE)
  expect_equal(unname(lengths(es)[3]<=5), TRUE)
  expect_equal(length(es), 3)
})

test_that("getNodeSimilarity works", {
  D1 <- rpois(10, 4); D2 <- rpois(10, 3)
  D <- list(D1, D2)
  ns <- getNodeSimilarity(D, mode="cor", method="spearman")
  expect_equal((ns>=-1)&(ns<=1), TRUE)
  ns <- getNodeSimilarity(D, mode="ksim", k.or.p=3)
  expect_equal((ns>=0)&(ns<=1), TRUE)
  ns <- getNodeSimilarity(D, mode="L2", scale=FALSE)
  expect_equal(ns>=0, TRUE)
  ns <- getNodeSimilarity(D, mode="L2", scale=TRUE)
  expect_equal((ns>=0)&(ns<=1), TRUE)
})

test_that("getEdgeSimilarityCorrected works", {
  A1 <- matrix(0, ncol=10, nrow=10); A2 <- matrix(0, ncol=10, nrow=10)
  A1[upper.tri(A1)] <- rbinom(45, 1, .2); A1 <- A1+t(A1)
  A2[upper.tri(A1)] <- rbinom(45, 1, .4); A2 <- A2+t(A2)
  colnames(A1) <- rownames(A1) <- LETTERS[1:10]
  colnames(A2) <- rownames(A2) <- LETTERS[6:15]
  es <- unname(getEdgeSimilarityCorrected(list(A1, A2), align=TRUE, type="random")[[2]])
  expect_equal((es>=-1)&(es<=1), TRUE)
  es2 <- unname(getEdgeSimilarityCorrected(list(A1, A2), align=TRUE, type="expected")[[2]])
  expect_equal((es2>=-1)&(es2<=1), TRUE)
})
