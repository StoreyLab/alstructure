# tests for error measurement computations

test_that("the function mean_ABS() is working properly",{
  A <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
  B <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)

  mean_error_true <- 1
  mean_error <- mean_ABS(A,B)

  expect_equal(mean_error, mean_error_true)
})

test_that("the function RMSE() is working properly",{
  A <- matrix(c(2,0,0,2), nrow = 2, ncol = 2)
  B <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)

  RMSE_true <- sqrt(10/4)
  RMSE_compute <- RMSE(A,B)

  expect_equal(RMSE_compute, RMSE_true)
})

test_that("the function binomial_likelihood() is working properly",{
  X <- matrix(c(0,0,0,0,NA,NA), nrow = 2, ncol = 3)
  F <- matrix(c(0.5,0.5,0.5,0.5, 1, 1), nrow = 2, ncol = 3)
  likelihood <- binomial_likelihood(X, F)

  expect_equal(likelihood$L, .25^4)
})
