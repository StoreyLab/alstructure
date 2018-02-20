# tests for estimation of Fhat

test_that("the function D_binomial() is working properly", {
  vals_Y <- c(1,1,0,1,2,1,2,1)
  Y <- matrix(vals_Y, nrow = 4, ncol = 2)

  vals_D <- c(0.75, 0, 0, 0.5)
  D_true <- matrix(vals_D, nrow = 2, ncol = 2)

  D <- D_binomial(Y)

  expect_is(D, "matrix")
  expect_equal(D, D_true)
})

test_that("the function estimate_F() is working properly", {
  Y <- matrix(rbinom(100, 2, 0.2), nrow = 20, ncol = 5)
  F_obj <- estimate_F(Y, d = 3)

  expect_is(F_obj$F_hat, "matrix")
  expect_equal(sum(F_obj$F_hat >= 0), 100)
  expect_equal(sum(F_obj$F_hat <= 1), 100)

  Y <- matrix(c(2,0,0,0,2,2,2,0,0,0,2,2,0,2,2,2,0,0,0,2,2,2,0,0),
              nrow = 6, ncol = 4)
  Fhat_true <- matrix(c(1,0,0,0,1,1,1,0,0,0,1,1,0,1,1,1,0,0,0,1,1,1,0,0),
                      nrow = 6, ncol = 4)
  F_obj <- estimate_F(Y, d = 2)

  expect_equal(F_obj$F_hat, Fhat_true)
})

test_that("the function estimate_F() is working properly with truncated_svd method", {
  Y <- matrix(rbinom(100, 2, 0.2), nrow = 20, ncol = 5)
  F_obj <- estimate_F(Y, d = 3, svd_method = "truncated_svd")

  expect_is(F_obj$F_hat, "matrix")
  expect_equal(sum(F_obj$F_hat >= 0), 100)
  expect_equal(sum(F_obj$F_hat <= 1), 100)

  Y <- matrix(c(2,0,0,0,2,2,2,0,0,0,2,2,0,2,2,2,0,0,0,2,2,2,0,0),
    nrow = 6, ncol = 4)
  Fhat_true <- matrix(c(1,0,0,0,1,1,1,0,0,0,1,1,0,1,1,1,0,0,0,1,1,1,0,0),
    nrow = 6, ncol = 4)
  F_obj <- estimate_F(Y, d = 2)

  expect_equal(F_obj$F_hat, Fhat_true)
})

test_that("the function estimate_d() is working properly", {

  m <- 1000; n <- 200; d <- 2; alpha <- c(0.1, 0.1); size = 2; seed = 1234
  data <- simulate_admixture(m, n, d, alpha, seed = seed)

  d_hat <- estimate_d(data$X)

  #expect_equal(k_hat, k)
})

test_that("the function lse() is working properly with truncated_svd method", {
  m <- 1000; n <- 200; d <- 2; alpha <- c(0.1, 0.1); size = 2; seed = 1234
  data <- simulate_admixture(m, n, d, alpha, seed = seed)

  rowspace_vanilla <- lse(data$X, d = d)
  rowspace_trunc <- lse(data$X, d = d, svd_method = "truncated_svd")

  expect_equal(dim(rowspace_vanilla$vectors), c(n, d))
  expect_equal(dim(rowspace_trunc$vectors), c(n, d))
  expect_equal(length(rowspace_vanilla$values), 2)
  expect_equal(length(rowspace_trunc$values), 2)

})

test_that("the function lse() is working properly with truncated_svd method", {
  m <- 1000; n <- 200; d <- 2; alpha <- c(0.1, 0.1); size = 2; seed = 1234
  data <- simulate_admixture(m, n, d, alpha, seed = seed)

  rowspace_vanilla <- lse(data$X, d = d)
  rowspace_trunc <- lse(data$X, d = d, svd_method = "base")

  expect_equal(dim(rowspace_vanilla$vectors), c(n, d))
  expect_equal(dim(rowspace_trunc$vectors), c(n, d))
  expect_equal(length(rowspace_vanilla$values), 2)
  expect_equal(length(rowspace_trunc$values), 2)

})
