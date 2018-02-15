# tests for factorization algorithms

test_that("the function projsplx() is working properly", {
  v <- c(2, 0, 0)
  v_proj <- projsplx(v)
  v_proj_true <- c(1, 0, 0)

  expect_equal(v_proj, v_proj_true)

  v <- c(10, 10, 10)
  v_proj <- projsplx(v)
  v_proj_true <- c(1/3, 1/3, 1/3)

  expect_equal(v_proj, v_proj_true)
})

test_that("the function factor_F is properly functioning", {
  set.seed(12345)
  m <- 10; n <- 3; k <- 2; max_iters <- 1000; tol = 0.0001; alpha <- runif(k)
  P <- matrix(runif(m * k), nrow = m, ncol = k)
  Q <- t(gtools::rdirichlet(n, alpha))
  F <- P %*% Q

  hats <- factor_F(F, k, tol = tol, max_iters = max_iters)
  P_hat <- hats$P_hat; Q_hat <- hats$Q_hat; F_hat <- P_hat %*% Q_hat

  expect_equal(F_hat, F, tolerance = 1e-3)
  expect_equal(sum(Q_hat), n)
  expect_equal(sum(P > 1), 0)
  expect_equal(sum(P < 0), 0)
})

test_that("the function best_perm() is working properly", {
  Q <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
  Q_hat <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

  perm_mat <- best_perm(Q, Q_hat)
  perm_mat_true <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)

  expect_equal(perm_mat, perm_mat_true)
  expect_is(perm_mat, "matrix")
})


test_that("the impute_mean function is working properly", {
  X <- matrix(c(0,2, NA, 2, NA, NA, NA, 6, NA), nrow = 3, ncol = 3, byrow = TRUE)
  X_impute_mean <- impute_mean(X)
  X_true <- matrix(c(0,2, 1, 2, 2, 2, 6, 6, 6), nrow = 3, ncol = 3, byrow = TRUE)
  expect_equal(X_impute_mean, X_true)
})
