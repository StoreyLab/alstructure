# tests main function

test_that("the function alstructure is working properly", {
  set.seed(12344)
  m <- 100; n <- 3; d <- 2; max_iters <- 1000; tol = 0.00001; alpha <- runif(d)
  data <- simulate_admixture(m, n, d, alpha, seed = NA)
  X <- data$X; Q <- data$Q
  factors <- alstructure(X = X, tol = tol, max_iters = max_iters, d_hat = d)

  expect_equal(names(factors), c("P_hat", "Q_hat", "rowspace", "iter","tol" ))
})

test_that("the function alstructure is working properly with trunc.svd", {
  set.seed(12344)
  m <- 100; n <- 3; d <- 2; max_iters <- 1000; tol = 0.00001; alpha <- runif(d)
  data <- simulate_admixture(m, n, d, alpha, seed = NA)
  X <- data$X; Q <- data$Q
  factors <- alstructure(X = X, tol = tol, max_iters = max_iters,
                         d_hat = d, svd_method = "truncated_svd")
})
