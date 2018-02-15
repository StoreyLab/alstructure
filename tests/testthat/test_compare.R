# tests for comparison functions

test_that("the make_bed() function is working properly",{
  m <- 101; n <- 24; d <- 3; alpha <- c(0.1, 0.1, 0.1); seed <- 1234

  data <- simulate_admixture(m, n, d, alpha, seed = seed)

  out_file <- "/Users/irineocabreros/Desktop/Research/Storey/cabreros/F=PQ/simulations/sample_beds/Y_sample"

  make_bed(data$X, out_file, B = 1)
})

test_that("the order_Q() frunction is working properly with var_explained method", {
  m <- 1000; n <- 100; d <- 10; real_inds <- c(1, 6, 8); alpha <- rep(0.1, d)

  P <- matrix(0, nrow = m, ncol = d)
  P[, real_inds] <- 1
  Q <- t(alstructure:::rdirichlet(n, alpha))
  F <- P %*% Q
  X <- matrix(rbinom(n * m, 2, c(F)), ncol = n, nrow = m)
  Q_space <- lse(X, d)

  Q_order <- order_pops(P = P, Q = Q, method = "var_explained", Q_space = Q_space)

  guess_inds <- c(which(Q_order$perm_mat[1, ] == 1),
    which(Q_order$perm_mat[2, ] == 1),
    which(Q_order$perm_mat[3, ] == 1))

  expect_equal(length(setdiff(real_inds, guess_inds)),  0)
})

test_that("the order_Q() frunction throws correct error", {
  m <- 1000; n <- 100; d <- 10; real_inds <- c(1, 6, 8); alpha <- rep(0.1, d)

  P <- matrix(0, nrow = m, ncol = d)
  P[, real_inds] <- 1
  Q <- t(alstructure:::rdirichlet(n, alpha))
  F <- P %*% Q
  X <- matrix(rbinom(n * m, 2, c(F)), ncol = n, nrow = m)
  Q_space <- lse(X, d)

  expect_error(order_pops(P = P, Q = Q, method = "var_explained"), "Q_space must be supplied for var_explained method")

})

test_that("the order_Q() frunction works properly with ave_admixture method", {
  m <- 1000; n <- 100; d <- 3; real_inds <- c(2, 1, 3); alpha <- rep(0.1, d)

  P <- matrix(0, nrow = m, ncol = d)
  P[, real_inds] <- 1
  Q <- matrix(0, nrow = d, ncol = n)
  Q[1, ] <- rep(0.25, n)
  Q[2, ] <- rep(0.75, n)
  Q[3, ] <- rep(0, n)

  Q_order <- order_pops(P = P, Q = Q)

  guess_inds <- c(which(Q_order$perm_mat[1, ] == 1),
    which(Q_order$perm_mat[2, ] == 1),
    which(Q_order$perm_mat[3, ] == 1))

  expect_equal(guess_inds, real_inds)

})
