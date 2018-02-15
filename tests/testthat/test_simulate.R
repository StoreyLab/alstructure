# tests for simulated data generation

test_that("the function simulate_data() is working properly", {
  m <- 10; n <- 3; d <- 3; alpha <- c(0.1,0.1,0.1); seed = 1234;

  data_noseed <- simulate_admixture(m, n, d, alpha, seed = NA)
  data_1 <- simulate_admixture(m, n, d, alpha, seed = seed)
  data_2 <- simulate_admixture(m, n, d, alpha, seed = seed)

  expect_is(data_noseed, "list")
  expect_equal(data_1$F[1,1], data_2$F[1,1])
})
