# tests for data visualization functions

test_that("the function compare_q() is working properly", {

    m <- 200; n <- 10; d <- 3; alpha <- c(0.1,0.1,0.1); seed = 1234
    data <- simulate_admixture(m, n, d, alpha, seed = seed)

    output_dir <- "/Users/irineocabreros/Desktop/Research/Storey/cabreros/F=PQ/simulations/sample_plots/example.png"
    output_dir_held_out <- "/Users/irineocabreros/Desktop/Research/Storey/cabreros/F=PQ/simulations/sample_plots/example_holdouts.png"

    ########################
    # these tests are slow #
    ########################

    #compare_q(data$Q, data$Q, output_file = output_dir, dims = c(1,2), res = 500, rand_held_out = 0)

    #compare_q(data$Q, data$Q, output_file = output_dir_held_out, dims = c(1,2), res = 100, rand_held_out = 5)
})

