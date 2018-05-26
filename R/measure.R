# Functions below are used to measure the goodness of fit of the estimates.

#' Compute the mean absolute value difference between matrices \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{B}
#'
#' Compute the mean absolute difference value between matrices \eqn{\boldsymbol{A}}{A}
#'  and \eqn{\boldsymbol{B}}{B}. If \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{B}
#'  are \eqn{m \times n}{m x n} matrices, then the mean absolute value difference
#' is:
#' \deqn{\frac{1}{mn} \sum_{ij} |a_{ij} - b_{ij}|}{1 / mn sum(|a_ij - b_ij|)}
#' This function is useful for comparing the performance of \code{ALStructure}
#' fits of global ancestry to either alternative methods or the ground truth
#' (in the case of simulations).
#'
#' @param A An \eqn{m \times n}{m x n} matrix
#' @param B An \eqn{m \times n}{m x n} matrix
#'
#' @return The mean absolute value difference between \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{B}.
#'
#' @keywords internal
mean_ABS <- function(A, B){
  mean_error <- mean(abs(A-B))
}

#' Compute the root mean squared error (RMSE) between matrices \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{B}.
#'
#' Compute the root mean squared error (RMSE) between matrices \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{A}. If \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{B} are \eqn{m \times n}{m x n} matrices, then the RMSE is:
#' \deqn{\sqrt{\frac{1}{mn} \sum_{ij} |a_{ij} - b_{ij}|^2}}{sqrt(1 / mn sum((a_ij - b_ij)^2))}
#' This function is useful for comparing the performance of \code{ALStructure}
#' fits of global ancestry to either alternative methods or the ground truth
#' (in the case of simulations). This is also the function used in the convergence
#' criterion of \code{cALS} and \code{uALS}.
#'
#' @param A An \eqn{m \times n}{m x n} matrix
#' @param B An \eqn{m \times n}{m x n} matrix
#'
#' @return The RMSE between \eqn{\boldsymbol{A}}{A} and \eqn{\boldsymbol{B}}{B}.
#'
#' @keywords internal
RMSE <- function(A, B){
  RMSE <- sqrt(mean((A-B)^2))
}

#' Compute the binomial likelihood of the SNP matrix given individual-specific
#' allele frequencies.
#'
#' Compute the binomial likelihood of the data \eqn{\boldsymbol{X}}{X} given the
#' parameters \eqn{\boldsymbol{F}}{F}. The likelihood
#' is the probability of the data given the parameters, which under the admixture
#' model is:
#' \deqn{\prod_{ij} \binom{2}{x_{ij}} f_{ij}^{x_{ij}} (1 - f_{ij})^{2 - x_{ij}}}
#'
#' @param X The \eqn{m \times n}{m x n} SNP data matrix
#' @param F The \eqn{m \times n}{m x n} matrix of binomial parameters.
#'
#' @return The likelihood and log likelihood of the data given the parameters
#'
#' @keywords internal
binomial_likelihood <- function(X, F){
  m <- dim(X)[1]

  X1 <- X[1:floor(m / 3), ]
  X2 <- X[(floor(m / 3) + 1):(floor(2 * m / 3)), ]
  X3 <- X[(floor(2 * m / 3) + 1):m, ]
  X1_vec <- c(X1[!is.na(X1)])
  X2_vec <- c(X2[!is.na(X2)])
  X3_vec <- c(X3[!is.na(X3)])

  F1 <- c(F[1:floor(m / 3), ])
  F2 <- c(F[(floor(m / 3) + 1):(floor(2 * m / 3)), ])
  F3 <- c(F[(floor(2 * m / 3) + 1):m, ])
  F1_vec <- c(F1[!is.na(X1)])
  F2_vec <- c(F2[!is.na(X2)])
  F3_vec <- c(F3[!is.na(X3)])

  L1 <- dbinom(X1_vec, 2, F1_vec)
  L2 <- dbinom(X2_vec, 2, F2_vec)
  L3 <- dbinom(X3_vec, 2, F3_vec)
  L <- prod(L1) * prod(L2) * prod(L3)
  l <- sum(log(L1)) + sum(log(L2)) + sum(log(L3))
  vals <- list(L = L, l = l)
}




