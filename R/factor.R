# Functions below are used to factor the F_hat matrix using the ALS algorithms.

#' Project a vector onto the simplex
#'
#' Algorithm from (Y. Chen and Ye 2011) to project an \eqn{n}{n} dimensional vector onto the \eqn{n}{n} dimensional
#' simplex (i.e. the set of points in \eqn{x \in \mathcal{R}^n}{x in R^n} such that
#' \eqn{\sum_i x_i = 1}{x_1 + x_2 + ... + x_n= 1} and \eqn{x_i > 0}{x_i > 0} for all \eqn{i}{i} ). This is used to enforce the admixture constraints in
#' the \code{uALS} algorithm.
#'
#' @param y a \eqn{n}{n} dimensional vector
#'
#' @return a \eqn{n}{n} dimensional vector which is the projection of \eqn{y}{y} onto
#'         the \eqn{n}{n} dimensional simplex
#'
#' @references
#' Chen, Y., and X. Ye. 2011. “Projection Onto A Simplex.” ArXiv E-Prints, January.
#'
#' @keywords internal
projsplx <- function(y){
  m = length(y)
  bget = FALSE
  s <- sort(y, decreasing = TRUE)
  tmpsum = 0

  for (ii in 1:(m - 1)){
    tmpsum = tmpsum + s[ii]
    tmax = (tmpsum - 1) / ii
    if (tmax >= s[ii + 1]){
      bget = TRUE
      break
    }
  }

  if (!bget){
    tmax = (tmpsum + s[m] - 1)/m
  }

  x <- pmax(y - tmax, rep(0, m));
}

#' Impute an SNP matrix with NAs
#'
#' Impute an SNP matrix with NAs by replacing each NA with the row mean
#' (excluding NAs) from each row
#'
#' @param an \eqn{m \times n}{m x n} SNP matrix with NAs
#'
#' @return an \eqn{m \times n}{m x n} SNP matrix with mean-imputed values.
#'
#' @keywords internal
impute_mean <- function(X){
  # compute the number of
  m <- dim(X)[1]
  mu_row <- apply(X, 1, function(x) mean(x, na.rm = TRUE))
  for(i in 1:m){
    X[i, which(is.na(X[i, ]))] <- mu_row[i]
  }
  X
}

#' A fast algorithm for factoring \eqn{\boldsymbol{\hat{F}}}{F_hat}.
#'
#' An algorithm for finding approximate factors \eqn{\boldsymbol{\hat{P}}}{P_hat} and \eqn{\boldsymbol{\hat{Q}}}{Q_hat} of the matrix \eqn{\boldsymbol{\hat{F}}}{F_hat} that obey constraints of the admixture model:
#' \enumerate{
#'  \item \eqn{p_{ij} \in [0,1]  ~ \forall (i,j)}{p_ij in [0, 1] for all (i, j)}
#'  \item \eqn{q_{ij} \ge 0  ~ \forall (i,j)}{q_ij >= 0 for all (i, j)} and \eqn{\sum_i q_{ij} = 1 ~ \forall j}{q_1j + q_2j + ... + q_dj = 1 for all j}
#' }
#' This algorithm is described in Algorithm 2 in (Cabreros and Storey 2017). While
#' it lacks the theoretical guarantees of the \code{cALS} algorithm, it is much faster.
#'
#' @param F_hat The estimate of the matrix \eqn{\boldsymbol{F}}{F} to be factored
#' @param d The dimension of the latent space. This can be estimated by the
#'        function \code{d_estimate}.
#' @param tol The convergence criterion. If \eqn{RMSE(\boldsymbol{\hat{Q}}_t - \boldsymbol{\hat{Q}}_{t + 1})
#'        < tol}{RMSE(Q_t - Q_(t + 1)) < tol}, then the algorithm halts
#' @param max_iters The maximum number of iterations (repetitions of steps (6)
#'        and (7) in Algorithm 1) to be executed
#'
#' @return A list with the following elements:
#'    \describe{
#'    \item{P_hat}{: The estimated \eqn{\boldsymbol{P}}{P} matrix. Each
#'    column of \eqn{\boldsymbol{P}}{P} may be interpreted as a vector of
#'    allele frequencies for a specific ancestral population.}
#'    \item{Q_hat}{: The estimated \eqn{\boldsymbol{Q}}{Q} matrix. Each
#'    column of \eqn{\boldsymbol{Q}}{Q} may be interpreted as the admixture
#'    proportions of a specific individual.}
#'    }
#'
#'  @references
#'  Cabreros, I., and J. D. Storey. 2017. “A Nonparametric Estimator of Population Structure Unifying Admixture Models and Principal Components Analysis.” BioRxiv. Cold Spring Harbor Laboratory. doi:10.1101/240812.
#'
#' @export
factor_F <- function(F_hat, d, tol = 0.00001, max_iters = 1000){
  m <- dim(F_hat)[1]
  n <- dim(F_hat)[2]
  P_hat <- matrix(runif(m * d, min = 0, max = 1), nrow = m, ncol = d)

  # simple truncation algorithm for ALS
  Q_hat <- solve(t(P_hat)%*%P_hat, t(P_hat)%*%F_hat)
  P_hat <- t(solve(Q_hat%*%t(Q_hat), Q_hat%*%t(F_hat)))

  iter <- 1
  current_RMSE <- Inf
  Q_hat_old <- Inf*matrix(1, nrow = d, ncol = n)
  while ((iter < max_iters) && (current_RMSE > tol)){
    Q_hat <- solve(t(P_hat)%*%P_hat, t(P_hat)%*%F_hat)
    for(j in 1:n){
      Q_hat[,j] <- projsplx(Q_hat[,j])
    }
    P_hat <- t(solve(Q_hat%*%t(Q_hat), Q_hat%*%t(F_hat)))
    P_hat[P_hat > 1] <- 1; P_hat[P_hat < 0] <- 0;
    iter <- iter + 1
    current_RMSE <- RMSE(Q_hat, Q_hat_old)
    Q_hat_old <- Q_hat
  }

  result <- list(P_hat = P_hat, Q_hat = Q_hat, iter = iter, tol = current_RMSE)
  return(result)
}

#' A provable algorithm for factoring \eqn{\boldsymbol{\hat{F}}}{F_hat}.
#'
#' An algorithm for finding approximate factors \eqn{\boldsymbol{\hat{P}}}{P_hat} and \eqn{\boldsymbol{\hat{Q}}}{Q_hat} of the matrix \eqn{\boldsymbol{\hat{F}}}{F_hat} that obey constraints of the admixture model:
#' \enumerate{
#'  \item \eqn{p_{ij} \in [0,1]  ~ \forall (i,j)}{p_ij in [0, 1] for all (i, j)}
#'  \item \eqn{q_{ij} \ge 0  ~ \forall (i,j)}{q_ij >= 0 for all (i, j)} and \eqn{\sum_i q_{ij} = 1 ~ \forall j}{q_1j + q_2j + ... + q_dj = 1 for all j}
#' }
#' This algorithm is described in Algorithm 1 in (Cabreros and Storey 2017) and
#' has the provable guarantee of asymptotically converging to a stationary
#' point of the objective function.
#'
#' @param F_hat The estimate of the matrix \eqn{\boldsymbol{F}}{F} to be factored
#' @param d The dimension of the latent space. This can be estimated by the
#'        function \code{d_estimate}.
#' @param tol The convergence criterion. If \eqn{RMSE(\boldsymbol{\hat{Q}}_t - \boldsymbol{\hat{Q}}_{t + 1})
#'        < tol}{RMSE(Q_t - Q_(t + 1)) < tol}, then the algorithm halts
#' @param max_iters The maximum number of iterations (repetitions of steps (6)
#'        and (7) in Algorithm 1) to be executed
#' @param P_init Optional initialization of P. If blank, P is initialized
#'        with each entry an independent draw from Unif(0,1)
#' @param Q_init Optional initialization of Q. If blank, Q is initialized as an
#'        all zeros matrix.
#' @param Q_samples Optional number of subsampled of columns of Q.
#' When left unspecified, \code{Q_samples} is set to \eqn{n}{n}. When specified, the
#' \code{cALS} algorithm performs the optimization over a random subset of size
#' \code{Q_samples} during each iteration. This can substantially speed up computation.
#' Only available for \code{cALS} method.
#' @param P_samples Optional number of subsampled of rows of P.
#' When left unspecified, \code{P_samples} is set to \eqn{m}{m}. When specified, the
#' \code{cALS} algorithm performs the optimization over a random subset of size
#' \code{P_samples} during each iteration. This can substantially speed up computation.
#' Only available for \code{cALS} method.
#'
#' @return A list with the following elements:
#'    \describe{
#'    \item{P_hat}{: The estimated \eqn{\boldsymbol{P}}{P} matrix. Each
#'    column of \eqn{\boldsymbol{P}}{P} may be interpreted as a vector of
#'    allele frequencies for a specific ancestral population.}
#'    \item{Q_hat}{: The estimated \eqn{\boldsymbol{Q}}{Q} matrix. Each
#'    column of \eqn{\boldsymbol{Q}}{Q} may be interpreted as the admixture
#'    proportions of a specific individual.}
#'    }
#'
#'  @references
#'  Cabreros, I., and J. D. Storey. 2017. “A Nonparametric Estimator of Population Structure Unifying Admixture Models and Principal Components Analysis.” BioRxiv. Cold Spring Harbor Laboratory. doi:10.1101/240812.
#'
#' @keywords internal
cals <- function(F_hat, d, tol = 0.00001, max_iters = 1000,
                            P_init = NULL, Q_init = NULL, P_samples = NULL, Q_samples = NULL){
  F <- F_hat
  if (is.null(P_init)){
    P_init = matrix(runif(dim(F)[1] * d), nrow = dim(F)[1], ncol = d)
  }

  if (is.null(Q_init)){
    Q_init = matrix(0, nrow = d, ncol = dim(F)[2])
  }

  if (is.null(Q_samples)){
    Q_samples = dim(F)[2]
  }

  if (is.null(P_samples)){
    P_samples = dim(F)[1]
  }

  P = P_init; Q = Q_init

  m <- dim(F)[1]
  n <- dim(F)[2]

  params = list(OutputFlag = 0)

  model_square <- list()
  model_square$lb <- rep(0, d)
  model_square$ub <- rep(1, d)
  model_square$A <- matrix(0, nrow = 1, ncol = d) # fake linear constraint
  model_square$rhs <- 0
  model_square$sense <- '='

  model_triangle <- list()
  model_triangle$A <- matrix(1, nrow = 1, ncol = d)
  model_triangle$sense <- '='
  model_triangle$rhs <- 1
  model_triangle$lb <- rep(0, d)
  model_triangle$ub <- rep(1, d)

  obj_vals <- matrix(0, nrow = 2, ncol = max_iters)

  iter <- 1
  current_RMSE <- Inf
  Q_old <- Inf*matrix(1, nrow = d, ncol = n)
  while ((iter < max_iters) && (current_RMSE > tol)){
    # randomly sample rows of the P matrix and columns of the Q matrix for
    # efficient fitting
    P_samps <- sample(1:m, P_samples, replace = FALSE)
    Q_samps <- sample(1:n, Q_samples, replace = FALSE)
    model_triangle$Q <- t(P) %*% P
    for (j in Q_samps){
      model_triangle$obj <- -2 * t(P) %*% F[, j]
      result <- gurobi::gurobi(model_triangle, params)
      Q[, j] <- result$x
    }
    obj_vals[1,iter] <- result$objval

    model_square$Q <- Q %*% t(Q)
    for (j in P_samps){
      model_square$obj <- -2 * t(F[j, ] %*% t(Q))
      result <- gurobi::gurobi(model_square, params)
      P[j, ] <- result$x
    }
    obj_vals[2, iter] <- result$objval
    iter <- iter + 1
    current_RMSE <- RMSE(Q, Q_old)
    Q_old <- Q
  }

  final_result = list(P_hat = P, Q_hat = Q, obj_vals = obj_vals)
  return(final_result)
}

#' Finds the best permutation of the rows of \eqn{\boldsymbol{\hat{Q}}}{Q_hat}.
#'
#' Returns a  \eqn{(d \times d)}{d x d} permutation matrix \eqn{\boldsymbol{A}}{A} that
#' minimizes the Frobenius distance between \eqn{\boldsymbol{Q}}{Q} and
#' \eqn{\boldsymbol{A}\boldsymbol{\hat{Q}}}{Q_hat}:
#' \deqn{||\boldsymbol{Q} -
#' \boldsymbol{A\hat{Q}}||_F = \sqrt{\sum_{(i,j)} (q_{ij} - \boldsymbol{a}_{i \cdot} \boldsymbol{\hat{q}}_{\cdot j})^2 }}{|| Q - Q_hat ||_F}
#' This is useful when trying to make meaningful comparisons between the output
#' of two methods for estimating global ancestry, or for comparison to a
#' ground truth (when evaluating performance on simulated data). Since for any permutation matrix
#'  \eqn{\boldsymbol{A}}{A} we have the identity
#' \eqn{\boldsymbol{A}^T\boldsymbol{A} = \boldsymbol{I}}{A^T A = I}, the permuted matrices
#' \eqn{\boldsymbol{\tilde{P}} = \boldsymbol{P}\boldsymbol{A}^T}{PA^T} and
#' \eqn{\tilde{\boldsymbol{Q}} = \boldsymbol{A}\boldsymbol{Q}}{AQ} provide a
#' model with likelihood equal to that specified
#' by \eqn{\boldsymbol{P}}{P} and \eqn{\boldsymbol{Q}}{Q} under the admixture model.
#'
#' @param Q the fixed matrix to which \eqn{\boldsymbol{\hat{Q}}}{Q_hat} is compared
#' @param Q_hat the matrix for which we try to find the optimal row permutation
#'
#' @return a \eqn{d \times d}{d x d} permutation matrix
#'
#' @keywords internal
best_perm <- function(Q, Q_hat){

  d <- dim(Q)[1]
  perms <- gtools::permutations(d,d)
  dists <- numeric(dim(perms)[1])

  # compute the distance between Q and each permutation of Q_hat
  for (i in 1:dim(perms)[1]){
    perm <- matrix(0, ncol = d, nrow = d)
    for (j in 1:d){
      perm[j, perms[i,j]] <- 1
    }
    Q_hat_perm <- perm %*% Q_hat
    dists[i] <- norm(Q - Q_hat_perm)
  }

  ind <- which.min(dists)
  perm_mat <- matrix(0, nrow = d, ncol = d)

  for (j in 1:d){
    perm_mat[j, perms[ind, j]] <- 1
  }

  perm_mat
}
