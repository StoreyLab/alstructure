# produces summary log files for ALStructure algorithm execution

#' Function for writing summary of \code{alstructure}
#'
#' Writes a file "summary.txt" to the file specified by the argument
#' "out_folder" that includes details about the ALStructure algorithm execution
#' including method used and time elapsed during exectution.
#'
#' @param method one of \code{cALS} or \code{uALS}
#' @param run_time the amount of time ALStructure tood to run
#' @param out_folder folder to write summary function
#
#' @return None
#'
#' @keywords internal
alstructure_summary <- function(run_time, out_folder){
  summary_str <- paste0("time: ", run_time)

  log_file <- paste0(out_folder, "/summary.txt")
  writeLines(summary_str, con = log_file, sep = "\n")
}

#' Main function for execution of the \code{ALStructure} algorithm
#'
#' Computes global ancestry estimates under the admixture model given a SNP data matrix
#' \eqn{\boldsymbol{X}}{X}. This function is based on the
#' \code{ALStrcture} algorithm from (Cabreros and Storey 2017).
#'
#' @param X The \eqn{m \times n}{m by n} SNP data matrix.
#' @param d_hat Estimate of the latent space dimension \eqn{d}{d}. If left blank,
#'        this is estimated by the function \code{estimate_d()}
#' @param svd_method One of "vanilla" or "trunc." If "vanilla" is chosen, the
#' base \code{svd()} function is used. If "trunc" is used, the truncated svd algorithm
#' from the \code{lfa} package is used.
#' @param tol The convergence criterion. If \eqn{RMSE(\boldsymbol{\hat{Q}}_t - \boldsymbol{\hat{Q}}_{t + 1})
#'        < tol}{RMSE(Q_t - Q_(t + 1)) < tol}, then the algorithm halts
#' @param max_iters The maximum number of iterations (repetitions of steps (6)
#'        and (7) in Algorithm 1) to be executed
#' @param order_method One of "ave_admixture" or "var_explained." If "ave_admixture," the \eqn{d}{d} populations are ordered by decreasing average admixture accross samples (i.e. \eqn{1 / n \sum_{j} q_{ij}}{1 / n (q_i1 + q_i2 + ... + q_in)}). If "var_explained", the \eqn{d}{d} populations are ordered be decreasing variation explained. Specifically, we compute a modified version of the \eqn{\mbox{eigen-}R^2}{eigen-R^2} statistic from (L. S. Chen and Storey 2008). The statistic is modified in the following ways: 1) we treat
#' rows of \eqn{\boldsymbol{Q}}{Q} as the response variables 2) we regress
#' each row of Q on the eigenvectors of \eqn{\boldsymbol{G}}{G} rather than the
#' eigenvectors of the data matrix itself 3) we take the weighted average only
#' over the top \eqn{d}{d} eigenvectors.
#' and columns of \eqn{\boldsymbol{P}}{P} are ordered by amount of variation
#' explained by each row of \eqn{\boldsymbol{Q}}{Q} by the function
#' \code{order_Q}
#' @param P_init Optional initialization of \eqn{\boldsymbol{P}}{P}. Only available for \code{cALS} method.
#' @param Q_init Optional initializtion of \eqn{\boldsymbol{Q}}{Q}. Only available for \code{cALS} method.
#'
#' @return A list with the following elements:
#'    \describe{
#'    \item{P_hat}{: The estimated \eqn{\boldsymbol{P}}{P} matrix. Each
#'    column of \eqn{\boldsymbol{P}}{P} may be interpreted as a vector of
#'    allele frequencies for a specific ancestral population.}
#'    \item{Q_hat}{: The estimated \eqn{\boldsymbol{Q}}{Q} matrix. Each
#'    column of \eqn{\boldsymbol{Q}}{Q} may be interpreted as the admixture
#'    proportions of a specific individual.}
#'    \item{rowspace}{: a list with the following elements:
#'    \describe{
#'        \item{vectors}{: The top \eqn{d}{d} eigenvectors of the matrix
#'        \eqn{\boldsymbol{G}}{G} sorted by decreasing eigenvalue. These vectors
#'    approximate the subspace spanned by the rows of \eqn{\boldsymbol{Q}}{Q}.}
#'        \item{values}{: The top \eqn{d}{d} eigenvalues of the matrix
#'        \eqn{\boldsymbol{G}}{G} sorted by decreasing eigenvalue.}}}
#'    }
#'
#' @references
#' Cabreros, I., and J. D. Storey. 2017. “A Nonparametric Estimator of Population Structure Unifying Admixture Models and Principal Components Analysis.” BioRxiv. Cold Spring Harbor Laboratory. doi:10.1101/240812.
#'
#' Hao, W., M. Song and J. D. Storey. 2015. “lfa: Logistic Factor Analysis for Categorical Data.” R package version 1.8.0, https://github.com/StoreyLab/lfa.
#' @export
alstructure <- function(X, d_hat = NULL,
                        svd_method = "base", tol = 0.00001,
                        max_iters = 1000, order_method = "ave_admixture", P_init,
                        Q_init){

  run_time <- proc.time()

  # if the there is no value for d, approximate d
  if (is.null(d_hat)){
    d_hat <- estimate_d(S, s = 2)
  }

  # impute missing values for X
  X <- alstructure:::impute_mean(X)

  # estimate F_hat
  F_obj <- estimate_F(X, d = d_hat, svd_method = svd_method)
  F_hat <- F_obj$F_hat
  rowspace <- F_obj$rowspace

  # estimate the factors
  factors <- factor_F(F = F_hat, d = d_hat, tol = tol, max_iters = max_iters)

  # order the P and Q matrices
  if(order_method == "ave_admixture"){
    ordering <- order_pops(factors$P_hat, factors$Q_hat, method = "ave_admixture", Q_space = NULL)
    factors$Q_hat <- ordering$Q_ordered
    factors$P_hat <- ordering$P_ordered
  } else if(order_method == "var_explained"){
    ordering <- order_pops(factors$P_hat, factors$Q_hat, method = "var_explained", Q_space = rowspace)
    factors$Q_hat <- ordering$Q_ordered
    factors$P_hat <- ordering$P_ordered
  }

  run_time <- proc.time() - run_time

  vals = list(P_hat = factors$P_hat, Q_hat = factors$Q_hat, rowspace = rowspace, iter = factors$iter, tol = factors$tol)
  return(vals)
}
