# Functions below are used for obtaining the estimate hat(F) from the raw data
# X ~ F.

#' Estimates the \eqn{n \times n}{n x n} matrix \eqn{\boldsymbol{D}}{D} for Latent Subspace
#' Estimation.
#'
#' Estimates the \eqn{n \times n}{n x n} matrix \eqn{\boldsymbol{D}}{D} using the estimator
#' from Chen and Storey (2015) (Table 1). \eqn{\boldsymbol{D}}{D} is used to adjust
#' for heteroskedasticity when estimating the latent subspace. This function is
#' used by \code{F_estimate()}.
#'
#' @param X the \eqn{m \times n}{m x n} SNP data matrix
#'
#' @return the \eqn{n \times n}{n x n} matrix \eqn{\boldsymbol{\hat{D}}}{D_hat}
#'
#' @references
#' CHEN, X., and J. D. STOREY, 2015 Consistent Estimation of Low-Dimensional
#' Latent Structure in High-Dimensional Data. arXiv:1510.03497v1.
#'
#' @keywords internal
D_binomial <- function(X){
  s <- 2
  m <- dim(X)[1]; n <- dim(X)[2]
  delta_hat <- rep(0,n)

  v <- function(x,s) (s*x - x^2)/(s-1)
  for (i in 1:n){
    delta_hat[i] <- 1/m * sum(v(X[,i],s))
  }

  D <- diag(delta_hat)
}

#' Estimates the individual-specific allele frequency matrix \eqn{\boldsymbol{F}}{F}
#'
#' Estimates the \eqn{m \times n}{m x n} individual-specific allele frequency matrix
#' \eqn{\boldsymbol{F}}{F} using the method of Latent Subspace Estimation
#'  (X. Chen and Storey 2015) as described in (Cabreros and Storey 2017)
#' in Section 2.3.
#'
#' @param X The \eqn{m \times n}{m x n} SNP data matrix
#' @param d The rank of \eqn{\boldsymbol{F}}{F}. This can be estimated using the function
#'        \code{d_estimate()}.
#' @param svd_method One of "base" or "truncated_svd." If "base" is chosen, the
#' base \code{svd()} function is used. If "truncated_svd" is used, the truncated svd algorithm
#' from the \code{lfa} package is used.
#'
#' @return A list with the following elements:
#'    \describe{
#'    \item{F_hat}{: The \eqn{m \times n}{m x n} matrix \eqn{\boldsymbol{\hat{F}}}{F_hat}.}
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
#' Chen, X., and J. D. Storey. 2015. “Consistent Estimation of Low-Dimensional Latent Structure in High-Dimensional Data.” ArXiv E-Prints, October.
#'
#' @export
estimate_F <- function(X, d, svd_method = "base"){
  m <- dim(X)[1]
  # compute D (see Chen and Storey)
  D <- D_binomial(X)

  rowspace <- lse(X = X, d = d, svd_method = svd_method)

  F_hat <- 1/2 * X %*% rowspace$vectors[, 1:d] %*% t(rowspace$vectors[, 1:d])
  # simple truncation of out-of-bound numbers
  F_hat[F_hat < 0] <- 0; F_hat[F_hat > 1] <- 1

  vals <- list(F_hat = F_hat, rowspace = rowspace)
}

#' Estimate the latent space dimension
#'
#' Estimates the dimension of the rowspace of \eqn{\boldsymbol{Q}}{Q}
#'  (equivalently, the rank
#' of the matrix \eqn{\boldsymbol{F}}{F}). This estimate \eqn{\hat{d}}{d_hat} is based on the estimator from (Leek 2011), page 6.
#'
#' @param X the \eqn{m \times n}{m x n} SNP data matrix
#'
#' @return an estimate \eqn{\hat{d}}{d_hat} of the dimension of the latent space
#'         dimension.
#'
#' @references
#' Leek, J. T. 2011. “Asymptotic conditional singular value decomposition for high-dimensional genomic data.” Biometrics 67 (4): 344–52.
#'
#' @export
estimate_d <- function(X){
  m <- dim(X)[1]; n <- dim(X)[2]
  cm <- n*m^(-1/3)
  D <- D_binomial(X)
  Wm <- 1 / m * t(X) %*% X - D
  e <- eigen (Wm)
  d_hat <- sum(e$values > cm)
  d_hat
}

#' Estimates the latent subspace
#'
#' Estimates the rowspace of Q using the method of latent subspace estimation.
#' The function returns the top d eigenvalues and vectors of the matrix
#' \deqn{\boldsymbol{G} = \frac{1}{m} \boldsymbol{X}^T \boldsymbol{X} - \boldsymbol{D}}{G = 1/m X^T X - D}
#' where the matrix D is a diagonal matrix with each diagonal entry
#' \eqn{d_{ii}}{d_ii} an estimate of the average of the variances of the
#' random variables in the \eqn{i}{i} column of \eqn{\boldsymbol{X}}{X}. As is
#' proven in (X. Chen and Storey 2015), the span of the top \eqn{d}{d} eigenvectors
#' of \eqn{\boldsymbol{G}}{G} span the same space as the rows of
#' \eqn{\boldsymbol{Q}}{Q}. The eigenvectors are returned in order of decreasing
#' eigenvalue.
#'
#' @param X The \eqn{m \times n}{m x n} SNP data matrix
#' @param d The rank of \eqn{\boldsymbol{F}}{F}. This can be estimated using the
#'  function \code{d_estimate()}. When \eqn{d = n}{d = n}, all eigenvectors
#'  of \eqn{\boldsymbol{G}}{G} are returned.
#' @param svd_method One of "base" or "truncated_svd." If "base" is chosen, the
#' base \code{svd()} function is used. If "truncated_svd" is used, the truncated svd algorithm
#' from the \code{lfa} package is used.
#'
#' @return A list with the following elements:
#'    \describe{
#'    \item{vectors}{: The top \eqn{d}{d} eigenvectors of the matrix
#'    \eqn{\boldsymbol{G}}{G} sorted by decreasing eigenvalue. These vectors
#'    approximate the subspace spanned by the rows of \eqn{\boldsymbol{Q}}{Q}}.
#'    \item{values}{: The top \eqn{d}{d} eigenvalues of the matrix
#'    \eqn{\boldsymbol{G}}{G} sorted by decreasing eigenvalue.}
#'    }
#'
#' @references
#' Chen, X., and J. D. Storey. 2015. “Consistent Estimation of Low-Dimensional Latent Structure in High-Dimensional Data.” ArXiv E-Prints, October.
#'
#' @export
lse <- function(X, d, svd_method = "base"){
  m <- dim(X)[1]
  # compute D (see Chen and Storey)
  D <- D_binomial(X)

  if(svd_method == "base"){
    # find the rowspace
    rowspace <- eigen(1 / m * t(X) %*% X - D)
    vectors <- rowspace$vectors[, 1:d]
    values <- rowspace$values[1:d]
  } else if(svd_method == "truncated_svd"){
    # find the rowspace using truc.svd
    rowspace <- lfa:::trunc.svd(1 / m * t(X) %*% X - D, d = d)
    vectors <- rowspace$v[, 1:d]
    values <- rowspace$d[1:d]
  }

  vals <- list(vectors = vectors, values = values)
  return(vals)
}
