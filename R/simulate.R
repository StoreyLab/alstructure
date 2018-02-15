# Functions below are for creating and writing simumlated data sets according
# to the PSD model.

#' function for simulating dirichlet matrices
#'
#' copied from gtools library
#'
#' @keywords internal
rdirichlet <- function (n, alpha)
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

#' Simulates data from the PSD model
#'
#' Creates a data frame that contains the parameters of the admixture model
#' \eqn{(\boldsymbol{F},\boldsymbol{P},\boldsymbol{Q})}{(F, P, Q)} as well as a
#'   single draw \eqn{\boldsymbol{X}}{X} such that \eqn{x_{ij} \sim \mbox{Bernoulli}(f_{ij})}{x_ij ~ Bernoulli(f_ij)}.
#'   The \eqn{\boldsymbol{Q}}{Q} matrix is drawn from the Dirichlet distribution with parameter
#'   \eqn{\alpha}{alpha}, and the \eqn{\boldsymbol{P}}{P} matrix is simulated from the Balding-Nichols model.
#'   The parameter \eqn{\alpha}{alpha} is supplied by the user. The \eqn{m \times 2}{m x 2} matrix
#'   of Balding-Nichols parameters is optional. If \code{BN_params} is not
#'   supplied, the parameters are derived from a random sample of estimated
#'   Balding-Nichols parameters from the Human Genomes Diversity Project (HGDP)
#'   dataset. The Balding-Nichols parameter estimates are provided by (Gopalan et al. 2016), and included in this package in the object \code{hgdpBN}.
#'
#' @param m number of SNPs
#' @param n number of individuals
#' @param d number of groups
#' @param alpha dirichlet parameter; \code{length(alpha) = d}
#' @param BN_params a \eqn{m \times 2}{m x 2} matrix of parameters. The first column
#'        contains \eqn{\textrm{F}_{\textrm{ST}}}{F_ST} for each SNP, while the second
#'         column contains the allele frequency.
#'
#' @return a list with the following elements:
#'    \describe{
#'    \item{P}{: the \eqn{m \times d}{m x d} matrix of loadings}
#'    \item{Q}{: the \eqn{d \times n}{d x n} matrix of latent admixture components}
#'    \item{F}{: the \eqn{m \times n}{m x n} matrix \eqn{PQ}{PQ}}
#'    \item{X}{: a random draw such that \eqn{x_{ij} \sim \mbox{Bernoulli}(f_{ij}, 2)}{x_ij ~ Bernoulli(f_ij, 2)}}.
#'    }
#'
#' @references
#' Gopalan, P., W. Hao, D. M. Blei, and J. D. Storey. 2016. “Scaling probabilistic models of genetic variation to millions of humans.” Nature Genetics 48 (12): 1587–90.
#'
#' @export
simulate_admixture <- function(m, n, d, alpha, BN_params = NA, seed = NA){

  if(length(alpha) != d) message("length(alpha) must equal d")
  if(!is.na(seed)) set.seed(seed)

  if(is.na(BN_params)){
    inds <- sample(1:dim(hgdpBN)[1], m, replace = TRUE)
    BN_params <- hgdpBN[inds, ]
  }

  P = t(apply(BN_params, 1, BN, d))

  Q <- t(alstructure:::rdirichlet(n, alpha))
  F <- P %*% Q
  X <- matrix(rbinom(n * m, 2, c(F)), ncol = n, nrow = m)
  result <- list(P = P, Q = Q, F = F, X = X)
  return(result)
}

#' Simulates \eqn{\boldsymbol{P}}{P} from the BN model
#'
#' Simulates a matrix \eqn{\boldsymbol{P}}{P} from the Balding-Nichols (BN) model.
#' The Balding-Nichols distribution takes two parameters \eqn{(F, p)}{(F, p)}, such that
#' \deqn{\textrm{Balding-Nichols}(F, p) = \textrm{Beta}\left(\frac{1 - F}{F}p, \frac{1 - F}{F}(1 - p) \right)}{Balding-Nichols(F, p) = Beta(p(1 - F) / F, (1 - p)(1 - F) / F)}
#' The function requires an \eqn{m \times 2}{m x 2} matrix of parameters, with each row
#' a pair of parameters \eqn{(F_i, p_i)}{(F_i, p_i)} for each SNP.
#'
#' @param BN_params a \eqn{m \times 2}{m x 2} matrix of parameters. The first column
#'        contains \eqn{\textrm{F}_{\textrm{ST}}}{F_ST} for each SNP, while the second
#'         column contains the allele frequency (\eqn{F}{F} and \eqn{p}{p}, repectively).
#' @param d the dimension of the latent space.
#'
#' @return an \eqn{m \times d}{m x d} matrix \eqn{\boldsymbol{P}}{P}.
#'
#' @keywords internal
BN = function(BN_params, d){
  fst=BN_params[1]
  freq=BN_params[2]
  if(fst < 1e-6)
    p_ind = rep(freq,d)
  else{
    param1 = freq*(1-fst)/fst
    param2 = (1-freq)*(1-fst)/fst
    p_ind = rbeta(d, param1, param2)
    #regularity...could be inf loop
    #while((sum(p_ind<0.05)+sum(p_ind>0.95)) > 0){
    #  p_ind[(p_ind<0.05) | (p_ind>0.95)] = rbeta( sum((p_ind<0.05) |
    #        (p_ind>0.95)) , param1, param2)
    #}
  }

  p_ind
}
