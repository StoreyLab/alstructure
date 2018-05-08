# Functions here are used for the conversion of SNP data matrices to .bed
# formats amenable to other algorithms: fastSTRUCTURE, Admixture,
# and terastructure. We also include functions to simultaneously test many
# datasets using ALStructure and the competing algorithms.

#' Creates \code{.bed} files from matrix
#'
#' Creates \code{.bed}, \code{.fam}, and \code{.bim} files from an SNP matrix for the purposes of comparing \code{ALStructure} to other methods. The \code{.bim} and \code{.fam} files are artificial, however we include them because other methods may require them as input.
#'
#' @param X The \eqn{m \times n}{m x n} SNP matrix of 0, 1, and 2.
#' @param out_file The output file. No extension should be added. Three files,
#'        \code{out_file.bed}, \code{out_file.fam}, and \code{out_file.bim} will be created.
#' @param B The number of blocks to to write data in. Default is set to
#'        \code{ceil(m / 10)}.
#'
#' @return Three files written to \code{out_file}:
#'          \describe{
#'          \item{\code{out_file.bed}}{SNP matrix in \code{.bed} format}
#'          \item{\code{out_file.fam}}{artificial file}
#'          \item{\code{out_file.bim}}{artificial file}}
#' @keywords internal
make_bed <- function(X, out_file, B = min(ceiling(dim(X)[2]/10))){

  m <- nrow(X)
  n <- ncol(X)

  # number of zeros needed to append in order for n to be divisible by 4
  n_diff <- (4 - (n %% 4)) %% 4
  buff_mat <- matrix(0, nrow = m, ncol = n_diff)
  X <- cbind(X, buff_mat)

  f <- file(sprintf("%s.bed", out_file), open = "wb")
  close(f)
  f <- file(sprintf("%s.bed", out_file), open = "ab")

  # required for PLINK formatting
  writeBin(as.raw(c(108, 27, 1)), f)

  # breaks SNP matrix into chunks of size B to write binary file
  splits <- split(1:m, ceiling(seq_along(1:m)/B))
  for(i in 1:length(names(splits))){

    # look at ith B-block of rows
    x <- X[splits[[i]], ]
    intsnp <- as.vector(t(x))

    # order matters
    intsnp[intsnp == 0] <- 3
    intsnp[intsnp == 2] <- 0
    intsnp[intsnp == 1] <- 2
    intsnp[is.na(intsnp)] <- 1

    # no missing data
    intsnp <- matrix(intsnp, 4)
    intsnp[2,] <- intsnp[2,] * 4
    intsnp[3,] <- intsnp[3,] * 16
    intsnp[4,] <- intsnp[4,] * 64
    intsnp <- colSums(intsnp)
    writeBin(as.raw(intsnp), f)
  }

  close(f)

  # write phony .fam and .bim files required for fastSTRUCTURE
  trait <- 1

  out <- cbind(1:n, 1:n, 4, 5, 1, trait)
  out <- format(out, scientific = FALSE)
  write.table(out, file = sprintf("%s.fam", out_file), quote = FALSE, sep = " ",
              row.names = FALSE, col.names = FALSE)

  out <- cbind(1, 1:m, 1, 1:m, 1, 1)
  out <- format(out, scientific=FALSE)
  write.table(out, file = sprintf("%s.bim", out_file), quote = FALSE, sep = " ",
              row.names = FALSE, col.names = FALSE)
}

#' Orders the \eqn{d}{d} populations
#'
#' Orders the \eqn{d}{d} populations according to one of two methods: "ave_admixture" or "var_explained." Function returns matrices \eqn{\boldsymbol{P}}{P} and \eqn{\boldsymbol{Q}}{Q} with permuted columns and rows according to the determined ordering of populations.
#'
#' @param P The \eqn{m \times d}{m x d} loadings matrix with columns to be
#'  ordered
#' @param Q The \eqn{d \times n}{d x n} admixture matrix with rows to be ordered
#' @param method One of "ave_admixture" or "var_explained." If "ave_admixture," the \eqn{d}{d} populations are ordered by decreasing average admixture accross samples (i.e. \eqn{1 / n \sum_{j} q_{ij}}{1 / n (q_i1 + q_i2 + ... + q_in)}). If "var_explained", the \eqn{d}{d} populations are ordered be decreasing variation explained. Specifically, we compute a modified version of the \eqn{\mbox{eigen-}R^2}{eigen-R^2} statistic from (L. S. Chen and Storey 2008). The statistic is modified in the following ways: 1) we treat
#' rows of \eqn{\boldsymbol{Q}}{Q} as the response variables 2) we regress
#' each row of Q on the eigenvectors of \eqn{\boldsymbol{G}}{G} rather than the
#' eigenvectors of the data matrix itself 3) we take the weighted average only
#' over the top \eqn{d}{d} eigenvectors.
#' @param Q_space Only required for "var_explained" methodThe list containing
#' the top \eqn{d}{d} eigenvectors and their
#' corresponding eigenvalues of the \eqn{\boldsymbol{G}}. These may be obtained
#' through the function \code{lse}.
#'
#' @return A list with the following elements:
#'    \describe{
#'    \item{P_ordered}{: The permuted \eqn{\boldsymbol{P}}{P} }
#'    \item{Q_ordered}{: The permuted \eqn{\boldsymbol{Q}}{Q} }
#'    \item{perm_mat}{: The permutation matrix \eqn{\boldsymbol{A}}{A} such that
#'    \eqn{\boldsymbol{PA}^T = \boldsymbol{P}_{\mbox{ord}}}{PA^T = P_ord} and \eqn{\boldsymbol{AQ} = \boldsymbol{Q}_{\mbox{ord}}}{AQ = Q_ord}.}
#'    }
#'
#' @references
#' Chen, L. S., and J. D. Storey. 2008. “Eigen-R2 for dissecting variation in high-dimensional studies.” Bioinformatics 24 (19): 2260–2.
#' @export
order_pops <- function(P, Q, method = "ave_admixture", Q_space = NULL){

  d <- dim(Q)[1]
  if(method == "ave_admixture"){
    n <- dim(Q)[2]
    ave_adx <- rep(0, d)
    ave_adx <- apply(Q, 1, sum)

    # sort in decreasing order
    row_order <- sort(ave_adx, index.return = TRUE, decreasing = TRUE)

    # create permutation matrix
    A <- matrix(0, nrow = d, ncol = d)
    for(i in 1:d) A[i, row_order$ix[i]] <- 1

    # order Q
    P_ordered <- P %*% t(A)
    Q_ordered <- A %*% Q

  } else if(method == "var_explained"){
    if(is.null(Q_space)){stop("Q_space must be supplied for var_explained method")}
    pi <- Q_space$values / sum(Q_space$values)
    eigen_r2 <- rep(0, d)
    R2 <- matrix(0, nrow = d, ncol = d)

    for(i in 1:d){
      R2 <- rep(0, d)
      for(j in 1:d){
        df <- data.frame(Qi = Q[i, ], vj = Q_space$vectors[, j])
        fit <- lm(Qi ~ vj, data = df)
        R2[j] <- summary(fit)$r.squared
      }
      eigen_r2[i] <- sum(pi * R2)
    }

    # sort in decreasing order
    row_order <- sort(eigen_r2, index.return = TRUE, decreasing = TRUE)

    # create permutation matrix
    A <- matrix(0, nrow = d, ncol = d)
    for(i in 1:d) A[i, row_order$ix[i]] <- 1

    # order Q
    P_ordered <- P %*% t(A)
    Q_ordered <- A %*% Q
  }

  vals <- list(P_ordered = P_ordered, Q_ordered = Q_ordered, perm_mat = A)
  return(vals)
}

