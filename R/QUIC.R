## QUIC function (from https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz)

#' @title QUIC function replicated from https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz
#'
#' @description QUadratic Inverse Covariance estimation
#' @description Estimates a sparse inverse covariance matrix using a combination of Newton's method and coordinate descent.
#' @description For details, please refer to https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz
#'
#' @param S Covariance matrix. A p by p symmetric matrix.
#' @param rho Regularization parameter. It can be a p by p matrix, a vector or scalar.
#' @param path If specified, then rho is scaled with the elements of path and the corresponding inverse covariance matrix estimation is carried out for each value.
#' @param tol Specifes the convergence tolerance.
#' @param msg Controls how verbose messages should be printed during execution. Valid value range: 0–4.
#' @param maxIter Specifies the maximum number of Newton iterations.
#' @param X.init The initial estimate for the regularized inverse covariance matrix.
#' @param W.init The inverse of initial estimate for the regularized inverse covariance matrix.
#'
#' @return
#' @export
#'
#' @author Matyas A. Sustik (package maintainer), Cho-Jui Hsieh, Inderjit S. Dhillon, Pradeep Ravikumar
#' @references Sparse Inverse Covariance Matrix Estimation Using Quadratic Approximation. Cho-Jui Hsieh, Matyas A. Sustik, Inderjit S. Dhillon, Pradeep Ravikumar, Advances in Neural Information Processing Systems, vol. 24, 2011, p. 2330–2338.
#' @references http://www.cs.utexas.edu/users/sustik/papers/invcov.pdf
#'
QUIC <- function(S, rho, path = NULL, tol = 1.0e-4, msg = 1, maxIter = 1000,
                 X.init = NULL, W.init = NULL) {
  ### $Id: QUIC.R,v 1.7 2012-05-01 02:12:19 sustik Exp $

  n <- nrow(S)
  if (is.null(path))
    npath <- 1
  else
    npath <- length(path)
  if (!is.matrix(rho) && length(rho) != 1 && length(rho) != n) {
    stop("Wrong number of elements in rho")
  }
  if (is.vector(rho)){
    rho <- matrix(sqrt(rho))%*%sqrt(rho)
  }
  if (length(rho) == 1){
    rho <- matrix(rho, ncol = n, nrow = n)
  }
  if (is.null(path)) {
    if (is.null(X.init)) {
      X <- diag(n)
      W <- diag(n)
    } else {
      X <- X.init
      W <- W.init
    }
  } else {
    if (is.null(X.init)) {
      X <- array(diag(n), c(n, n, npath))
      W <- array(diag(n), c(n, n, npath))
    } else {
      X <- array(0, c(n, n, npath))
      W <- array(0, c(n, n, npath))
      X[, , 1] <- X.init
      W[, , 1] <- W.init
    }
  }
  opt <- matrix(0, ncol = npath, nrow = 1)
  cputime <- matrix(0, ncol = npath, nrow = 1)
  iter <- matrix(0, ncol = npath, nrow = 1)
  dGap <- matrix(0, ncol = npath, nrow = 1)
  if (is.null(path))
    job <- "d"
  else
    job <- "p"

  storage.mode(job) <- "character"
  storage.mode(S) <- "double"
  storage.mode(rho) <- "double"
  storage.mode(npath) <- "integer"
  storage.mode(path) <- "double"
  storage.mode(tol) <- "double"
  storage.mode(msg) <- "integer"
  storage.mode(maxIter) <- "integer"
  storage.mode(X) <- "double"
  storage.mode(W) <- "double"
  storage.mode(opt) <- "double"
  storage.mode(cputime) <- "double"
  storage.mode(iter) <- "integer"
  storage.mode(dGap) <- "double"
  tmp<-.C("QUICR", job, n, S, rho, npath, path, tol, msg, maxIter,
          X = X, W = W, opt = opt, cputime = cputime, iter = iter,
          dGap = dGap)
  return (list(X = tmp$X, W = tmp$W, opt = tmp$opt, cputime = tmp$cputime,
               iter = tmp$iter, regloglik = -(n/2)*tmp$opt,
               dGap = tmp$dGap))
}
