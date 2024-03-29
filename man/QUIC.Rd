% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QUIC.R
\name{QUIC}
\alias{QUIC}
\title{QUIC function replicated from https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz}
\usage{
QUIC(
  S,
  rho,
  path = NULL,
  tol = 1e-04,
  msg = 1,
  maxIter = 1000,
  X.init = NULL,
  W.init = NULL
)
}
\arguments{
\item{S}{Covariance matrix. A p by p symmetric matrix.}

\item{rho}{Regularization parameter. It can be a p by p matrix, a vector or scalar.}

\item{path}{If specified, then rho is scaled with the elements of path and the corresponding inverse covariance matrix estimation is carried out for each value.}

\item{tol}{Specifes the convergence tolerance.}

\item{msg}{Controls how verbose messages should be printed during execution. Valid value range: 0–4.}

\item{maxIter}{Specifies the maximum number of Newton iterations.}

\item{X.init}{The initial estimate for the regularized inverse covariance matrix.}

\item{W.init}{The inverse of initial estimate for the regularized inverse covariance matrix.}
}
\value{

}
\description{
QUadratic Inverse Covariance estimation

Estimates a sparse inverse covariance matrix using a combination of Newton's method and coordinate descent.

For details, please refer to https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz
}
\references{
Sparse Inverse Covariance Matrix Estimation Using Quadratic Approximation. Cho-Jui Hsieh, Matyas A. Sustik, Inderjit S. Dhillon, Pradeep Ravikumar, Advances in Neural Information Processing Systems, vol. 24, 2011, p. 2330–2338.

http://www.cs.utexas.edu/users/sustik/papers/invcov.pdf
}
\author{
Matyas A. Sustik (package maintainer), Cho-Jui Hsieh, Inderjit S. Dhillon, Pradeep Ravikumar
}
