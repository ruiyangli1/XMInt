## data generation function
#' Title Data generation function
#'
#' @title Data generation function
#' @description This function generates the simulated data.
#'
#' @param N sample size
#' @param V number of potential mediators
#' @param es effect size
#' @param seed seed. Default is 1234.
#'
#' @return dataset with 3 lists, containing the exposure X, the outcome Y and the potential V mediators M
#' @export
#'
#' @examples
#' dat_gen(100, 50, 1)
dat_gen <- function(N,V,es,seed = 1234){
  a = c(rep(es,3),rep(0,V - 3))
  b1 = b2 = rep(0,V)
  b1[1:3] <- es
  b2[1] <- es

  set.seed(seed)
  X = rnorm(N, mean = 0, sd = 1)
  M = X %*% t(a) + matrix(rnorm(N*V, mean = 0, sd = 1), nrow = N, ncol = V)
  Y = X + M %*% b1 + (X*M) %*% b2 + rnorm(N, mean = 0, sd = 1)

  X = as.numeric(scale(X))
  Y = as.numeric(scale(Y))
  M = as.matrix(scale(M))
  colnames(M) = paste0("M",1:V)

  #data.frame(cbind(X,Y,M))
  list(X = X, Y = Y, M = M)

}

