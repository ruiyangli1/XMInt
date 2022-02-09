## data generation function
#'
#' @title Data generation function
#' @description This function generates the simulation data (exposure X, outcome Y, and potential mediators M).
#' @description For each subject i = 1, ..., N:
#' - X_i ~ i.i.d N(0, 1)
#' - M_{i,v} = a_v X_i + e{1_{i,v}}, where e{1_{i,v}} ~ i.i.d N(0, 1), v = 1, ..., V
#' - Y_i = X_i + sum_v b{1_{i,v}} M_{i,v} + sum_v b{2_{i,v}} X_i x M_{i,v} + e{2_{i,v}}, where e{2_{i,v}} ~ i.i.d. N(0, 1)
#' @description The first three M variables (M1,M2,M3) are set to be the true mediators (i.e., having non-zero a and b_1 coefficients), X x M1 is set to be the true exposure-by-mediator interaction term (i.e., having non-zero b_2 coefficients), and all other coefficients are set to be 0.
#'
#' @param N sample size
#' @param V number of potential mediators
#' @param es effect size, representing the value of a, b1, b2 of the truth
#' @param seed seed. Default is 1234.
#'
#' @return The resulting dataset has 3 lists: X, Y and M
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

