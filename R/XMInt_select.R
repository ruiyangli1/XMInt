## algorithm function to select mediators/interaction

#' @title Model selection results
#' @description This function gives the model selection results.
#'
#' @param X exposure
#' @param Y outcome
#' @param M mediators
#' @param K lambda sequence length (default = 20)
#' @param zeta min/max lambda ratio (default = 0.05)
#' @param max.factor enlarging factor (default = 1.5) to start with base model
#' @param hbic_plot plot HBIC curve (default = \code{FALSE})
#' @param coef_print print coefficients (default = \code{FALSE})
#'
#' @return \emph{selected_mediator}: selected mediator
#' @return \emph{selected_interaction}: selected interaction
#' @return hbic_plot: HBIC plot for final model selection
#' @return hbic: computed HBIC
#' @return lambda: lambda sequence used
#' @return coefficient: estimated coefficients
#' @export
#'
#' @examples
#' data = dat_gen(N = 400, V = 50, es = 1, seed = 1234)
#' X = data$X; Y = data$Y; M = data$M
#' XMInt_select(X,Y,M)
XMInt_select <- function(X,Y,M,
                         K = 20,
                         zeta = 0.05,
                         max.factor = 1.5,
                         hbic_plot = FALSE, coef_print = FALSE){

  hbic_calc <- function(fit, fit.n){

    # calculate last term (sum)
    sum1 = 0
    sum2 = 0
    for (j in 1:N) {
      s1 = as.numeric( t(M[j,] - fit$hata * X[j]) %*% data.matrix(data.frame(fit$Omega)) %*% (M[j,] - fit$hata * X[j]))
      sum1 = sum1 + s1
      s2 = as.numeric(t(M[j,] - fit.n$hata * X[j]) %*% data.matrix(data.frame(fit.n$Omega)) %*% (M[j,] - fit.n$hata * X[j]) )
      sum2 = sum2 + s2
    }

    # calculate difference in 2*log-likelihood (l(fit.n) - l(fit))
    lik.diff = -N*log(fit.n$sigmasq)+N*log(fit$sigmasq) +
      N*log(det(data.matrix(data.frame(fit.n$Omega))))-
      N*log(det(data.matrix(data.frame(fit$Omega)))) -
      1/fit.n$sigmasq * rowsum(((Y - X*fit.n$c - M %*% fit.n$hatb1 - (X*M) %*% fit.n$hatb2))^2, rep(1,N)) +
      1/fit$sigmasq * rowsum(((Y - X*fit$c - M %*% fit$hatb1 - (X*M) %*% fit$hatb2))^2, rep(1,N)) -
      sum2 + sum1

    # calculate hbic
    a = fit.n$Omega
    fit.n.covnpar = sum(a[lower.tri(a)] !=0)

    b = fit$Omega
    fit.covnpar = sum(b[lower.tri(b)] !=0)

    hbic = lik.diff - ((fit.n$nump+fit.n.covnpar) - (fit$nump+fit.covnpar))*log(N/(2*pi))
    names(hbic) = NULL

    #message("Null model's npar: --- coef: ", fit.n$nump, " // Omega: ", fit.n.covnpar, " // total: ",fit.n$nump+fit.n.covnpar)
    #message("Fitted model's npar: --- coef: ", fit$nump, " // Omega: ", fit.covnpar, " // total: ",fit$nump+fit.covnpar)

    return(hbic)
  }

  # data preparation

  X = as.numeric(scale(X))
  Y = as.numeric(scale(Y))
  M = as.matrix(scale(M))

  M_name = colnames(M)

  N = as.numeric(nrow(M))
  V = as.numeric(ncol(M))

  colnames(M) = paste0("M",1:V)
  I = X*M
  colnames(I) = paste0("X*", colnames(M))



  # Initialization

  I_update = I
  penalty = c(0,rep(1,ncol(M)),rep(1,ncol(I_update)),rep(1,ncol(M))) # penalty (c,b1,b2,a)

  X_full = cbind(X,M,I)

  lambda_max = 1/N*max(abs(t(X_full) %*% Y))

  fit = try(model_estimate(X,M,Y,I_update, lambda1 = lambda_max,lambda2 = exp(-1),alpha = 1,penalty.factor = penalty,Omega.out = TRUE))
  if (fit$nump > 1) {
    cat("Max lambda needs to be enlarged. \n")
  }
  while (fit$nump > 1) {
    lambda_max = max.factor*lambda_max
    fit = try(model_estimate(X,M,Y,I_update, lambda1 = lambda_max,lambda2 = exp(-1),alpha = 1,penalty.factor = penalty,Omega.out = TRUE))
  }
  cat("Maximum lambda has been identified. \n")

  lambda_min = zeta * lambda_max
  lambda_seq = exp(seq(log(lambda_max), log(lambda_min), length = K))


  # Initialization (cont'd)

  hbic = NULL
  coef = NULL
  penalty_hist = NULL


  # Path-building
  # ~~~ iteration starts: on lambda_seq[i]

  ## null model
  fit.n = try(model_estimate(X, M, Y, I_update, lambda1 = exp(1), lambda2 = exp(-1), alpha = 1, Omega.out = TRUE))

  for (i in 1:K) {

    ## interaction model
    fit = try(model_estimate(X,M,Y,I_update, lambda1 = lambda_seq[i],lambda2 = exp(-1),alpha = 1,penalty.factor = penalty,Omega.out = TRUE))

    a.name = names(fit$hata[fit$hata != 0])
    b1.name = names(fit$hatb1[fit$hatb1 != 0])
    b2.name = names(fit$hatb2[fit$hatb2 != 0])

    # update penalty (c,b1,b2,a)
    ## not penalize a, b1 for M identified as mediators or from interaction
    ## not penalize b2 for M identified from interaction

    ## identify the location of mediators
    M.id = intersect(a.name,b1.name)
    M.id.loc = which(colnames(M) %in% M.id)

    ## identify the location of interactions
    I.id.loc <- gsub("[^0-9]", "", b2.name) # remove everything except 0-9
    I.id.loc <- as.numeric(I.id.loc)
    ### corresponding location of M
    M.id.loc.from.I = which(colnames(M) %in% b2.name)

    #message("M: ", paste0(M.id," "))
    #message("I: ", paste0(b2.name," "))

    ## penalty
    penalty = c(1,rep(1,ncol(M)),rep(1,max(0,ncol(I_update))),rep(1,ncol(M)))
    penalty[c(1,
              1 + union(M.id.loc, M.id.loc.from.I),
              1 + ncol(M) + I.id.loc,
              1 + ncol(M) + max(0,ncol(I_update)) + union(M.id.loc, M.id.loc.from.I))] = 0

    ## hbic
    hbic[i] = hbic_calc(fit,fit.n)

    ## coefficients
    coef[[paste0("lambda", i)]] = data.frame(
      a = fit$hata,
      b1 = fit$hatb1,
      b2 = fit$hatb2,
      c = fit$c)

    ## store penalty for re-estimation
    penalty_hist[[i]] = penalty

    # ~~ iteration ends, goes back
    cat(paste("Lambda", i, "was finished.\n"))

  }


  # selected optimal lambda (that gives min hbic)
  hbic_min = which.min(hbic)

  # visualize hbic
  if (hbic_plot == FALSE) {
    hbic_plt = "HBIC plot is not printed."
  } else {
    hbic_plt = ggplot(data = data.frame(x = 1:length(hbic),hbic = hbic),aes(x = x, y = hbic)) +
      geom_point() +
      labs(x = "lambda index", y = "HBIC") +
      geom_point(aes(x = hbic_min, y = min(hbic)),
                 pch = 5, size = 4) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }

  # SELECTED coeff results
  result_optimal = coef[[paste0("lambda",hbic_min)]]
  rownames(result_optimal) = M_name

  # re-estimation
  fit = try(model_estimate(X,M,Y,I_update, lambda1 = lambda_seq[hbic_min],lambda2 = exp(-1),alpha = 1,penalty.factor = penalty_hist[[hbic_min]],Omega.out = TRUE))
  result_optimal_re = data.frame(
    a = fit$hata,
    b1 = fit$hatb1,
    b2 = fit$hatb2,
    c = fit$c)
  rownames(result_optimal_re) = M_name
  coeff_m_re = subset(result_optimal_re, a != 0 | b1 != 0 | b2 != 0)

  coeff_m = subset(result_optimal, a != 0 | b1 != 0 | b2 != 0)
  int_selected = rownames(subset(coeff_m, b2 != 0)) # selected interaction
  med_selected = rownames(subset(coeff_m, (a != 0 & b1 != 0) | b2!=0)) # selected M

  if (coef_print == FALSE) {
    coefficient = "Coefficients are not printed."
  } else {
    coefficient = subset(coeff_m_re, a!=0&b1!=0)
  }

  return(list(selected_mediator = med_selected, selected_interaction = int_selected,
              hbic_plot = hbic_plt, hbic = hbic,
              lambda = lambda_seq,
              coefficient = coefficient))
}
