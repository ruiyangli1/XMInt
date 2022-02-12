## estimation function

#' @title Model estimation
#' @description This function gives the model estimates.
#'
#' @param X one-dimensional exposure
#' @param M multivariate mediators
#' @param Y one-dimensional outcome
#' @param I_update interaction term
#' @param X.scale whether scale X (default = \code{TRUE})
#' @param tol convergence criterion (default = -10^(-10))
#' @param max.iter maximum iteration (default = 100)
#' @param lambda1 tuning parameter for regression coefficient L1 penalization
#' @param lambda2 tuning parameter for covariance matrix
#' @param alpha alpha in glmnet() (default = 1: lasso penalty)
#' @param penalty.factor penalty factor vector, in the order of (c,b1,b2,a)
#' @param verbose print progress (default = \code{TRUE})
#' @param Omega.out output Omega estimates (default = \code{FALSE})
#' @param non.zeros.stop when to stop the regularization path (default = ncol(M))
#'
#' @return c: direct effect estimate
#' @return hatb1: path b1 (M->Y given X) estimates
#' @return hatb2: path b2 (X*M->Y) estimates
#' @return hata: path a (X->M) estimates
#' @return nump: number of selected paths + 1 direct effect
#' @return Omega: estimated covariance matrix of the mediators
#' @return sigmasq: estimated variance of the outcome
#' @export
#'
#' @examples
#' data = dat_gen(200,100,es = 1)
#' X = data$X; Y = data$Y; M = data$M; I = X*M
#' model_estimate(X, M, Y, I, lambda1 = 0.2, lambda2 = 0.1, alpha = 1, Omega.out = F)

model_estimate =
  function(X, M, Y,
           I_update,
           X.scale = TRUE,
           tol = 10^(-10),
           max.iter = 10,
           lambda1 = exp(seq(1,-3,length = 10)), ###will use lambda_seq[i] in ramp step
           lambda2 = exp(-1), ###for omega estimation
           alpha = 1,
           penalty.factor = c(1,rep(1,ncol(M)*2),rep(1,max(0,ncol(I_update)))), ###here set to penalize (c,b1,b2,a) to be consistent w/ beta coeff order, will use specific penalty in ramp step
           verbose = FALSE,
           Omega.out = FALSE,
           non.zeros.stop = ncol(M)){

    ## Center all values, and also make their scales to be 1. In this context, all coefficients will be described in terms of correlation or partial correlations.
    N = nrow(M)
    V = ncol(M)

    Y.sd = as.vector(sqrt(var(Y)))
    X.sd = as.vector(sqrt(var(X)))
    M.sd = sqrt(apply(M,2,var))

    Y = scale(Y,center = FALSE,scale = TRUE)

    if (X.scale == TRUE) {
      X = matrix(scale(X,center = TRUE,scale = TRUE),N,1)
    }

    M = scale(M, center = FALSE,scale = TRUE)
    XM = I_update

    if (ncol(X) > 1) {stop("X has more than 1 column. Stop.")}


    ## Initialization###
    ## OLS Estimation ###
    U = cbind(X,M,XM)
    tXX = t(X) %*% X
    tUY = t(U) %*% Y
    tMX = t(M) %*% X

    ## from https://github.com/seonjoo/smm
    sqrtmat.comp <- function(mat, thresh = 10^(-20), K = NULL){
      if (is.null(K)) {K = ncol(mat)}
      if (ncol(mat) > 200) {
        eigenmat = rsvd(mat, k = K)
      } else {eigenmat = svd(mat, nv = K, nu = K)}
      #  print(paste('Dimension of mat:', dim(mat)))
      ncomp = sum(eigenmat$d > thresh)
      #print(ncomp)
      # print(eigenmat$d)
      if (ncomp < 2) {
        sqmat = as.matrix(eigenmat$v[,1]) %*% sqrt(eigenmat$d[1]) %*% t(as.matrix(eigenmat$v[,1]))
      }else{sqmat = eigenmat$v[,1:ncomp] %*% diag(sqrt(eigenmat$d[1:ncomp])) %*% t(eigenmat$v[,1:ncomp])
      }

      return(sqmat)
    }

    ginv.largep <- function(x.c, sqrtmat = TRUE, sqrtinvmat = TRUE, thresh = 10^{-20}){
      xxt.inv = MASS::ginv(x.c %*% t(x.c))
      tmp = xxt.inv %*% x.c
      sqrt.mat = sqrt.invmat = NULL
      if (sqrtinvmat == TRUE) {
        sqrt.mat = t(sqrtmat.comp(xxt.inv, thresh = thresh) %*% x.c) %*% x.c
      }
      if (sqrtinvmat == TRUE) {
        sqrt.invmat = t(sqrtmat.comp(xxt.inv,thresh = thresh) %*% x.c) %*% xxt.inv %*% x.c
      }
      return(list(inv = t(tmp) %*% tmp, sqrtinv = sqrt.invmat, sqrtmat = sqrt.mat))
    }

    tUU = ginv.largep(U, sqrtmat = TRUE, sqrtinvmat = TRUE)



    ## Iterative Update
    lam1 = rep(sort(lambda1,decreasing = TRUE), each = length(lambda2))
    lam2 = rep(lambda2, length(lambda1))

    myfunc <- function(j,k, #~~~in our case j=1,b/c lam1,lam2 is just a # ; k=1 b/c just 1 alpha value
                       gamma_new = rep(0, V + max(0,ncol(I_update)) + 1),
                       alpha_new = rep(0,V)){

      if (verbose == TRUE) {print(paste("Lambda1=",lam1[j]))}


      iter = 0
      err = 1000
      allzero.count = 0
      sigma2penalty = matrix(1,V,V); diag(sigma2penalty) <- 0

      while (err > tol & iter < max.iter & allzero.count < 4) {

        alpha_old = alpha_new
        gamma_old = gamma_new
        beta_old = c(gamma_old, alpha_old)

        sigma1 = mean((Y - U %*% gamma_old)^2)
        tmp = M - matrix(X, N, 1) %*% matrix(alpha_old, 1, V)
        Sigma2 = t(tmp) %*% tmp/N

        Omega = QUIC(Sigma2, rho = sigma2penalty*lam2[j], msg = 0)#Inverse matrix of the covariance matrix of M
        Omega.sqrtmat = try(t(base::chol(Omega$X)), TRUE)
        if (is.matrix(Omega.sqrtmat) == FALSE) {
          tmp.omega.1 = base::chol(Omega$X, pivot = TRUE)
          Omega.sqrtmat = t(tmp.omega.1[,order(attr(tmp.omega.1, 'pivot'))])
        }

        #sqrtmat.comp(Omega$X)
        Omega.sqrtmat.inv = try(t(base::chol(Omega$W)), TRUE) #sqrtmat.comp(Omega$W)
        if (is.matrix(Omega.sqrtmat.inv) == FALSE) {
          tmp.omega.2 = base::chol(Omega$W,pivot = TRUE)
          Omega.sqrtmat.inv = t(tmp.omega.2[,order(attr(tmp.omega.2, 'pivot'))])
        }

        Asqmat = Matrix::bdiag(1/sqrt(sigma1) * tUU$sqrtmat, sqrt(as.numeric(tXX)) * Omega.sqrtmat)
        Asqmat.inv = Matrix::bdiag(sqrt(sigma1) * tUU$sqrtinv, 1/sqrt(as.numeric(tXX)) * Omega.sqrtmat.inv)
        C = Asqmat.inv %*% rbind(tUY/sigma1, Omega$X %*% tMX)

        fit = glmnet(as.matrix(Asqmat), as.matrix(C),
                     lambda = lam1[j],
                     alpha = 1,
                     penalty.factor = penalty.factor)

        beta_new = fit$beta
        if (all(beta_new[-1] == 0)) {allzero.count = allzero.count + 1}
        gamma_new = beta_new[1:(V + max(0,ncol(I_update)) + 1)] #this includes c,b1,b2
        alpha_new = beta_new[c(1:V) + max(0,ncol(I_update)) + V + 1] #this includes a

        err = sqrt(sum((beta_old[-1] - c(gamma_new[-1],alpha_new))^2))
        iter = iter + 1
        if (verbose == TRUE) {print(c(iter, err))}

      }



      return(list(betahat = beta_new,gammahat = gamma_new,alphahat = alpha_new,
                  Omegahat = Omega$X,sigmasq = sigma1,
                  lambda = lam1[j]))
    }


    zzz <- list()

    ## when the algorithm selects too many parameters, we stop there.
    for (k in 1:length(alpha)) {
      zzz[[k]] <- list()
      j = 0
      nonzeros = 0
      while (j < length(lam1) & nonzeros < non.zeros.stop) {
        j = j + 1
        re <- c();
        zzz[[k]][[j]] <- NULL
        gamma_init = rep(0,V + max(0,ncol(I_update)) + 1)
        alpha_init = rep(0,V)
        if (j > 1) {if (is.null(zzz[[k]][[j - 1]]) == FALSE) {
          gamma_init = zzz[[k]][[j - 1]]$betahat[1:(max(0,ncol(I_update)) + V + 1)]
          alpha_init = zzz[[k]][[j - 1]]$betahat[1:(V) + (max(0,ncol(I_update)) + V + 1)]
        }}

        try(re <- myfunc(j = j,k = k,gamma_new = gamma_init,alpha_new = alpha_init))
        zzz[[k]][[j]] <- re
        nonzeros = sum(re$betahat != 0)
      }
    }

    betaest = do.call(cbind,lapply(zzz,function(x0){
      do.call(cbind,lapply(x0,function(xx){
        xx$betahat
      }))
    }))
    lam1s = unlist(lapply(zzz,function(x0){unlist(lapply(x0,function(xx){xx$lambda}))}))
    sigmasqs = unlist(lapply(zzz,function(x0){unlist(lapply(x0,function(xx){xx$sigmasq}))}))

    cest = betaest[1,]
    nump = apply(as.matrix(betaest),2,function(x){sum(abs(x) > 0)})
    names(nump) = NULL

    if (Omega.out == FALSE) {
      Omegas = NULL
    } else {
      Omegas = lapply(zzz,function(x0){lapply(x0,function(xx){xx$Omegahat})
      })}

    hatb1 = betaest[(1:V) + 1,]*Y.sd/M.sd
    names(hatb1) = colnames(M)

    hatb2 = betaest[-c(1:(V + 1),(nrow(betaest) - V + 1):nrow(betaest)),]
    names(hatb2) = gsub("\\D", "", colnames(I_update))
    interaction = rep(0,ncol(M))
    names(interaction) = 1:ncol(M)
    interaction[which(names(interaction) %in% names(hatb2))] <- hatb2
    hatb2 = interaction*Y.sd/M.sd
    names(hatb2) = paste0("M",1:ncol(M))

    hata = betaest[(1:V) + max(0,ncol(I_update)) + V + 1,]*M.sd/X.sd
    names(hata) = colnames(M)

    return(list(
      c = cest,
      hatb1 = hatb1,
      hatb2 = hatb2,
      hata = hata,
      nump = nump,
      Omega = Omegas,
      sigmasq = sigmasqs
    ))
  }
