# code-utils-comparison-methods-linear.R

library(freebird)
#' Run Zhou et al. (2020)'s high-dimensional mediation test via `freebird::hilma()`
#'
#' @param X The n by q exposure matrix.
#' @param Y The n-dimensional outcome vector.
#' @param M The n by p mediator matrix (standardized internally via `scale()`).
#' @param S The n by s confounding variables matrix, or NULL if no confounders.
#' @return A list with:
#'   - ZWZbhat: estimated overall indirect (beta) effect for the first exposure
#'   - ZWZ_var_beta: estimated variance of ZWZbhat
#'   - ZWZ_Sn: test statistic for the indirect effect
#'   - ZWZ_alpha1: estimated direct effect for the first exposure
#'   - ZWZ_var_alpha1: estimated variance of ZWZ_alpha1
#'   - ZWZ_Tn: test statistic for the direct effect
#'   - ZWZ_Sn_pvalue: p-value for the indirect effect test
#'   - ZWZ_Tn_pvalue: p-value for the direct effect test
#'   - ZWZ_time: elapsed run time in seconds
Zhou2020 <- function(X, Y, M, S=NULL){
  ZWZ_start = Sys.time()
  # Tuning parameters lam_list should change base on different settings c1 or c2 
  # to achieve the best result
  ZWZ = hilma(Y,scale(M),cbind(X,S),mediation_setting='incomplete')
  ZWZ_time = unclass(difftime(Sys.time(),ZWZ_start, units = "secs"))[1]
  return(list(ZWZbhat = ZWZ$beta_hat[1],
              ZWZ_var_beta = ZWZ$sigma_beta_hat[1,1],
              ZWZ_Sn = ZWZ$teststat_beta[1],
              ZWZ_alpha1 = ZWZ$alpha1_hat[1], 
              ZWZ_var_alpha1 = ZWZ$sigma_alpha1_hat[1,1],
              ZWZ_Tn = ZWZ$teststat_alpha1[1],
              ZWZ_Sn_pvalue = ZWZ$pvalue_beta_hat[1], 
              ZWZ_Tn_pvalue = ZWZ$pvalue_alpha1_hat[1],
              ZWZ_time = ZWZ_time))
}


library(glmnet)
library(stats)

#' Guo et al. (2022, JoE; 2022, JASA)'s high-dimensional mediation test for
#' a continuous (gaussian) outcome
#'
#' Selects mediators via SCAD-penalized regression (tuned by HBIC over
#' `lamb_grid`/`lamb_grid0`) and then performs de-biased inference on the
#' indirect (beta) and direct (alpha1) effects.
#'
#' @param X The n by q exposure matrix.
#' @param Y The n-dimensional outcome vector; must be continuous (not
#'   binary/count), else the function stops/warns.
#' @param M The n by p mediator matrix; columns must have nonzero variance.
#' @param Z The n by s confounding variables matrix, or NULL if no confounders.
#' @param Y_family Outcome family; only "gaussian" is supported here.
#' @param scale Logical; if TRUE, X, M, Z are standardized and Y is centered
#'   before estimation.
#' @param lamb_grid Grid of SCAD tuning parameters used to select active
#'   mediators (model with X, M, and Z).
#' @param lamb_grid0 Grid of SCAD tuning parameters used for the reduced model
#'   (Y regressed on M and Z only, i.e. without X) used to build the LRT
#'   reference model.
#' @return A list with:
#'   - stat_HDGMM: test statistic for the indirect (mediation) effect
#'   - pval_HDGMM: p-value for the indirect effect (chi-square with q df)
#'   - beta_hat: estimated indirect effect(s) (0 if no mediators selected)
HDGMM_linear <- function(X, Y, M, Z = NULL,
                         Y_family = "gaussian",
                         scale = TRUE,
                         lamb_grid = seq(0.05, 1, length.out=20),
                         lamb_grid0 = lamb_grid) {
  
  Y_family <- match.arg(Y_family)
  
  q <- ncol(X)
  pp <- ncol(M)
  n  <- nrow(M)
  
  # sanity check: Y should be continuous
  if (!is.numeric(Y)) {
    stop("`Y` must be a numeric vector for `HDGMM()`.")
  }
  
  # check for binary outcome
  y_unique <- unique(Y)
  if (length(y_unique) <= 2) {
    stop("`Y` appears to be binary. Use the logistic version of HDGMM instead.")
  }
  
  # check for count outcome (all nonnegative integers)
  if (all(Y >= 0) && all(abs(Y - round(Y)) < .Machine$double.eps^0.5)) {
    warning("`Y` appears to be a count variable. Consider using the Poisson version of HDGMM.")
  }
  
  
  # sanity check: design matrix (X, Z) should not be singular
  if (is.null(Z)) {
    W <- X
    s <- 0
  } else {
    W <- cbind(X, Z)
    s <- ncol(Z)
  }
  
  qrW <- qr(W)
  
  if (qrW$rank < ncol(W)) {
    stop(
      "The design matrix formed by `X` and `Z` is singular (contains collinear columns). ",
      "Please remove redundant variables."
    )
  }
  
  # sanity check: mediators should not contain constant columns
  m_var <- apply(M, 2, stats::var)
  
  if (any(m_var == 0)) {
    const_cols <- which(m_var == 0)
    
    stop(
      "The mediator matrix `M` contains constant columns. ",
      "Columns with zero variance: ",
      paste(const_cols, collapse = ", "),
      ". Please remove these mediators."
    )
  }
  
  
  
  if (scale) {
    X <- base::scale(X)
    M <- base::scale(M)
    if (!is.null(Z)) {
      Z <- base::scale(Z)
    }
    Y <- Y - mean(Y)
  }
  
  S <- Z
  
  ## HDGMM
  hbic= c()
  results =lapply(lamb_grid, HBIC_calc, xx=X,yy=Y,mm=M, S=S, n_imp = 0)
  ngrid = length(lamb_grid)
  for( ii in 1: ngrid){
    hbic[ii] = results[[ii]]$BIC
  }
  
  id = which(hbic==min(hbic))
  id = utils::tail(id,1)
  result = results[[id]]
  alpha0_hat = result$alpha0
  alpha1_hat = result$alpha1
  alpha2_hat = result$alpha2
  
  hbic0= c()
  if(length(S) ==0){
    intcpt = matrix(rep(1, n), ncol = 1)
    results =lapply(lamb_grid0, HBIC_calc, xx=intcpt,yy=Y,mm=M,n_imp = 0)
  }else{
    results =lapply(lamb_grid0, HBIC_calc, xx=S,yy=Y,mm=M,n_imp = 0)
  }
  ngrid0 = length(lamb_grid0)
  for( ii in 1: ngrid0){
    hbic0[ii] = results[[ii]]$BIC
  }
  id = which(hbic0==min(hbic0))
  id = utils::tail(id,1)
  alpha0_tld = results[[id]]$alpha0
  alpha2_tld = results[[id]]$alpha1
  
  A = which(alpha0_hat!=0)
  A_tld = which(alpha0_tld!=0)
  
  
  if(length(A)==0)
  {
    stat_HDGMM <- 0
    pval_HDGMM <- 1-pchisq(stat_HDGMM, df=q)
    beta_hat <- 0
    output <- list(stat_HDGMM = stat_HDGMM,
                   pval_HDGMM = pval_HDGMM,
                   beta_hat = beta_hat)
  } else{
    
    M_A = as.matrix(M[,A])
    
    if(length(A_tld)==0)
    {
      M_tld = NULL
      alpha0_tld_input <- rep(0, length(A))
    }else{
      M_tld =  as.matrix(M[,A_tld])
      alpha0_tld_input <- alpha0_tld[A_tld]
    }
    
    
    # HDGMM
    HDGMM = mediation_inference_continuous(X, Y, M_A, S=S, M_tld = M_tld,
                                           alpha0_hat = alpha0_hat[A],
                                           alpha0_tld = alpha0_tld_input,
                                           alpha1_hat = alpha1_hat,
                                           alpha2_hat = alpha2_hat,
                                           alpha2_tld = alpha2_tld)
    
    # Output: HDGMM
    beta_hat <- HDGMM$beta_hat
    stat_HDGMM <- HDGMM$Sn
    pval_HDGMM <- 1-pchisq(stat_HDGMM, df=q)
    
    output <- list(stat_HDGMM = stat_HDGMM,
                   pval_HDGMM = pval_HDGMM,
                   beta_hat = beta_hat)
  }
  
  output
  
}


# Internal helper:
#' Compute de-biased test statistics and covariances for the indirect (beta)
#' and direct (alpha1) effects, given a set of selected mediators
#'
#' @param X The n by q exposure matrix.
#' @param Y The n-dimensional outcome vector.
#' @param M The n by p matrix of *selected* mediators (columns in active set A).
#' @param S The n by s confounding variables matrix, or NULL if no confounders.
#' @param M_tld The n by p_tld matrix of mediators selected under the reduced
#'   (no-X) model; if NULL, falls back to using M and re-estimating alpha0_tld.
#' @param alpha0_hat Estimated mediator coefficients (on M) from the full
#'   model; if NULL, all coefficients (alpha0/alpha1/alpha2) are re-estimated
#'   via OLS on Z = cbind(M, X[, S]).
#' @param alpha0_tld Estimated mediator coefficients from the reduced (no-X)
#'   model; if NULL, re-estimated via OLS of Y on MS = cbind(M_tld, S).
#' @param alpha1_hat Estimated exposure coefficients from the full model.
#' @param alpha2_hat Estimated confounder coefficients from the full model
#'   (unused if S is NULL).
#' @param alpha2_tld Estimated confounder coefficients from the reduced model
#'   (unused if S is NULL).
#' @return A list with:
#'   - Sn: Wald test statistic for the indirect effect beta_hat
#'   - Tn: likelihood-ratio-type test statistic for the direct effect alpha1
#'   - beta_hat: estimated indirect effect (total effect gamma_hat minus
#'     direct effect alpha1_hat)
#'   - alpha0_hat, alpha1_hat, alpha2_hat: (re-)estimated coefficients
#'   - B: bias-correction matrix used in the sandwich variance formulas
#'   - var_beta: estimated covariance matrix of beta_hat
#'   - var_alpha1_hat: estimated covariance matrix of alpha1_hat
mediation_inference_continuous <-function(X, Y, M, S = NULL, M_tld = NULL,
                                          alpha0_hat = NULL, alpha0_tld = NULL,
                                          alpha1_hat = NULL, alpha2_hat=NULL,
                                          alpha2_tld = NULL){
  p = ncol(M)
  q = ncol(X)
  n = nrow(X)
  if(is.null(S)){
    ## no confounders
    Z = cbind(M,X)
    s = 0
    if(is.null(M_tld)){MS = M }else{ MS = M_tld}
    if(is.null(alpha0_tld)){alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y}
    RSS02 = t(Y - MS%*% alpha0_tld) %*% (Y - MS%*%alpha0_tld)
  }
  else{
    ## confounders
    Z = cbind(M,X,S)
    s = ncol(S)
    if(is.null(M_tld)){MS = cbind(M,S) }else{ MS = cbind(M_tld,S)}
    if(is.null(alpha0_tld)){alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y}
    RSS02 = t(Y - MS%*% c(alpha0_tld, alpha2_tld)) %*% (Y - MS%*% c(alpha0_tld,alpha2_tld))
  }
  if(is.null(alpha0_hat)){
    alpha_rf = solve(t(Z)%*%Z)%*%t(Z)%*%Y
    alpha0_hat = alpha_rf[1:p]
    alpha1_hat = alpha_rf[(p+1):(p+q)]
    if(s >0){alpha2_hat = alpha_rf[(p+q+1):(p+q+s)]}
  }else{
    alpha_rf=c(alpha0_hat,alpha1_hat,alpha2_hat)
  }
  
  res = Y - Z%*%alpha_rf
  RSS12 = as.numeric(t(res) %*% (res)) # Test direct effect
  
  # RSS01 = t(Y - X%*%alpha1_hat) %*% (Y - X%*%alpha1_hat) # Test indirect effect
  # RSS11 = t(Y-X%*%gamma_hat) %*% (Y-X%*%gamma_hat)
  
  df = p+q + s
  sigma1_hat = RSS12/(n - df)
  Sigma_MM = t(M)%*%M /n
  if(s == 0){
    ## no confounders
    gamma_hat = solve(t(X)%*%X)%*%t(X)%*%Y
    sigmaT_hat = t(Y-X%*%gamma_hat) %*% (Y-X%*%gamma_hat)/(n-q)
    beta_hat = gamma_hat -alpha1_hat
    sigma2_hat = pmax(0,(sigmaT_hat - sigma1_hat))
    
    #tmp1 = cbind(t(X)%*%X,t(X)%*%M)
    #tmp2 = cbind(t(M)%*%X,t(M)%*%M )
    #Sigma_hat = rbind(tmp1,tmp2)/n
    invXX = solve(t(X)%*%X/n)
    Sigma_MX =t(M)%*%X/n
    
    B = invXX %*%t(Sigma_MX) %*%solve(Sigma_MM -  Sigma_MX%*%invXX %*% t(Sigma_MX)) %*%Sigma_MX %*%invXX
    var_alpha1_hat = sigma1_hat*(invXX + B)
    cov_beta_hat = sigma2_hat * invXX + sigma1_hat * B
  }else{
    ## confounders
    V = cbind(X,S)
    gamma_hat = solve(t(V)%*%V)%*%t(V)%*%Y
    sigmaT_hat = t(Y-V%*%gamma_hat) %*% (Y-V%*%gamma_hat)/(n-q-s)
    beta_hat = gamma_hat[1:q] -alpha1_hat
    sigma2_hat = pmax(0,(sigmaT_hat - sigma1_hat))
    Sigma_VV = t(V)%*%V/n
    invVV = solve(Sigma_VV)
    Sigma_MV =t(M)%*%V/n
    Sigma_VM =t(Sigma_MV)
    B = invVV %*% Sigma_VM %*%solve(Sigma_MM - Sigma_MV %*% invVV %*% Sigma_VM)%*%Sigma_MV%*%invVV
    var_alpha1_hat = sigma1_hat*(invVV + B)[1:q,1:q]
    cov_beta_hat = sigma2_hat * invVV + sigma1_hat * B
    cov_beta_hat = cov_beta_hat[1:q, 1:q]
  }
  # Test for beta
  # Wald's test
  Sn = n*t(beta_hat) %*% solve(cov_beta_hat) %*% beta_hat
  # Test for alpha1
  # LRT
  Tn1 = (n-df) * (RSS02-RSS12)/RSS12
  #Tn2 = n*log(RSS02/RSS12)
  return(list(Sn = Sn, Tn = Tn1,
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B,
              var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat))
}


# Internal helper:
#' @keywords internal
#' @noRd
deSCAD <- function(z,lamb,a=3.7){
  # First order derivative of SCAD penalty
  # tuning parameter "a" (or "gamma") use the default value 3.7
  #' @param z a coefficient value (or vector) at which to evaluate the
  #'   derivative of the SCAD penalty
  #' @param lamb the SCAD tuning parameter lambda
  #' @param a the SCAD shape parameter (default 3.7, as recommended in the
  #'   original SCAD literature)
  #' @return the SCAD derivative evaluated at z, used as a penalty.factor
  #'   weight in the local linear approximation (LLA) step of `LLA_h1`
  return(1*(z<=lamb)+pmax((a*lamb-z),0)/((a-1)*lamb)*(lamb<z))
}

# Internal helper:
#' @keywords internal
#' @noRd
#' Estimate mediator/exposure/confounder coefficients via a two-step local
#' linear approximation (LLA) to the SCAD penalty (Lasso warm start, then one
#' SCAD-reweighted Lasso refit)
#'
#' @param X The n by q exposure matrix.
#' @param Y The n-dimensional outcome vector.
#' @param M The n by p mediator matrix.
#' @param lamb The tuning parameter lambda for the SCAD/Lasso penalty.
#' @param n Sample size (number of rows of X/M/Y).
#' @param p Number of mediators (columns of M).
#' @param q Number of exposures (columns of X).
#' @param n_imp Number of "important" leading mediator columns to leave
#'   unpenalized (default 0, i.e. all mediators penalized).
#' @param S The n by s confounding variables matrix, or NULL if no confounders.
#' @return A numeric vector of length (p + q + s): the estimated coefficients
#'   on M, X, and S (in that order) from the second-step (SCAD-reweighted)
#'   Lasso fit.
LLA_h1<- function(X,Y,M,lamb,n,p,q, n_imp = 0,S = NULL){
  if(length(S) == 0){
    s = 0
    V = X
  }else{
    s = ncol(S)
    V = cbind(X,S)
  }
  # Step 1 using Lasso
  w1 = matrix(0,nrow = (p + q + s),ncol=1)
  w1[1:(p-n_imp)] = 1 ## do not penalize the last n_imp coefficients
  alpha_int = coef(glmnet(cbind(M,V),Y,family = 'gaussian', alpha=1,
                          lambda = lamb, penalty.factor=w1,intercept = FALSE))[-1]
  # Step 2 using linear approximation of SCAD
  w2 = matrix(0,nrow = (p+q + s),ncol=1)
  for(j in 1:(p- n_imp)){
    w2[j] = deSCAD(alpha_int[j],lamb)
  }
  alpha = coef(glmnet(cbind(M,V),Y,family = 'gaussian', alpha=1,
                      lambda = lamb, penalty.factor=w2, intercept = FALSE))
  return(alpha[-1])
}


# Internal helper:
HBIC_calc <- function(lamb, xx,yy,mm,S = NULL, n_imp = 0){
  # Calculate HBIC for a specific tuning parameter lambda
  #' @param lamb a float value of tuning parameter lambda
  #' @param xx The n by q exposure matrix. q can be 1, and q < n is required
  #' @param yy The n-dimensional outcome vector.
  #' @param mm The n by p mediator matrix. p can be larger than n.
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required.
  #' @param n_imp an int for important mediators that will not be penalized
  #' @return
  #'        A list of class `HBIC_calc'
  #'        - BIC: HBIC score
  #'        - alpha0: estimated alpha0,
  #'        - alpha1: estimated alpha1,
  #'        - alpha2: estimated alpha2,
  #'        - sigma1_hat: estimated sigma1
  #'
  n = nrow(xx)
  p = ncol(mm)
  q = ncol(xx)
  if(is.null(S)){
    s = 0
    result <- LLA_h1(xx,yy,mm,lamb,n,p,q, n_imp = n_imp)
    alpha0 = result[1:p]
    alpha1 = result[(p+1):(p+q)]
    alpha2 = NULL
    tmp = yy - mm%*%alpha0 - xx%*% alpha1
  }else{
    s = ncol(S)
    result <- LLA_h1(xx,yy,mm,lamb,n,p,q,n_imp = n_imp, S = S)
    alpha0 = result[1:p]
    alpha1 = result[(p+1):(p+q)]
    alpha2 = result[(p+q+1):(p+q+s)]
    
    tmp = yy - mm%*%alpha0 - xx%*% alpha1 - S %*% alpha2
  }
  
  df = length(which(alpha0!= 0))+q + s
  sigma_hat = t(tmp)%*%tmp/n
  BIC = log(sigma_hat) + df*log(log(n))*log(p+q + s)/n
  #obj = objective(xx,yy,M,alpha0,alpha1,lamb)
  return(list(BIC=BIC,alpha0=alpha0,alpha1 = alpha1,
              alpha2 = alpha2, sigma1_hat = sigma_hat))
}

