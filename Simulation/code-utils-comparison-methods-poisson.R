# code-utils-comparison-methods-poisson.R

# Djordjilovic et al. (2019)'s method - "globaltest" R package 
library(globaltest)
Vera2019_Poisson <- function(X,Y,M, S=NULL){
  Vera_start = Sys.time()
  colnames(M) <- paste0("M", 1:ncol(M))
  n = dim(M)[1]
  
  if(is.null(S))
  {
    MY_test <- gt(response= Y, 
                  alternative = M, 
                  null = ~X, model = 'poisson')
    
    # Null model contain intercept
    XM_test <- gt(response = X,
                  alternative = M, model ='linear')
  }else{
    # Null model contain intercept and X
    MY_test <- gt(response= Y, 
                  alternative = M, 
                  null = ~X+S, model = 'poisson')
    
    # Null model contain intercept
    XM_test <- gt(response = X,
                  alternative = M, 
                  null = ~S, model ='linear')
    
  }
  pval_Vera = pmax(MY_test@result[1], XM_test@result[1])
  Vera_time = unclass(difftime(Sys.time(),Vera_start, units = "secs"))[1]
  return(list(pval_Vera=pval_Vera, Vera_time=Vera_time))
}


## Guo et al. (2023, JBES)
library(glmnet)
library(stats)
library(ncvreg)
HDGMM_poisson <- function(X, Y, M, Z = NULL,
                          Y_family = "poisson",
                          scale = TRUE,
                          lamb_grid = seq(0.05, 1, length.out=20),
                          lamb_grid0 = lamb_grid) {
  
  Y_family <- match.arg(Y_family)
  
  q <- ncol(X)
  pp <- ncol(M)
  p <- ncol(M)
  n  <- nrow(M)
  phi0 <- 1
  
  # sanity check: Y should be a count outcome for Poisson model
  if (!is.numeric(Y)) {
    stop("`Y` must be a numeric vector for `HDGMM_poisson()`.")
  }
  
  y_unique <- unique(Y)
  y_unique <- y_unique[!is.na(y_unique)]
  
  # check nonnegative
  if (any(Y < 0, na.rm = TRUE)) {
    stop("`Y` must be nonnegative for `HDGMM_poisson()`.")
  }
  
  # check integer-valued (within tolerance)
  tol <- sqrt(.Machine$double.eps)
  if (any(abs(Y - round(Y)) > tol, na.rm = TRUE)) {
    stop("`Y` must contain integer counts for `HDGMM_poisson()`.")
  }
  
  # check variation (not all identical)
  if (length(y_unique) == 1) {
    stop("`Y` has no variation (all values are identical). Poisson model is not identifiable.")
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
    Y <- Y
  }
  
  S <- Z
  
  ## HDMM
  ngrid = length(lamb_grid)
  HBIC_rcd= rep(0,ngrid)
  alpha_rcd = matrix(0,ncol = (p+q+s), nrow = ngrid)
  for (j in 1:ngrid){
    #HBIC_poi = HBIC_poisson(X,Y,M,lamb_grid[j],S)
    HBIC_poi = HBIC_poisson_V2(X,Y,M,lamb_grid[j],S)
    HBIC_rcd[j] = HBIC_poi$HBIC
    alpha_rcd[j,] = HBIC_poi$alpha_hat
    #print(which(HBIC_poi$alpha_hat!=0))
  }
  id = which(HBIC_rcd == min(HBIC_rcd))[1]
  lamb = lamb_grid[id]
  alpha_hat = alpha_rcd[id,]
  alpha0_hat = alpha_hat[1:p]
  alpha1_hat = alpha_hat[(p+1):(p+q)]
  alpha2_hat <- if (s > 0) alpha_hat[(p + q + 1):(p + q + s)] else NULL
  A = which(alpha0_hat!=0)
  
  
  if(length(A)==0)
  {
    stat_HDGMM <- 0
    pval_HDGMM <- 1-pchisq(stat_HDGMM, df=q)
    beta_hat <- 0
    output <- list(stat_HDGMM = stat_HDGMM,
                   pval_HDGMM = pval_HDGMM,
                   beta_hat = beta_hat)
  }
  
  else{
    
    M_A = as.matrix(M[,A]) # SCAD-selected mediators
    
    
    # HDMM
    HDGMM = Testing_poisson(X,Y,M_A,S,phi0,
                            alpha0_hat=alpha0_hat[A],
                            alpha1_hat=alpha1_hat,
                            alpha2_hat=alpha2_hat)
    
    # Output: HDMM
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
dedeb_poi<-function(z){
  return(exp(z))
}


# Internal helper:
deSCAD <- function(z,lamb,a=3.7){
  return(1*(z<=lamb)+pmax((a*lamb-z),0)/((a-1)*lamb)*(lamb<z))
}



# Internal helper:
ZouAlgo_h1 <- function(X,Y,M,w,lamb){
  X = as.matrix(X) # unpenalized block: can be X or cbind(X, S)
  M = as.matrix(M) # penalized mediator block
  n = nrow(X)
  p = ncol(M)
  q = ncol(X)
  alpha_int = matrix(NA,ncol=1,nrow=(p+q))
  U = which(w == 0)
  V = which(w!=0)
  Xt = sqrt(2)*cbind(M,X)
  Xts = Xt
  for(j in 1:length(V)){
    Xts[,V[j]] = Xt[,V[j]] * lamb/w[V[j]]
  }
  Xus = as.matrix(Xts[,U])
  Xvs = as.matrix(Xts[,V])
  Pu = Xus%*%solve(t(Xus)%*%Xus)%*%t(Xus)
  Qu = diag(n) - Pu
  Ys = sqrt(2)*Y
  Yss1 = sqrt(2)*Qu%*%Y
  Xvss1 = Qu%*%Xvs
  reg1 = glmnet(Xvss1,Yss1,family = "gaussian",alpha=1,lambda=lamb)
  Betavs = matrix(reg1$beta,ncol=1)
  Betaus = solve(t(Xus)%*%Xus)%*%t(Xus)%*%(Ys - Xvs%*%Betavs)
  alpha_int[U] = Betaus
  alpha_int[V] = Betavs*lamb/w[V]
  return(alpha_int)
}



# Internal helper:
poisson_h1 <- function(X,Y,M,lamb,S=NULL){
  p = ncol(M)
  q = ncol(X)
  n = nrow(M)
  
  if(is.null(S)){
    s = 0
    Z = as.matrix(cbind(M,X))
  }else{
    s = ncol(S)
    Z = as.matrix(cbind(M,X,S))
  }
  
  # Step 1 using Lasso
  # penalize M, do not penalize X or S
  w1 = vector(mode = "double",length= (p+q+s))
  w1[1:p] = 1
  res1 = glmnet(Z, Y,family = "poisson", alpha=1, lambda=lamb,penalty.factor = w1)
  alpha_int = as.matrix(res1$beta)
  
  # Step 2 using linear approximation of SCAD
  w2 =vector(mode = "double",length= (p+q+s))
  for(j in 1:p){
    w2[j] = deSCAD(alpha_int[j],lamb)
  }
  xs = Z
  ys = Y
  for (i in 1:nrow(M)){
    b_dede = sqrt(as.numeric(exp(t(Z[i,])%*%alpha_int)))
    xs[i,] = b_dede%*%Z[i,]
    ys[i] = b_dede* t(Z[i,])%*% alpha_int + (Y[i] - exp(t(Z[i,])%*%alpha_int))/b_dede
  }
  alpha = ZouAlgo_h1(xs[,(p+1):(p+q+s)],ys, xs[,1:p],w2,lamb)
  return(alpha)
}

# Internal helper:
HBIC_poisson_V2 <-function(X,Y,M,lamb,S = NULL){
  M = as.matrix(M)
  X = as.matrix(X)
  p = ncol(M)
  q = ncol(X)
  n = nrow(M)
  #Cn = log(n*log(n))/2
  if(is.null(S))
  {
    s <- 0
    Z = as.matrix(cbind(M,X))
    alpha_hat = poisson_h1(X,Y,M,lamb)
  }else{
    s <- ncol(S)
    Z <- as.matrix(cbind(M, X, S))
    alpha_hat <- poisson_h1(X, Y, M, lamb, S = S)
  }
  
  eta <- Z%*%alpha_hat
  ETA_MAX <- log(.Machine$double.xmax) - 1
  eta_safe <- pmin(eta, ETA_MAX)
  logLkhd <- mean(Y*eta_safe - exp(eta_safe))
  
  tol <- 10e-12
  alpha0_hat <- alpha_hat[1:p, , drop = FALSE]
  df <- sum(abs(alpha0_hat) > tol) + q + s
  
  Cn = log(log(n))
  
  HBIC <- -logLkhd + df * Cn * log(p + q + s) / n
  
  return(list(HBIC=HBIC,alpha_hat=alpha_hat))
}

# Internal helper:
Testing_poisson <- function(X,Y,M_A,S=NULL, phi0=1,
                            alpha0_hat=NULL, alpha1_hat=NULL, alpha2_hat=NULL){
  Y_family = "poisson"
  M_A = as.matrix(M_A)
  pm = ncol(M_A)
  q = ncol(X)
  n = nrow(M_A)
  
  if (!is.null(S)) {
    s <- ncol(S)
  } else {
    s <- 0
  }
  
  if (is.null(alpha0_hat) || is.null(alpha1_hat) || (!is.null(S) && is.null(alpha2_hat))) {
    if (is.null(S)) {
      refit <- glm(Y ~ 0 + M_A + X, family = Y_family)
      cf <- coef(refit)
      alpha0_hat <- as.matrix(cf[1:pm])
      alpha1_hat <- as.matrix(cf[(pm + 1):(pm + q)])
      alpha2_hat <- NULL
    } else {
      refit <- glm(Y ~ 0 + M_A + X + S, family = Y_family)
      cf <- coef(refit)
      alpha0_hat <- as.matrix(cf[1:pm])
      alpha1_hat <- as.matrix(cf[(pm + 1):(pm + q)])
      alpha2_hat <- as.matrix(cf[(pm + q + 1):(pm + q + s)])
    }
  } else {
    alpha0_hat <- as.matrix(alpha0_hat)
    alpha1_hat <- as.matrix(alpha1_hat)
    if (!is.null(S)) alpha2_hat <- as.matrix(alpha2_hat)
  }
  
  
  if(is.null(S))
  {
    
    beta_hat = solve(t(X)%*%X) %*% t(X) %*% M_A %*% alpha0_hat
    hatsigma2 = sum((M_A %*% alpha0_hat - X %*% beta_hat)^2)/(n-q)
    z = X %*% alpha1_hat + M_A %*% alpha0_hat
    bbb = diag(c(dedeb_poi(z))) # the second derivative of b(z)
    Sigmaxx = t(X) %*% bbb %*% X/n
    Sigmaxm = t(X) %*% bbb %*% M_A/n
    Sigmamm = t(M_A) %*% bbb %*% M_A/n
    tldSigmaxx = t(X)%*%X/n
    tldSigmaxm = t(X) %*%M_A/n
    Sigmaxx_inv = solve(Sigmaxx)
    Sigmamm.x = Sigmamm - t(Sigmaxm) %*% Sigmaxx_inv %*% Sigmaxm
    Sigmamm.x_inv = solve(Sigmamm.x)
    tldSigmaxx_inv = solve(tldSigmaxx)
    var_alpha = phi0*(Sigmaxx_inv +Sigmaxx_inv %*%Sigmaxm %*% Sigmamm.x_inv %*% t(Sigmaxm) %*%Sigmaxx_inv)
    var_beta = hatsigma2 * tldSigmaxx_inv + phi0* tldSigmaxx_inv %*% tldSigmaxm %*% Sigmamm.x_inv %*% t(tldSigmaxm) %*%tldSigmaxx_inv
    
  } else {
    P_S <- S %*% solve(t(S) %*% S) %*% t(S)
    M_S <- diag(n) - P_S
    
    beta_hat <- solve(t(X) %*% M_S %*% X) %*% t(X) %*% M_S %*% M_A %*% alpha0_hat
    
    gamma_hat <- solve(t(S) %*% S) %*% t(S) %*% (M_A %*% alpha0_hat - X %*% beta_hat)
    resid_vec <- M_A %*% alpha0_hat - X %*% beta_hat - S %*% gamma_hat
    hatsigma2 <- sum(resid_vec^2) / (n - q - s)
    
    z = X %*% alpha1_hat + M_A %*% alpha0_hat + S %*% alpha2_hat
    bbb = diag(c(dedeb_poi(z)))
    
    tldSigmaxx <- t(X) %*% M_S %*% X / n
    tldSigmaxm <- t(X) %*% M_S %*% M_A / n
    z <- X %*% alpha1_hat + M_A %*% alpha0_hat + S %*% alpha2_hat
    Sigmaxx <- t(X) %*% bbb %*% X / n
    Sigmaxm <- t(X) %*% bbb %*% M_A / n
    Sigmamm <- t(M_A) %*% bbb %*% M_A / n
    tldSigmaxx_inv <- solve(tldSigmaxx)
    Sigmaxs <- t(X) %*% bbb %*% S / n
    Sigmass <- t(S) %*% bbb %*% S / n
    Sigmasm <- t(S) %*% bbb %*% M_A / n
    Sigma_xs_xs <- rbind(cbind(Sigmaxx, Sigmaxs),cbind(t(Sigmaxs), Sigmass))
    Sigma_xs_m <- rbind(Sigmaxm,Sigmasm)
    
    Sigmamm.xs <- Sigmamm - t(Sigma_xs_m) %*% solve(Sigma_xs_xs) %*% Sigma_xs_m
    Sigmamm.xs_inv <- solve(Sigmamm.xs)
    
    ## X block adjusted for S
    Sigmaxx.s <- Sigmaxx - Sigmaxs %*% solve(Sigmass) %*% t(Sigmaxs)
    Sigmaxm.s <- Sigmaxm - Sigmaxs %*% solve(Sigmass) %*% Sigmasm
    Sigmaxx.s_inv <- solve(Sigmaxx.s)
    
    var_alpha <- phi0 * (
      Sigmaxx.s_inv +
        Sigmaxx.s_inv %*% Sigmaxm.s %*% Sigmamm.xs_inv %*%
        t(Sigmaxm.s) %*% Sigmaxx.s_inv
    )
    
    var_beta <- hatsigma2 * tldSigmaxx_inv +
      phi0 * tldSigmaxx_inv %*% tldSigmaxm %*% Sigmamm.xs_inv %*%
      t(tldSigmaxm) %*% tldSigmaxx_inv
  }
  
  # Test on indirect effect
  Sn =n * t(beta_hat) %*% solve(var_beta) %*% beta_hat
  
  # Test on direct effect
  ## constrained model removes X but keeps M_A and S
  if (is.null(S)) {
    refit_tld = glm(Y ~ 0 + M_A,family = Y_family)
    alpha0_tld = as.matrix(coef(refit_tld))
    Tn = 2*n*(Ln_poi(X,Y,M_A,alpha0_hat,alpha1_hat) - Ln_poi(X,Y,M_A,alpha0_tld,as.matrix(integer(q))))/phi0
  } else {
    refit_tld <- glm(Y ~ 0 + M_A + S, family = Y_family)
    cf_tld <- coef(refit_tld)
    alpha0_tld <- as.matrix(cf_tld[1:pm])
    alpha2_tld <- as.matrix(cf_tld[(pm + 1):(pm + s)])
    Tn <- 2*n*(
      Ln_poi(X, Y, M_A, alpha0_hat, alpha1_hat, S = S, a2 = alpha2_hat) -
        Ln_poi(X, Y, M_A, alpha0_tld, as.matrix(integer(q)), S = S, a2 = alpha2_tld)
    ) / phi0
  }
  return(list(beta_hat = beta_hat, alpha0_hat = alpha0_hat, alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat,
              var_beta = var_beta, var_alpha1 = var_alpha, Sn = Sn, Tn = Tn, sigma2 = hatsigma2))
}

# Internal helper:
Ln_poi <- function(X,Y,M_A,a0,a1,S = NULL, a2= NULL){
  z = X %*% a1 + M_A %*% a0
  if (!is.null(S)) {
    z <- z + S %*% a2
  }
  return(mean(diag(c(Y))%*%z - exp(z)))
}
