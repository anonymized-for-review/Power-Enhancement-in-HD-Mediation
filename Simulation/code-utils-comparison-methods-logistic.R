# code-utils-comparison-methods-logistic.R

# Djordjilovic et al. (2019)'s method - "globaltest" R package 
library(globaltest)
Vera2019_GLM <- function(X,Y,M, S=NULL){
  Vera_start = Sys.time()
  colnames(M) <- paste0("M", 1:ncol(M))
  n = dim(M)[1]
  if(is.null(S))
  {
    MY_test <- gt(response= Y, 
                  alternative = M, 
                  null = ~X, model = 'logistic')
    
    # Null model contain intercept
    XM_test <- gt(response = X,
                  alternative = M, model ='linear')
  }else{
    # Null model contain intercept and X
    MY_test <- gt(response= Y, 
                  alternative = M, 
                  null = ~X+S, model = 'logistic')
    
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
library(stats)
library(glmnet)
library(ncvreg)
HDGMM_logistic <- function(X, Y, M, Z = NULL,
                           Y_family = "binomial",
                           scale = TRUE,
                           lamb_grid = seq(0.05, 1, length.out=20),
                           lamb_grid0 = lamb_grid) {
  
  Y_family <- match.arg(Y_family)
  
  q <- ncol(X)
  pp <- ncol(M)
  p <- ncol(M)
  n  <- nrow(M)
  phi0 <- 1
  
  # sanity check: Y should be dichotomous (0/1)
  if (!is.numeric(Y)) {
    stop("`Y` must be a numeric vector for `PEHDGMM_logistic()`.")
  }
  
  y_unique <- unique(Y)
  y_unique <- y_unique[!is.na(y_unique)]
  
  if (length(y_unique) != 2) {
    stop("`Y` must be binary with two distinct values for `PEHDGMM_logistic()`.")
  }
  
  if (!all(y_unique %in% c(0, 1))) {
    stop(
      "`Y` must be coded as 0/1 for `PEHDGMM_logistic()`. ",
      "Detected values: ", paste(y_unique, collapse = ", ")
    )
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
  HBIC_rcd <- rep(0, ngrid)
  alpha_rcd <- matrix(0, ncol = (p + q + s), nrow = ngrid)
  # penalize mediators only; do not penalize X or S
  w <- rep(0, p + q + s)
  w[1:p] <- 1
  
  for (j in 1:ngrid){
    result = HBIC_bino(X=X,Y=Y,M=M, S=S,w=w, lamb=lamb_grid[j])
    HBIC_rcd[j] = result$HBIC
    alpha_rcd[j,] = result$alpha
    #print(which(result$alpha !=0))
  }
  # Select the lambda and estimated coefficients on smallest HBIC
  id = which(HBIC_rcd == min(HBIC_rcd))
  id = utils::tail(id,1)
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
    selected_mediators <- numeric(0)
    output <- list(stat_HDGMM = stat_HDGMM,
                   pval_HDGMM = pval_HDGMM,
                   beta_hat = beta_hat)
  }
  
  else{
    
    M_A = as.matrix(M[,A]) # SCAD-selected mediators
    
    
    # HDMM
    HDGMM = Testing_bino(X,Y,M_A,S,phi0,a0_hat=alpha0_hat[A],a1_hat=alpha1_hat,a2_hat=alpha2_hat)
    
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
SCAD <- function(z, l1, l2, gamma,v) {
  if (abs(z) <= l1)
  {return(0)}
  else if (abs(z) <= (l1*(1+l2)+l1))
  {return(sign(z)*(abs(z)-l1)/(v*(1+l2)))}
  else if (abs(z) <= gamma*l1*(1+l2))
  {return(sign(z)*(abs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)))}
  else
  {return(z/(v*(1+l2)))}
}

# Internal helper:
MCP <- function(z, l1, l2, gamma, v){
  if (abs(z) <= l1) return(0)
  else if (abs(z) <= gamma*l1*(1+l2)) return(sign(z)*(abs(z)-l1)/(v*(1+l2-1/gamma)))
  else return(z/(v*(1+l2)))
}


# Internal helper:
p_binomial <- function(eta) {
  if (eta > 10) {
    return(1)
  } else if (eta < -10) {
    return(0)
  } else {
    return(exp(eta)/(1+exp(eta)))
  }
}

# Internal helper:
dedeb <- function(z){
  return(diag(c(exp(z)/(1+exp(z))^2)))
}

# Internal helper:
MedReg <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"),
                   gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100,
                   lambda, eps=1e-4, max_iter=10000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(x)),
                   warn=TRUE, returnX, Intercept = FALSE){
  x = std(X)
  p = ncol(x)
  n = nrow(x)
  ns = attr(x, "nonsingular")
  # Coersion
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  #if (class(x) != "matrix") {
  #  tmp <- try(x <- model.matrix(~0+., data=x), silent=TRUE)
  #  if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  #}
  if (storage.mode(x)=="integer") storage.mode(x) <- "double"
  if (!is.numeric(y)) {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }
  if (storage.mode(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"
  L = length(lambda)
  a = matrix(0,nrow = 1,ncol = p)
  b = matrix(0,nrow = L,ncol = p)
  lam = lambda
  b0 = matrix(0,nrow = 1,ncol = L)
  Dev = matrix(0, nrow = 1, ncol = L)
  e1 = matrix(0,nrow = 1, ncol = p)
  e2 = e1
  Eta = matrix(0,nrow = n, ncol = L)
  eta = matrix(0,nrow = n,ncol = 1)
  iter = matrix(0,nrow = 1,ncol = L)
  tot_iter = 0
  r = matrix(0,nrow = n,ncol = 1)
  s = matrix(0,nrow = n,ncol = 1)
  w = matrix(0,nrow = n,ncol = 1)
  z = matrix(0,nrow = 1,ncol = p)
  user = as.integer(any(penalty.factor==0))
  if (family=="binomial") {
    ## Initialization
    ybar = mean(y)
    b0[1] = log(ybar/(1-ybar))
    a0 = b0[1]
    nullDev = sum( -y*log(ybar) - (1-y)*log(1-ybar))
    s = y - ybar
    eta = rep(a0,n)
    z = t(s) %*% x
    
    ## If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
    if (user) {
      lstart = 1
    } else {
      lstart = 2
      Dev[1] = nullDev
      Eta[1,] = eta
    }
    for (l in lstart:L) {
      if (l != 1) {
        ## Assign a, a0
        a0 = b0[l-1]
        a =  b[(l-1),]
        
        ## Check dfmax
        nv = length(which(a!=0))
        if ((nv > dfmax) | (tot_iter == max_iter)) {
          for (ll in l:L) iter[ll] = NA
          break
        }
        
        ## Determine eligible set
        if (penalty== "MCP") cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1])
        if (penalty== "SCAD") cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1])
        for (j in 1:p) if (abs(z[j]) > (cutoff * alpha * penalty.factor[j])) e2[j] = 1
      } else {
        
        ## Determine eligible set
        lmax = max(abs(z))
        if (penalty == "MCP") cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax)
        if (penalty == "SCAD") cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax)
        for (j in 1:p){
          if (abs(z[j]) > (cutoff * alpha *  penalty.factor[j])){
            e2[j] = 1
          }
        }
      }
      
      
      while (tot_iter < max_iter) {
        while (tot_iter < max_iter) {
          maxChange = 1e-3
          while ((tot_iter < max_iter) & (maxChange > eps)) {
            iter[l] = iter[l] + 1
            tot_iter = tot_iter + 1
            Dev[l] = 0
            for (i in 1:n) {
              pi = p_binomial(eta[i])
              w[i] = pmax(pi*(1-pi), 0.0001)
              s[i] = y[i] - pi
              r[i] = s[i]/w[i]
              if (y[i]==1) Dev[l] = Dev[l] - log(pi)
              if (y[i]==0) Dev[l] = Dev[l] - log(1-pi)
            }
            if (Dev[l]/nullDev < 0.01) {
              if (warn) warning("Model saturated; exiting...");
              tot_iter = max_iter
              break
            }
            
            ## Intercept
            if(Intercept){
              xwr = t(w)%*%r
              xwx = sum(w)
              b0[l] = xwr/xwx + a0
              for (i in 1:n) {
                si = b0[l] - a0
                r[i] =  r[i] -si
                eta[i] = eta[i]+ si
              }
              maxChange = abs(si)*xwx/n
            }else{
              maxChange = 0
            }
            ## Covariates
            for (j in 1:p) {
              if (e1[j]) {
                
                ## Calculate u, v
                xwr = sum(x[,j] * r * w)
                xwx = sum(x[,j]^2*w)
                u = xwr/n + (xwx/n)*a[j]
                v = xwx/n
                
                ## Update b_j
                l1 = lam[l] * penalty.factor[j] * alpha
                l2 = lam[l] * penalty.factor[j] * (1-alpha)
                if (penalty == "MCP") b[l,j] = MCP(u, l1, l2, gamma, v)
                if (penalty =="SCAD") b[l,j] = SCAD(u, l1, l2, gamma, v)
                
                ## Update r
                shift = b[l,j] - a[j]
                if (shift !=0) {
                  for (i in 1:n) {
                    si = shift*x[i,j]
                    r[i] =r[i]- si
                    eta[i] =eta[i]+ si
                  }
                  if (abs(shift)*sqrt(v) > maxChange) {
                    maxChange = abs(shift)*sqrt(v)}
                }
              }
            }
            
            ## Check for convergence
            a0 = b0[l]
            a = b[l,]
            
          }
          ## Scan for violations in strong set
          violations = 0
          for (j in 1:p) {
            if (e1[j]==0 & e2[j]==1) {
              z[j] = t(x[,j])%*%s/n
              l1 = lam[l] * penalty.factor[j] * alpha
              if (abs(z[j]) > l1) {
                e1[j] = 1
                e2[j] = 1
                violations = violations + 1
              }
            }
          }
          if (violations==0) break
        }
        
        ## Scan for violations in rest
        violations = 0
        for (j in 1:p) {
          if (e2[j]==0) {
            z[j] = t(x[,j])%*%s/n
            l1 = lam[l] * penalty.factor[j] * alpha
            if (abs(z[j]) > l1) {
              e1[j] = 1
              e2[j] = 1
              violations = violations+1
            }
          }
        }
        if (violations==0) {
          Eta[,l] = eta
          break
        }
      }
    }
    ## Unstandardize
    b =matrix(b, p, L,byrow=TRUE)
    
    if(Intercept){
      beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
      bb <- b/attr(x, "scale")[ns]
      beta[ns+1,] <- bb
      beta[1,] <- b0 - attr(x, "center")[ns]%*% bb
    }else{
      beta = b/attr(x, "scale")
    }
  }
  return(list(b = b,beta = beta,Dev=Dev, Eta = Eta, iter = iter))
}


# Internal helper:
HBIC_bino <- function(X,Y,M,S = NULL,w,lamb,
                      Y_family = "binomial",
                      penalty_type = "SCAD"){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
  
  if (is.null(S)) {
    s <- 0
    x <- as.matrix(cbind(M, X))
  } else {
    S <- as.matrix(S)
    if (nrow(S) != n) stop("S must have the same number of rows as X.")
    s <- ncol(S)
    x <- as.matrix(cbind(M, X, S))
  }
  
  #x = as.matrix(cbind(M,X))
  #p = ncol(x)
  #n = nrow(x)
  
  result <- MedReg(
    x, Y,
    family = Y_family,
    penalty = penalty_type,
    lambda = lamb,
    penalty.factor = w,
    alpha = 1,
    Intercept = FALSE
  )
  
  alpha <- as.vector(result$beta)
  
  alpha0 <- alpha[1:p]
  alpha1 <- alpha[(p + 1):(p + q)]
  alpha2 <- if (s > 0) alpha[(p + q + 1):(p + q + s)] else NULL
  
  logLkhd <- Ln(X, Y, M, alpha0, alpha1, S = S, alpha2 = alpha2)
  
  # mediator coefficients are sparse; X and S are always counted
  df <- sum(abs(alpha0) > 0) + q + s
  
  # p + q + s is the total number of regression coefficients
  bic <- -n * logLkhd + df * log(log(n)) * log(p + q + s)
  
  return(list(
    HBIC = bic,
    logLkhd = logLkhd,
    alpha = alpha,
    alpha0 = alpha0,
    alpha1 = alpha1,
    alpha2 = alpha2
  ))
}


# Internal helper:
Testing_bino <-function(X,Y,M_A, S = NULL, phi0=1,
                        a0_hat=NULL, a1_hat=NULL,a2_hat = NULL,
                        Y_family = "binomial"){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M_A)
  
  if (!is.null(S)) {
    s <- ncol(S)
  } else {
    s <- 0
  }
  
  # if(is.null(a0_hat)){
  # #Ocl = glm(Y~0 + M_A + X,family = binomial(link = "logit")) # Refit the model with selected subset of variable
  # Ocl = glm(Y~0 + M_A + X,family = Y_family) # Refit the model with selected subset of variable
  # a0 = Ocl$coefficients[1:p]
  # a1 = Ocl$coefficients[(p+1):(p+q)]
  # }
  # else{a0 = a0_hat; a1 = a1_hat}
  
  # Refit logistic model if coefficients are not supplied
  if (is.null(a0_hat) || is.null(a1_hat) || (!is.null(S) && is.null(a2_hat))) {
    dat <- data.frame(Y = Y, M_A, X)
    if (!is.null(S)) dat <- cbind(dat, S)
    
    # no intercept, consistent with your original setup
    Ocl <- glm(Y ~ 0 + ., family = Y_family, data = dat)
    
    cf <- coef(Ocl)
    a0 <- cf[1:p]
    a1 <- cf[(p + 1):(p + q)]
    a2 <- if (!is.null(S)) cf[(p + q + 1):(p + q + s)] else NULL
  } else {
    a0 <- a0_hat
    a1 <- a1_hat
    a2 <- if (!is.null(S)) a2_hat else NULL
  }
  
  
  if(is.null(S))
  {
    beta_hat = solve(t(X)%*%X) %*% t(X) %*% M_A %*% a0
    hatsigma2 = sum((M_A %*% a0 - X %*% beta_hat)^2)/(n-1)
    z = X %*% a1 + M_A %*% a0
    bbb = dedeb(z) # the second derivative of b(z)
    Sigmaxx = t(X) %*% bbb %*% X/n
    Sigmaxm = t(X) %*% bbb %*% M_A/n
    Sigmamm = t(M_A) %*% bbb %*% M_A/n
    tldSigmaxx = t(X)%*%X/n
    tldSigmaxm = t(X) %*%M_A/n
    Sigmaxx_inv = solve(Sigmaxx)
    Sigmamm.x = Sigmamm - t(Sigmaxm) %*% Sigmaxx_inv %*% Sigmaxm
    Sigmamm.x_inv = solve(Sigmamm.x)
    tldSigmaxx_inv = solve(tldSigmaxx)
    # Calculate the variance of \hat{alpha_1} and \hat{\beta_1}
    var_alpha = phi0*(Sigmaxx_inv +Sigmaxx_inv %*%Sigmaxm %*% Sigmamm.x_inv %*% t(Sigmaxm) %*%Sigmaxx_inv)
    var_beta = hatsigma2 * tldSigmaxx_inv + phi0* tldSigmaxx_inv %*% tldSigmaxm %*% Sigmamm.x_inv %*% t(tldSigmaxm) %*%tldSigmaxx_inv
    
  } else {
    
    # Residual-maker for S
    P_S <- S %*% solve(t(S) %*% S) %*% t(S)
    M_S <- diag(n) - P_S
    
    # beta_hat: coefficient of X in regression of M_A %*% a0 on (X, S)
    beta_hat <- solve(t(X) %*% M_S %*% X) %*% t(X) %*% M_S %*% M_A %*% a0
    
    # coefficient of S in same regression
    gamma_hat <- solve(t(S) %*% S) %*% t(S) %*% (M_A %*% a0 - X %*% beta_hat)
    
    # residual variance
    resid_vec <- M_A %*% a0 - X %*% beta_hat - S %*% gamma_hat
    hatsigma2 <- sum(resid_vec^2) / (n - ncol(X) - ncol(S))
    
    # logistic linear predictor
    z <- X %*% a1 + M_A %*% a0 + S %*% a2
    
    # second derivative matrix
    bbb <- dedeb(z)
    
    # weighted blocks
    Sigmaxx <- t(X)   %*% bbb %*% X   / n
    Sigmaxs <- t(X)   %*% bbb %*% S   / n
    Sigmaxm <- t(X)   %*% bbb %*% M_A / n
    
    Sigmass <- t(S)   %*% bbb %*% S   / n
    Sigmasm <- t(S)   %*% bbb %*% M_A / n
    
    Sigmamm <- t(M_A) %*% bbb %*% M_A / n
    
    # unweighted partialed-out matrices
    tldSigmaxx.s <- t(X) %*% M_S %*% X   / n
    tldSigmaxm.s <- t(X) %*% M_S %*% M_A / n
    tldSigmaxx.s_inv <- solve(tldSigmaxx.s)
    
    # block matrices for (X, S)
    Sigma_xs_xs <- rbind(
      cbind(Sigmaxx, Sigmaxs),
      cbind(t(Sigmaxs), Sigmass)
    )
    
    Sigma_xs_m <- rbind(
      Sigmaxm,
      Sigmasm
    )
    
    # mediator block adjusted for X and S
    Sigmamm.xs <- Sigmamm - t(Sigma_xs_m) %*% solve(Sigma_xs_xs) %*% Sigma_xs_m
    Sigmamm.xs_inv <- solve(Sigmamm.xs)
    
    # X block adjusted for S
    Sigmaxx.s <- Sigmaxx - Sigmaxs %*% solve(Sigmass) %*% t(Sigmaxs)
    Sigmaxm.s <- Sigmaxm - Sigmaxs %*% solve(Sigmass) %*% Sigmasm
    Sigmaxx.s_inv <- solve(Sigmaxx.s)
    
    # variance of alpha1
    var_alpha <- phi0 * (
      Sigmaxx.s_inv +
        Sigmaxx.s_inv %*% Sigmaxm.s %*% Sigmamm.xs_inv %*%
        t(Sigmaxm.s) %*% Sigmaxx.s_inv
    )
    
    # variance of beta_hat
    var_beta <- hatsigma2 * tldSigmaxx.s_inv +
      phi0 * tldSigmaxx.s_inv %*% tldSigmaxm.s %*% Sigmamm.xs_inv %*%
      t(tldSigmaxm.s) %*% tldSigmaxx.s_inv
  }
  
  
  Sn <- n * t(beta_hat) %*% solve(var_beta) %*% beta_hat
  
  # constrained model: drop X, keep M_A and S
  if (is.null(S)) {
    refit_tld <- glm(Y ~ 0 + M_A, family = Y_family)
    alpha0_tld <- as.matrix(coef(refit_tld))
    alpha2_tld <- NULL
    Tn <- 2 * n * (Ln(X, Y, M_A, a0, a1) -
                     Ln(X, Y, M_A, alpha0_tld, as.matrix(integer(q)))) / phi0
  } else {
    dat_tld <- data.frame(Y = Y, M_A, S)
    refit_tld <- glm(Y ~ 0 + ., family = Y_family, data = dat_tld)
    cf_tld <- coef(refit_tld)
    alpha0_tld <- as.matrix(cf_tld[1:p])
    alpha2_tld <- as.matrix(cf_tld[(p + 1):(p + s)])
    
    Tn <- 2 * n * (
      Ln(X, Y, M_A, a0, a1, S = S, alpha2 = a2) -
        Ln(X, Y, M_A, alpha0_tld, as.matrix(integer(q)), S = S, alpha2 = alpha2_tld)
    ) / phi0
  }
  
  # # Test of direct effect likelihood ratio test
  # Tn = 2*n*(Ln(X,Y,M_A,a0,a1) - Ln(X,Y,M_A,alpha0_tld,as.matrix(integer(q))))/phi0
  # #Tn2 = 2 * n * (logLik(Ocl) - logLik(refit_tld))/phi0
  # return(list(beta_hat = beta_hat, alpha0_hat = a0, alpha1_hat = a1,
  #             var_beta = var_beta, var_alpha1 = var_alpha,
  #             Sn = Sn, Tn = Tn, sigma2 = hatsigma2))
  #
  
  return(list(
    beta_hat = beta_hat,
    alpha0_hat = a0,
    alpha1_hat = a1,
    alpha2_hat = a2,
    var_beta = var_beta,
    var_alpha1 = var_alpha,
    Sn = Sn,
    Tn = Tn,
    sigma2 = hatsigma2
  ))
}

# Internal helper:
Ln <- function(X, Y, M, alpha0, alpha1, S = NULL, alpha2 = NULL) {
  eta <- (M%*%alpha0 + X%*% alpha1)
  if (!is.null(S)) {
    eta <- eta + S %*% alpha2
  }
  return(mean(diag(Y)%*% eta - log(1 + exp(eta))))
}
