##############################################################
##  Simulation studies of PE-HDGMM for continuous outcomes  ##
##############################################################

parid <- 21 # taking values in c(1:21)
setting_id <- 2 # taking values in {1,2}
includeCovariate <- TRUE # taking values in {TRUE, FALSE}
excludeComparison <- TRUE
# taking values in {TRUE, FALSE}
# logical flag controlling whether slow comparison methods are skipped

Y_family <- "gaussian"
penalty_type <- "SCAD"

library(POEM)

source("code-utils-comparison-methods-linear.R")

n <- 300
q <- 1
p <- 500
c2 <- 0.5
repnum <- 5

## Signal pattern setup 
{
  if(setting_id==1) {
  alpham <- c(1,0.8,0.6,0.4,0.2,rep(0,p-5)) 
  }
  if(setting_id==2) {
    alpham <- c(1,-0.5, 0.4, -0.3, rep(0,p-4))  
  }
}


gen_error<-function(N,p,rho){
  X = matrix(NA,N,p)
  X[,1] = rnorm(N)
  for(ii in 2:p){
    X[,ii] = rho*X[,(ii-1)] + sqrt(1-rho^2)*rnorm(N)
  }
  return(X)
}


tau <- c(0.1, 0.2, 0.3, 0.4, 0.5, rnorm(p-5, mean=0, sd=0.1))
alphax <- c2

# Setup HBIC for different 
ngrid = 20 # ngrid for lambda
# The range of tuning parameters, search the best HBIC on an equal space grid
# The grid range may be adjusted for different settings to achieve the best result
lamb_grid = seq(0.2,0.39,length.out = ngrid) # Under full model
lamb_grid0 = seq(0.27,0.46,length.out = ngrid) # Under reduce model H0: \alpha_1=0

output_summary_tab <- c() 
output_pvals <- c() 
c1_list <- seq(-1, 1, by=0.1)
for (c1 in c1_list[parid])
{
  print(paste0("c1:", c1))
  Gamma <- c1* tau
  beta <- t(Gamma) %*% alpham

  ans <- c()
  for(ite in 1:repnum)
  {
  print(ite)
    
  ## Generate samples: 
  ## X: a matrix of covariates, n by q
  ## Y: a vector of responses, n by 1
  ## M: a matrix of mediators, n by p
  ## S: a matrix of confounders, n by s 
  
  if(!includeCovariate)
  {
  ## no confounders
  s <- 0
  X <- matrix(rnorm(n), nrow=n, ncol=q) ## n by q
  M <- X%*%t(Gamma) +  gen_error(n,p,rho = 0.5) ## n by p
  Y <- as.vector(M %*% alpham + X%*% alphax + rnorm(n,0,sd=0.5)) ## n by 1
  S = NULL
  
  ## PEHDGMM
  output_PEHDMM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M),  
                                  scale=FALSE, lambda_grid = lamb_grid, 
                                  lambda_grid_reduced =lamb_grid0)
  pval_PEHDMM <- output_PEHDMM$value[output_PEHDMM$term == "pval_pe"]
  
  
  ## Guo et al. (2022, JoE; 2022, JASA)
  output_HDMM <- HDGMM_linear(scale(X), Y-mean(Y), scale(M),  
                              scale=FALSE,  lamb_grid=lamb_grid, lamb_grid0=lamb_grid0)
  pval_HDMM <- output_HDMM$pval_HDGMM
  
  ## Zhou et al. (2020)
  if(excludeComparison)
    {
    output_zhou <- NULL
    pval_zhou <- NA
    } else{
    output_zhou <- Zhou2020(scale(X), Y-mean(Y), scale(M))
    pval_zhou <- output_zhou$ZWZ_Sn_pvalue
    }
  
  
  } else {
  ## with confounders
  s <- 5
  X <- matrix(rnorm(n), nrow=n, ncol=q) ## n by q
  S = matrix(rnorm(n*s,mean=0,sd=1),ncol=s)
  alphaz = runif(s, min=-0.5, max=0.5)
  Gammaz = matrix(runif(s, min=-0.5, max=0.5),nrow=s, ncol=p)
  M <- X%*%t(Gamma) + S %*% Gammaz + gen_error(n,p,rho = 0.5) ## n by p
  Y <- as.vector(M %*% alpham + X%*% alphax + S%*% alphaz + rnorm(n,0,sd=0.5)) ## n by 1
  
  ## PEHDGMM
  output_PEHDMM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), scale(S),
                                       scale=FALSE, lambda_grid = lamb_grid, 
                                       lambda_grid_reduced =lamb_grid0)
  pval_PEHDMM <- output_PEHDMM$value[output_PEHDMM$term == "pval_pe"]
  
  
  ## Guo et al. (2022, JoE; 2022, JASA)
  output_HDMM <- HDGMM_linear(scale(X), Y-mean(Y), scale(M), scale(S), 
                              scale=FALSE,  lamb_grid=lamb_grid, lamb_grid0=lamb_grid0)
  pval_HDMM <- output_HDMM$pval_HDGMM
  
  ## Zhou et al. (2020)
  if(excludeComparison)
    {
    output_zhou <- NULL
    pval_zhou <- NA
    } else{
    output_zhou <- Zhou2020(scale(X), Y-mean(Y), scale(M), scale(S))
    pval_zhou <- output_zhou$ZWZ_Sn_pvalue
    }
  
  }

  ## Combine results
  pvals <- c(pval_PEHDMM, pval_HDMM, pval_zhou)
  ans <- rbind(ans, pvals)
  }

  res <- apply(ans<0.05,2,function(x) mean(x, na.rm=T))
  output <- c(n=n, p=p, c1 = c1, c2 = c2, res)
  output_summary_tab <- rbind(output_summary_tab, output)
  ans_tab <- cbind(rep(n,nrow(ans)), rep(p,nrow(ans)), rep(c1,nrow(ans)), rep(c2, nrow(ans)), ans)
  output_pvals <- rbind(output_pvals, ans_tab)
}


colnames(output_summary_tab) <- c("n", "p", "c1", "c2", "pval_PEHDMM", "pval_HDMM", "pval_zhou")

colnames(output_pvals) <- c("n", "p", "c1", "c2", "pval_PEHDMM", "pval_HDMM", "pval_zhou")


write.csv(output_summary_tab, 
          file = paste0("summary_LM_Basic_pvals_setting", setting_id, "_parId", parid, ".csv"), 
          row.names = F)
write.csv(output_pvals, 
          file =  paste0("simu_LM_Basic_pvals_setting", setting_id, "_parId", parid, ".csv"),
          row.names = F)
