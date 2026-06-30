##########################################################
##  Simulation studies of PE-HDGMM for Poisson outcomes ##
##########################################################

parid <- 21 # taking values in c(1:21)
setting_id <- 2 # taking values in {1,2}
includeCovariate <- FALSE # taking values in {TRUE, FALSE}

library(POEM)
source("code-utils-comparison-methods-poisson.R")

n <- 300
q <- 1
p <- 500
c2 <- 0.4
repnum <- 1000

Y_family <- "poisson"
penalty_type <- "SCAD"

## Signal pattern setup 
{
  if(setting_id==1) {
    alpham = matrix(c(0.9,0.8,0,0,0.7,rep(0,p-5)),ncol=1) # setting 1
  }
  if(setting_id==2) {
    alpham = matrix(c(0,0,0,0.8,-0.7,rep(0,p-5)),ncol=1) # setting 2
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

gamma = matrix(NA,nrow = p,ncol = 1)
gamma[1:5] = c(0.2,0.4,0.6,0.8,1)/2
gamma[6:p] = rnorm((p-5),mean=0,sd=0.01)
alphax = c2

phi0 = 1
ngrid = 15
lamb_grid = seq(0.5, 2.5,length.out = ngrid)

output_summary_tab <- c() 
output_pvals <- c() 

c1_list <- seq(-1, 1, by=0.1)
for (c1 in c1_list[parid])
{
  print(paste0("c1:", c1))
  beta0 = t(alpham)%*%gamma
  gamma_true=c1*gamma[1:5]
  beta_true = c1 * beta0
  Gamma = c1*t(gamma)

  ans <- c() 
  
  for(ite in 1:repnum)
  {
    print(ite)
    
    tryCatch({
   
    ## NOTE: 
    ## X: a matrix of covariates, n by q
    ## Y: a vector of responses, n by 1
    ## M: a matrix of mediators, n by p
    ## S: a matrix of confounders, n by s 
  
    if(!includeCovariate)
    {
      ## no confounders 
      s <- 0 
      X = as.matrix(rnorm(n,mean=0,sd=1),ncol=q)
      M = X%*%Gamma +gen_error(n,p,rho=0.5)
      Z = M%*%alpham + X%*%alphax
      if(max(Z) > 5){
        if(max(Z)==min(Z)){
          Z_scaled <- rep(0, length(Z))
        }else{
          Z_scaled <- -5 + (Z - min(Z)) * (10) / (max(Z) - min(Z))
        }
      }else{
        Z_scaled <- Z
      }
      lambda = exp(Z_scaled)
      Y = rpois(n, lambda)
      S = NULL
    } else{
        ## with confounders
        s <- 5
        X = as.matrix(rnorm(n,mean=0,sd=1),ncol=q)
        S = matrix(rnorm(n*s,mean=0,sd=1),ncol=s)
        alphaz = runif(s, min=-0.5, max=0.5)
        Gammaz = matrix(runif(s, min=-0.5, max=0.5),nrow=s, ncol=p)
        M = X%*%Gamma + S %*% Gammaz +gen_error(n,p,rho=0.5)
        Z = M%*%alpham + X%*%alphax+ S%*% alphaz
        if(max(Z) > 5){
          if(max(Z)==min(Z)){
            Z_scaled <- rep(0, length(Z))
          }else{
            Z_scaled <- -5 + (Z - min(Z)) * (10) / (max(Z) - min(Z))
          }
        }else{
          Z_scaled <- Z
        }
        
        lambda = exp(Z_scaled)
        Y = rpois(n, lambda)
    }
      
  ## PEHDGMM
  output_PEHDGMM <- pe_mediation_poisson(scale(X), Y, scale(M), S,  
                                        scale=FALSE, lambda_grid = lamb_grid)
  pval_PEHDGMM <- output_PEHDGMM$value[output_PEHDGMM$term == "pval_pe"]
  
  ## HDGMM (Guo et al., 2023)
  output_HDGMM <- HDGMM_poisson(scale(X), Y, scale(M), S, 
                                scale=FALSE, lamb_grid=lamb_grid)
  pval_HDGMM <- output_HDGMM$pval_HDGMM
  

  ## Djordjilovic et al. (2019)'s method
  output_Vera <- Vera2019_Poisson(scale(X), Y, scale(M), S=S)
  pval_Vera <- output_Vera$pval_Vera
  
  ## Combine results
  pvals <- c(pval_PEHDGMM, pval_HDGMM, pval_Vera)
  ans <- rbind(ans, pvals)
  },
  error = function(e) {
    print(paste("An error occurred at iteration", ite, "with message:", e$message))
  })
  }

  res <- apply(ans<0.05,2,mean)
  output <- c(c1 = c1, c2 = c2, res)
  output_summary_tab <- rbind(output_summary_tab, output)
  ans_tab <- cbind(rep(c1,nrow(ans)), rep(c2, nrow(ans)), ans)
  output_pvals <- rbind(output_pvals, ans_tab)

}

colnames(output_summary_tab) <- c("c1", "c2", "pval_PEHDGMM", "pval_HDGMM", "pval_Vera")

colnames(output_pvals) <- c("c1", "c2", "pval_PEHDGMM", "pval_HDGMM", "pval_Vera")

write.csv(output_summary_tab,
          file = paste0("summary_Poisson_pvals_setting", setting_id, "_parId", parid, ".csv"),
          row.names = F)
write.csv(output_pvals, 
          file =  paste0("simu_Poisson_pvals_setting", setting_id, "_parId", parid, ".csv"), 
          row.names = F)

