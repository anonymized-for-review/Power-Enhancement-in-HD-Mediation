###############################################################################
######      Real Data Analysis: PE-HDMM vs HDMM                          ######
######    Using Infant Mortality Rate (IMR) as the outcome variable      ######  
###############################################################################

###############################################################################
## Load Data ## 
rm(list=ls())
library(POEM)
load("data/data_GDPR_GHED_IMR.RData")
dim(data_gdpr) # 91 by 23
dim(data_imr)  # 91 by 23
dim(data_ghed) # 2002 by 62
## columns: country, code, region, income, year, 57 indicators 
## rows: 91 countries * 22 years 

## sanity check 
sum(data_gdpr$Country.Code != sort(data_gdpr$Country.Code))
sum(data_imr$Country.Code != sort(data_imr$Country.Code))
sum(data_ghed$code != sort(data_ghed$code))
###############################################################################

###############################################################################
## ----------------------------------------- ##
## Model Fitting 1: A Global Mediation Model ## 
## ----------------------------------------- ## 

## transforming the format of data_gdpr 
code_column <- c()
gdp_column <- c()
year_column <- c()
for(i in 1:nrow(data_gdpr))
{
  code_temp <- rep(data_gdpr[i,]$Country.Code,22)
  gdp_temp <- unlist(data_gdpr[i,2:ncol(data_gdpr)])
  code_column <- c(code_column, code_temp)
  gdp_column <- c(gdp_column, gdp_temp)
  year_column <- c(year_column, 2000:2021)
}
data_gdpr_mat <- data.frame(code = code_column,
                            year = year_column, 
                            gdp = unname(gdp_column))
dim(data_gdpr_mat) # 2002 by 3

## transforming the format of data_imr
code_column <- c()
imr_column <- c()
year_column <- c()
for(i in 1:nrow(data_imr))
{
  code_temp <- rep(data_imr[i,]$Country.Code,22)
  imr_temp <- unlist(data_imr[i,2:ncol(data_imr)])
  code_column <- c(code_column, code_temp)
  imr_column <- c(imr_column, imr_temp)
  year_column <- c(year_column, 2000:2021)
}
data_imr_mat <- data.frame(code = code_column,
                           year = year_column, 
                           imr = unname(imr_column))
dim(data_imr_mat)  # 2002 by 3

## preprocessing data_ghed
data_ghed$country <- as.factor(data_ghed$country)
data_ghed$code <- as.factor(data_ghed$code)
data_ghed$region <- as.factor(data_ghed$region)
data_ghed$income <- as.factor(data_ghed$income)
data_gdpr_mat$code <- as.factor(data_gdpr_mat$code)
data_imr_mat$code <- as.factor(data_imr_mat$code)

## model fitting 
X <- matrix(data_gdpr_mat$gdp,ncol=1)
Y <- data_imr_mat$imr
M <- data_ghed[,6:ncol(data_ghed)]
S <- data_ghed[,c("region", "income", "year")]
## S contains categorical variables --> one-hot encoding
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]

ngrid = 100 ## set up HBIC for different 
lamb_grid = seq(0.1,5,length.out = ngrid) 
lamb_grid0 = seq(0.1,5,length.out = ngrid) 

# HDMM & PE-HDMM
output_PE <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                 scale=FALSE, lambda_grid = lamb_grid, 
                                 lambda_grid_reduced =lamb_grid0)
output_PE$value[output_PE$term == "pval_hdmm"]  # 0.0409  
output_PE$value[output_PE$term == "pval_pe"]    # 0  
colnames(M)[attr(output_PE, "active_mediators")] # "gge_gdp

###############################################################################

###############################################################################
## ------------------------------------------ ##
## Model Fitting 2: Regional Mediation Models ## 
## ------------------------------------------ ##
unique(data_ghed$region)
country_list_AFR  <- unique(data_ghed[data_ghed$region == 'AFR', 'code'])
country_list_AMR  <- unique(data_ghed[data_ghed$region == 'AMR', 'code'])
country_list_EMR  <- unique(data_ghed[data_ghed$region == 'EMR', 'code'])
country_list_EUR  <- unique(data_ghed[data_ghed$region == 'EUR', 'code'])
country_list_SEAR <- unique(data_ghed[data_ghed$region == 'SEAR', 'code'])
country_list_WPR  <- unique(data_ghed[data_ghed$region == 'WPR', 'code'])
length(country_list_AFR)  # 32 
length(country_list_AMR)  # 28
length(country_list_EMR)  #  8
length(country_list_EUR)  #  9
length(country_list_SEAR) #  5
length(country_list_WPR)  #  9

data_gdpr_mat_AFR  <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_AFR, ]
data_gdpr_mat_AMR  <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_AMR, ]
data_gdpr_mat_EMR  <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_EMR, ]
data_gdpr_mat_EUR  <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_EUR, ]
data_gdpr_mat_SEAR <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_SEAR, ]
data_gdpr_mat_WPR  <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_WPR, ]

data_imr_mat_AFR  <- data_imr_mat[data_imr_mat$code %in% country_list_AFR,]
data_imr_mat_AMR  <- data_imr_mat[data_imr_mat$code %in% country_list_AMR,]
data_imr_mat_EMR  <- data_imr_mat[data_imr_mat$code %in% country_list_EMR,]
data_imr_mat_EUR  <- data_imr_mat[data_imr_mat$code %in% country_list_EUR,]
data_imr_mat_SEAR <- data_imr_mat[data_imr_mat$code %in% country_list_SEAR,]
data_imr_mat_WPR  <- data_imr_mat[data_imr_mat$code %in% country_list_WPR,]

data_ghed_AFR  <- data_ghed[data_ghed$code %in% country_list_AFR,]
data_ghed_AMR  <- data_ghed[data_ghed$code %in% country_list_AMR,]
data_ghed_EMR  <- data_ghed[data_ghed$code %in% country_list_EMR,]
data_ghed_EUR  <- data_ghed[data_ghed$code %in% country_list_EUR,]
data_ghed_SEAR <- data_ghed[data_ghed$code %in% country_list_SEAR,]
data_ghed_WPR  <- data_ghed[data_ghed$code %in% country_list_WPR,]

## region == 'AFR' (32 countries)
X <- matrix(data_gdpr_mat_AFR$gdp,ncol=1); dim(X) # 704
Y <- data_imr_mat_AFR$imr; length(Y)
M <- data_ghed_AFR[,6:ncol(data_ghed_AFR)]; dim(M)
S <- data_ghed_AFR[,c("income", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_region_AFR_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                              scale = FALSE,
                                              lambda_grid = lamb_grid,
                                              lambda_grid_reduced = lamb_grid0)
output_region_AFR_POEM$value[output_region_AFR_POEM$term == "pval_hdmm"] # 0.4672 
output_region_AFR_POEM$value[output_region_AFR_POEM$term == "pval_pe"] # 0.4672 
colnames(M)[attr(output_region_AFR_POEM, "active_mediators")] # null  

## region == 'AMR' (28 countries)
X <- matrix(data_gdpr_mat_AMR$gdp,ncol=1); dim(X) # 616
Y <- data_imr_mat_AMR$imr; length(Y)
M <- data_ghed_AMR[,6:ncol(data_ghed_AMR)]; dim(M)
S <- data_ghed_AMR[,c("income", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_region_AMR_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                              scale = FALSE,
                                              lambda_grid = lamb_grid,
                                              lambda_grid_reduced = lamb_grid0)
output_region_AMR_POEM$value[output_region_AMR_POEM$term == "pval_hdmm"] # 0.4129
output_region_AMR_POEM$value[output_region_AMR_POEM$term == "pval_pe"] # 0.4129 
colnames(M)[attr(output_region_AMR_POEM, "active_mediators")] # null 

## region == 'EMR' (8 countries)
X <- matrix(data_gdpr_mat_EMR$gdp,ncol=1); dim(X) # 176
Y <- data_imr_mat_EMR$imr; length(Y)
M <- data_ghed_EMR[,6:ncol(data_ghed_EMR)]; dim(M)
S <- data_ghed_EMR[,c("income", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_region_EMR_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                              scale = FALSE,
                                              lambda_grid = lamb_grid,
                                              lambda_grid_reduced = lamb_grid0)
output_region_EMR_POEM$value[output_region_EMR_POEM$term == "pval_hdmm"] # 0.0304
output_region_EMR_POEM$value[output_region_EMR_POEM$term == "pval_pe"]   # 0.0304 
colnames(M)[attr(output_region_EMR_POEM, "active_mediators")] # null 
 
## region == 'EUR' (9 countries)
X <- matrix(data_gdpr_mat_EUR$gdp,ncol=1); dim(X) # 198
Y <- data_imr_mat_EUR$imr; length(Y)
M <- data_ghed_EUR[,6:ncol(data_ghed_EUR)]; dim(M)
S <- data_ghed_EUR[,c("income", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_region_EUR_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                              scale = FALSE,
                                              lambda_grid = lamb_grid,
                                              lambda_grid_reduced = lamb_grid0)
output_region_EUR_POEM$value[output_region_EUR_POEM$term == "pval_hdmm"]  # 0.0035 
output_region_EUR_POEM$value[output_region_EUR_POEM$term == "pval_pe"] # 0
colnames(M)[attr(output_region_EUR_POEM, "active_mediators")] # "gge_gdp" 

## region == 'SEAR' (5 countries)
X <- matrix(data_gdpr_mat_SEAR$gdp,ncol=1); dim(X) # 110
Y <- data_imr_mat_SEAR$imr; length(Y)
M <- data_ghed_SEAR[,6:ncol(data_ghed_SEAR)]; dim(M)
S <- data_ghed_SEAR[,c("income", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_region_SEAR_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                               scale = FALSE,
                                               lambda_grid = lamb_grid,
                                               lambda_grid_reduced = lamb_grid0)
output_region_SEAR_POEM$value[output_region_SEAR_POEM$term == "pval_hdmm"] # 0.3087 
output_region_SEAR_POEM$value[output_region_SEAR_POEM$term == "pval_pe"]   # 0
colnames(M)[attr(output_region_SEAR_POEM, "active_mediators")] # "ext_usd2021_pc"

## region == 'WPR' (9 countries)
X <- matrix(data_gdpr_mat_WPR$gdp,ncol=1); dim(X) # 198
Y <- data_imr_mat_WPR$imr; length(Y)
M <- data_ghed_WPR[,6:ncol(data_ghed_WPR)]; dim(M)
S <- data_ghed_WPR[,c("income", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_region_WPR_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                              scale = FALSE,
                                              lambda_grid = lamb_grid,
                                              lambda_grid_reduced = lamb_grid0)
output_region_WPR_POEM$value[output_region_WPR_POEM$term == "pval_hdmm"] # 0.8721  
output_region_WPR_POEM$value[output_region_WPR_POEM$term == "pval_pe"]   # 0  
colnames(M)[attr(output_region_WPR_POEM, "active_mediators")]  # "chi_che"  "pvtd_gdp" 
###############################################################################

###############################################################################
## ------------------------------------------------ ##
## Model Fitting 3: Group Data by Income Categories ## 
## ------------------------------------------------ ##
unique(data_ghed$income)
country_list_High  <- unique(data_ghed[data_ghed$income == 'High', 'code'])
country_list_UM  <- unique(data_ghed[data_ghed$income == 'Upper-middle', 'code']) 
country_list_LM  <- unique(data_ghed[data_ghed$income == 'Lower-middle', 'code'])
country_list_Low  <- unique(data_ghed[data_ghed$income == 'Low', 'code'])
length(country_list_High) # 22
length(country_list_UM)   # 26
length(country_list_LM)   # 27
length(country_list_Low)  # 16 

data_gdpr_mat_High <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_High, ]
data_gdpr_mat_UM   <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_UM, ]
data_gdpr_mat_LM   <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_LM, ]
data_gdpr_mat_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% country_list_Low, ]

data_imr_mat_High <-  data_imr_mat[data_imr_mat$code %in% country_list_High, ]
data_imr_mat_UM   <-  data_imr_mat[data_imr_mat$code %in% country_list_UM, ]
data_imr_mat_LM   <-  data_imr_mat[data_imr_mat$code %in% country_list_LM, ]
data_imr_mat_Low  <-  data_imr_mat[data_imr_mat$code %in% country_list_Low, ]

data_ghed_High <- data_ghed[data_ghed$code %in% country_list_High,]
data_ghed_UM   <- data_ghed[data_ghed$code %in% country_list_UM,]
data_ghed_LM   <- data_ghed[data_ghed$code %in% country_list_LM,]
data_ghed_Low  <- data_ghed[data_ghed$code %in% country_list_Low,]

## income == 'Low' (16 countries)
X <- matrix(data_gdpr_mat_Low$gdp,ncol=1); dim(X) # 352
Y <- data_imr_mat_Low$imr; length(Y)
M <- data_ghed_Low[,6:ncol(data_ghed_Low)]; dim(M)
#  S <- data_ghed_Low[,c("region", "year")]; dim(S)
#  S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
## Error:  contrasts can be applied only to factors with 2 or more levels
## Soln: remove "region" from S
S <- matrix(data_ghed_Low[,c("year")], ncol=1); dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_income_Low_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                              scale = FALSE,
                                              lambda_grid = lamb_grid,
                                              lambda_grid_reduced = lamb_grid0)
output_income_Low_POEM$value[output_income_Low_POEM$term == "pval_hdmm"] # 0.3157
output_income_Low_POEM$value[output_income_Low_POEM$term == "pval_pe"]   # 0 
colnames(M)[attr(output_income_Low_POEM, "active_mediators")] # "pvtd_usd2021"

## income == 'Lower-middle' (27 countries)
X <- matrix(data_gdpr_mat_LM$gdp,ncol=1); dim(X) # 594
Y <- data_imr_mat_LM$imr; length(Y)
M <- data_ghed_LM[,6:ncol(data_ghed_LM)]; dim(M)
S <- data_ghed_LM[,c("region", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_income_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                             scale = FALSE,
                                             lambda_grid = lamb_grid,
                                             lambda_grid_reduced = lamb_grid0)
output_income_LM_POEM$value[output_income_LM_POEM$term == "pval_hdmm"] # 0.9546
output_income_LM_POEM$value[output_income_LM_POEM$term == "pval_pe"]   # 0.9546 
colnames(M)[attr(output_income_LM_POEM, "active_mediators")] # null 

## income == 'Upper-middle' (26 countries)
X <- matrix(data_gdpr_mat_UM$gdp,ncol=1); dim(X) # 572
Y <- data_imr_mat_UM$imr; length(Y)
M <- data_ghed_UM[,6:ncol(data_ghed_UM)]; dim(M)
S <- data_ghed_UM[,c("region", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_income_UM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                             scale = FALSE,
                                             lambda_grid = lamb_grid,
                                             lambda_grid_reduced = lamb_grid0)
output_income_UM_POEM$value[output_income_UM_POEM$term == "pval_hdmm"] # 0.4945 
output_income_UM_POEM$value[output_income_UM_POEM$term == "pval_pe"]   # 0.4945 
colnames(M)[attr(output_income_UM_POEM, "active_mediators")] # null

## income == 'High' (22 countries)
X <- matrix(data_gdpr_mat_High$gdp,ncol=1); dim(X) # 484
Y <- data_imr_mat_High$imr; length(Y)
M <- data_ghed_High[,6:ncol(data_ghed_High)]; dim(M)
S <- data_ghed_High[,c("region", "year")]; dim(S)
S <- model.matrix(lm(Y~., data=data.frame(Y, S)))[,-1]; dim(S)
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_income_High_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                               scale = FALSE,
                                               lambda_grid = lamb_grid,
                                               lambda_grid_reduced = lamb_grid0)
output_income_High_POEM$value[output_income_High_POEM$term == "pval_hdmm"] # 0.0238
output_income_High_POEM$value[output_income_High_POEM$term == "pval_pe"] # 0
colnames(M)[attr(output_income_High_POEM, "active_mediators")] # "oops_che" "shi_che"
###############################################################################

###############################################################################
## --------------------------------------------- ##
## Model Fitting 4: Group Data by Region*Income  ## 
## --------------------------------------------- ##

length(intersect(country_list_AFR, country_list_High)) #  0
length(intersect(country_list_AFR, country_list_UM))   #  4
length(intersect(country_list_AFR, country_list_LM))   # 12
length(intersect(country_list_AFR, country_list_Low))  # 16 

length(intersect(country_list_AMR, country_list_High)) #  9
length(intersect(country_list_AMR, country_list_UM))   # 15
length(intersect(country_list_AMR, country_list_LM))   #  4
length(intersect(country_list_AMR, country_list_Low))  #  0 

length(intersect(country_list_EMR, country_list_High)) #  2
length(intersect(country_list_EMR, country_list_UM))   #  1
length(intersect(country_list_EMR, country_list_LM))   #  5
length(intersect(country_list_EMR, country_list_Low))  #  0 

length(intersect(country_list_EUR, country_list_High)) #  8
length(intersect(country_list_EUR, country_list_UM))   #  0
length(intersect(country_list_EUR, country_list_LM))   #  1
length(intersect(country_list_EUR, country_list_Low))  #  0

length(intersect(country_list_SEAR, country_list_High)) # 0
length(intersect(country_list_SEAR, country_list_UM))   # 2
length(intersect(country_list_SEAR, country_list_LM))   # 3
length(intersect(country_list_SEAR, country_list_Low))  # 0

length(intersect(country_list_WPR, country_list_High)) #  3
length(intersect(country_list_WPR, country_list_UM))   #  4
length(intersect(country_list_WPR, country_list_LM))   #  2
length(intersect(country_list_WPR, country_list_Low))  #  0

data_gdpr_mat_AFR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AFR, country_list_High), ]
data_gdpr_mat_AFR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AFR, country_list_UM), ]
data_gdpr_mat_AFR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AFR, country_list_LM), ]
data_gdpr_mat_AFR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AFR, country_list_Low), ]

data_gdpr_mat_AMR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AMR, country_list_High), ]
data_gdpr_mat_AMR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AMR, country_list_UM), ]
data_gdpr_mat_AMR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AMR, country_list_LM), ]
data_gdpr_mat_AMR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_AMR, country_list_Low), ]

data_gdpr_mat_EMR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EMR, country_list_High), ]
data_gdpr_mat_EMR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EMR, country_list_UM), ]
data_gdpr_mat_EMR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EMR, country_list_LM), ]
data_gdpr_mat_EMR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EMR, country_list_Low), ]

data_gdpr_mat_EUR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_High), ]
data_gdpr_mat_EUR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_UM), ]
data_gdpr_mat_EUR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_LM), ]
data_gdpr_mat_EUR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_Low), ]

data_gdpr_mat_EUR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_High), ]
data_gdpr_mat_EUR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_UM), ]
data_gdpr_mat_EUR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_LM), ]
data_gdpr_mat_EUR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_EUR, country_list_Low), ]

data_gdpr_mat_SEAR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_SEAR, country_list_High), ]
data_gdpr_mat_SEAR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_SEAR, country_list_UM), ]
data_gdpr_mat_SEAR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_SEAR, country_list_LM), ]
data_gdpr_mat_SEAR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_SEAR, country_list_Low), ]

data_gdpr_mat_WPR_High  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_WPR, country_list_High), ]
data_gdpr_mat_WPR_UM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_WPR, country_list_UM), ]
data_gdpr_mat_WPR_LM  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_WPR, country_list_LM), ]
data_gdpr_mat_WPR_Low  <- data_gdpr_mat[data_gdpr_mat$code %in% intersect(country_list_WPR, country_list_Low), ]

data_imr_mat_AFR_High  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AFR, country_list_High), ]
data_imr_mat_AFR_UM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AFR, country_list_UM), ]
data_imr_mat_AFR_LM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AFR, country_list_LM), ]
data_imr_mat_AFR_Low  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AFR, country_list_Low), ]

data_imr_mat_AMR_High  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AMR, country_list_High), ]
data_imr_mat_AMR_UM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AMR, country_list_UM), ]
data_imr_mat_AMR_LM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AMR, country_list_LM), ]
data_imr_mat_AMR_Low  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_AMR, country_list_Low), ]

data_imr_mat_EMR_High  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EMR, country_list_High), ]
data_imr_mat_EMR_UM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EMR, country_list_UM), ]
data_imr_mat_EMR_LM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EMR, country_list_LM), ]
data_imr_mat_EMR_Low  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EMR, country_list_Low), ]

data_imr_mat_EUR_High  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EUR, country_list_High), ]
data_imr_mat_EUR_UM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EUR, country_list_UM), ]
data_imr_mat_EUR_LM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EUR, country_list_LM), ]
data_imr_mat_EUR_Low  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_EUR, country_list_Low), ]

data_imr_mat_SEAR_High  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_SEAR, country_list_High), ]
data_imr_mat_SEAR_UM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_SEAR, country_list_UM), ]
data_imr_mat_SEAR_LM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_SEAR, country_list_LM), ]
data_imr_mat_SEAR_Low  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_SEAR, country_list_Low), ]

data_imr_mat_WPR_High  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_WPR, country_list_High), ]
data_imr_mat_WPR_UM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_WPR, country_list_UM), ]
data_imr_mat_WPR_LM  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_WPR, country_list_LM), ]
data_imr_mat_WPR_Low  <- data_imr_mat[data_imr_mat$code %in% intersect(country_list_WPR, country_list_Low), ]

data_ghed_AFR_High  <- data_ghed[data_ghed$code %in% intersect(country_list_AFR, country_list_High), ]
data_ghed_AFR_UM  <- data_ghed[data_ghed$code %in% intersect(country_list_AFR, country_list_UM), ]
data_ghed_AFR_LM  <- data_ghed[data_ghed$code %in% intersect(country_list_AFR, country_list_LM), ]
data_ghed_AFR_Low  <- data_ghed[data_ghed$code %in% intersect(country_list_AFR, country_list_Low), ]

data_ghed_AMR_High  <- data_ghed[data_ghed$code %in% intersect(country_list_AMR, country_list_High), ]
data_ghed_AMR_UM  <- data_ghed[data_ghed$code %in% intersect(country_list_AMR, country_list_UM), ]
data_ghed_AMR_LM  <- data_ghed[data_ghed$code %in% intersect(country_list_AMR, country_list_LM), ]
data_ghed_AMR_Low  <- data_ghed[data_ghed$code %in% intersect(country_list_AMR, country_list_Low), ]

data_ghed_EMR_High  <- data_ghed[data_ghed$code %in% intersect(country_list_EMR, country_list_High), ]
data_ghed_EMR_UM  <- data_ghed[data_ghed$code %in% intersect(country_list_EMR, country_list_UM), ]
data_ghed_EMR_LM  <- data_ghed[data_ghed$code %in% intersect(country_list_EMR, country_list_LM), ]
data_ghed_EMR_Low  <- data_ghed[data_ghed$code %in% intersect(country_list_EMR, country_list_Low), ]

data_ghed_EUR_High  <- data_ghed[data_ghed$code %in% intersect(country_list_EUR, country_list_High), ]
data_ghed_EUR_UM  <- data_ghed[data_ghed$code %in% intersect(country_list_EUR, country_list_UM), ]
data_ghed_EUR_LM  <- data_ghed[data_ghed$code %in% intersect(country_list_EUR, country_list_LM), ]
data_ghed_EUR_Low  <- data_ghed[data_ghed$code %in% intersect(country_list_EUR, country_list_Low), ]

data_ghed_SEAR_High  <- data_ghed[data_ghed$code %in% intersect(country_list_SEAR, country_list_High), ]
data_ghed_SEAR_UM  <- data_ghed[data_ghed$code %in% intersect(country_list_SEAR, country_list_UM), ]
data_ghed_SEAR_LM  <- data_ghed[data_ghed$code %in% intersect(country_list_SEAR, country_list_LM), ]
data_ghed_SEAR_Low  <- data_ghed[data_ghed$code %in% intersect(country_list_SEAR, country_list_Low), ]

data_ghed_WPR_High  <- data_ghed[data_ghed$code %in% intersect(country_list_WPR, country_list_High), ]
data_ghed_WPR_UM  <- data_ghed[data_ghed$code %in% intersect(country_list_WPR, country_list_UM), ]
data_ghed_WPR_LM  <- data_ghed[data_ghed$code %in% intersect(country_list_WPR, country_list_LM), ]
data_ghed_WPR_Low  <- data_ghed[data_ghed$code %in% intersect(country_list_WPR, country_list_Low), ]

## region=='AFR' & income=='Low'
X <- matrix(data_gdpr_mat_AFR_Low$gdp,ncol=1); dim(X) # 352
Y <- data_imr_mat_AFR_Low$imr; length(Y)
M <- data_ghed_AFR_Low[,6:ncol(data_ghed_AFR_Low)]; dim(M)
S <- matrix(data_ghed_AFR_Low[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # TRUE
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_AFR_Low_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                           scale = FALSE,
                                           lambda_grid = lamb_grid,
                                           lambda_grid_reduced = lamb_grid0)
output_AFR_Low_POEM$value[output_AFR_Low_POEM$term == "pval_hdmm"] # 0.3157 
output_AFR_Low_POEM$value[output_AFR_Low_POEM$term == "pval_pe"]   # 0
colnames(M)[attr(output_AFR_Low_POEM, "active_mediators")] # "pvtd_usd2021"

## region=='AFR' & income=='Lower-middle'
X <- matrix(data_gdpr_mat_AFR_LM$gdp,ncol=1); dim(X) # 264
Y <- data_imr_mat_AFR_LM$imr; length(Y)
M <- data_ghed_AFR_LM[,6:ncol(data_ghed_AFR_LM)]; dim(M)
S <- matrix(data_ghed_AFR_LM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 264 by 56
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_AFR_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_AFR_LM_POEM$value[output_AFR_LM_POEM$term == "pval_hdmm"] # 0.5947
output_AFR_LM_POEM$value[output_AFR_LM_POEM$term == "pval_pe"]   # 0.5947
colnames(M)[attr(output_AFR_LM_POEM, "active_mediators")] # null 

## region=='AFR' & income=='Upper-middle'
X <- matrix(data_gdpr_mat_AFR_UM$gdp,ncol=1); dim(X) # 88
Y <- data_imr_mat_AFR_UM$imr; length(Y) 
M <- data_ghed_AFR_UM[,6:ncol(data_ghed_AFR_UM)]; dim(M)
S <- matrix(data_ghed_AFR_UM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 88 by 56
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_AFR_UM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_AFR_UM_POEM$value[output_AFR_UM_POEM$term == "pval_hdmm"] # 0.0080 
output_AFR_UM_POEM$value[output_AFR_UM_POEM$term == "pval_pe"]   # 0
colnames(M)[attr(output_AFR_UM_POEM, "active_mediators")] # "pvtd_ncu2021_pc" 

## region=='AFR' & income=='High' -- no samples
X <- matrix(data_gdpr_mat_AFR_High$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_AFR_High$imr; length(Y) # 0

## region=='AMR' & income=='Low' -- no samples
X <- matrix(data_gdpr_mat_AMR_Low$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_AMR_Low$imr; length(Y)  

## region=='AMR' & income=='Lower-middle'
X <- matrix(data_gdpr_mat_AMR_LM$gdp,ncol=1); dim(X)  # 88
Y <- data_imr_mat_AMR_LM$imr; length(Y) 
M <- data_ghed_AMR_LM[,6:ncol(data_ghed_AMR_LM)]; dim(M)
S <- matrix(data_ghed_AMR_LM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # TRUE
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_AMR_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_AMR_LM_POEM$value[output_AMR_LM_POEM$term == "pval_hdmm"] # 0.1295
output_AMR_LM_POEM$value[output_AMR_LM_POEM$term == "pval_pe"]   # 0.1295 
colnames(M)[attr(output_AMR_LM_POEM, "active_mediators")] # null

## region=='AMR' & income=='Upper-middle'
X <- matrix(data_gdpr_mat_AMR_UM$gdp,ncol=1); dim(X) # 330
Y <- data_imr_mat_AMR_UM$imr; length(Y) 
M <- data_ghed_AMR_UM[,6:ncol(data_ghed_AMR_UM)]; dim(M)
S <- matrix(data_ghed_AMR_UM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # TRUE
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_AMR_UM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_AMR_UM_POEM$value[output_AMR_UM_POEM$term == "pval_hdmm"] # 0.0426 
output_AMR_UM_POEM$value[output_AMR_UM_POEM$term == "pval_pe"]   # 0 
colnames(M)[attr(output_AMR_UM_POEM, "active_mediators")] # "gge_gdp"
 
## region=='AMR' & income=='High'
X <- matrix(data_gdpr_mat_AMR_High$gdp,ncol=1); dim(X) # 198
Y <- data_imr_mat_AMR_High$imr; length(Y) 
M <- data_ghed_AMR_High[,6:ncol(data_ghed_AMR_High)]; dim(M)
S <- matrix(data_ghed_AMR_High[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # TRUE
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_AMR_High_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                            scale = FALSE,
                                            lambda_grid = lamb_grid,
                                            lambda_grid_reduced = lamb_grid0)
output_AMR_High_POEM$value[output_AMR_High_POEM$term == "pval_hdmm"] # 0.0606
output_AMR_High_POEM$value[output_AMR_High_POEM$term == "pval_pe"]   # 0.0606 
colnames(M)[attr(output_AMR_High_POEM, "active_mediators")]  # null 
 
## region=='EMR' & income=='Low' -- no sample
unique(data_ghed_EMR_Low$country) 
X <- matrix(data_gdpr_mat_EMR_Low$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_EMR_Low$imr; length(Y) # 0

## region=='EMR' & income=='Lower-middle' 
unique(data_ghed_EMR_LM$country) 
X <- matrix(data_gdpr_mat_EMR_LM$gdp,ncol=1); dim(X) # 110
Y <- data_imr_mat_EMR_LM$imr; length(Y) 
M <- data_ghed_EMR_LM[,6:ncol(data_ghed_EMR_LM)]; dim(M)
S <- matrix(data_ghed_EMR_LM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # TRUE
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_EMR_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_EMR_LM_POEM$value[output_EMR_LM_POEM$term == "pval_hdmm"] # 0.0884
output_EMR_LM_POEM$value[output_EMR_LM_POEM$term == "pval_pe"]   # 0.0884 
colnames(M)[attr(output_EMR_LM_POEM, "active_mediators")]  # null

## region=='EMR' & income=='Upper-middle'
unique(data_ghed_EMR_UM$country)
X <- matrix(data_gdpr_mat_EMR_UM$gdp,ncol=1); dim(X) # 22
Y <- data_imr_mat_EMR_UM$imr; length(Y) # 22 
M <- data_ghed_EMR_UM[,6:ncol(data_ghed_EMR_UM)]; dim(M)
S <- matrix(data_ghed_EMR_UM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_EMR_UM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_EMR_UM_POEM$value[output_EMR_UM_POEM$term == "pval_hdmm"] # 0.0015
output_EMR_UM_POEM$value[output_EMR_UM_POEM$term == "pval_pe"]   # 0.0015
colnames(M)[attr(output_EMR_UM_POEM, "active_mediators")]  # null

## region=='EMR' & income=='High'
unique(data_ghed_EMR_High$country) 
X <- matrix(data_gdpr_mat_EMR_High$gdp,ncol=1); dim(X) # 44
Y <- data_imr_mat_EMR_High$imr; length(Y) # 44
M <- data_ghed_EMR_High[,6:ncol(data_ghed_EMR_High)]; dim(M)
S <- matrix(data_ghed_EMR_High[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0 ## FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 44 by 44
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_EMR_High_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                            scale = FALSE,
                                            lambda_grid = lamb_grid,
                                            lambda_grid_reduced = lamb_grid0)
output_EMR_High_POEM$value[output_EMR_High_POEM$term == "pval_hdmm"] # 0.6096 
output_EMR_High_POEM$value[output_EMR_High_POEM$term == "pval_pe"]   # 0.6096 
colnames(M)[attr(output_EMR_High_POEM, "active_mediators")] # null

## region=='EUR' & income=='Low' -- no samples
X <- matrix(data_gdpr_mat_EUR_Low$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_EUR_Low$imr; length(Y)  

## region=='EUR' & income=='Lower-middle'
X <- matrix(data_gdpr_mat_EUR_LM$gdp,ncol=1); dim(X) # 22
Y <- data_imr_mat_EUR_LM$imr; length(Y)  
M <- data_ghed_EUR_LM[,6:ncol(data_ghed_EUR_LM)]; dim(M)
S <- matrix(data_ghed_EUR_LM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 22 by 54
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_EUR_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_EUR_LM_POEM$value[output_EUR_LM_POEM$term == "pval_hdmm"] # 0.0003
output_EUR_LM_POEM$value[output_EUR_LM_POEM$term == "pval_pe"]   # 0.0003 
colnames(M)[attr(output_EUR_LM_POEM, "active_mediators")] # null

## region=='EUR' & income=='Upper-middle'
X <- matrix(data_gdpr_mat_EUR_UM$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_EUR_UM$imr; length(Y)  

## region=='EUR' & income=='High' 
X <- matrix(data_gdpr_mat_EUR_High$gdp,ncol=1); dim(X) # 176
Y <- data_imr_mat_EUR_High$imr; length(Y)  # 176
M <- data_ghed_EUR_High[,6:ncol(data_ghed_EUR_High)]; dim(M)
S <- matrix(data_ghed_EUR_High[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 176 by 56
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_EUR_High_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                            scale = FALSE,
                                            lambda_grid = lamb_grid,
                                            lambda_grid_reduced = lamb_grid0)
output_EUR_High_POEM$value[output_EUR_High_POEM$term == "pval_hdmm"] # 0.0017
output_EUR_High_POEM$value[output_EUR_High_POEM$term == "pval_pe"]   # 0 
colnames(M)[attr(output_EUR_High_POEM, "active_mediators")] # "shi_che" "gge_gdp" "gghed_ppp_pc" 

## region=='SEAR' & income=='Low' -- no samples
X <- matrix(data_gdpr_mat_SEAR_Low$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_SEAR_Low$imr; length(Y)  # 0

## region=='SEAR' & income=='Lower-middle' 
X <- matrix(data_gdpr_mat_SEAR_LM$gdp,ncol=1); dim(X) # 66
Y <- data_imr_mat_SEAR_LM$imr; length(Y) 
M <- data_ghed_SEAR_LM[,6:ncol(data_ghed_SEAR_LM)]; dim(M)
S <- matrix(data_ghed_SEAR_LM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_SEAR_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                           scale = FALSE,
                                           lambda_grid = lamb_grid,
                                           lambda_grid_reduced = lamb_grid0)
output_SEAR_LM_POEM$value[output_SEAR_LM_POEM$term == "pval_hdmm"]  # 0.9963
output_SEAR_LM_POEM$value[output_SEAR_LM_POEM$term == "pval_pe"]    # 0.9963 
colnames(M)[attr(output_SEAR_LM_POEM, "active_mediators")] # null

## region=='SEAR' & income=='Upper-middle'
X <- matrix(data_gdpr_mat_SEAR_UM$gdp,ncol=1); dim(X) # 44
Y <- data_imr_mat_SEAR_UM$imr; length(Y) 
M <- data_ghed_SEAR_UM[,6:ncol(data_ghed_SEAR_UM)]; dim(M)
S <- matrix(data_ghed_SEAR_UM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # TRUE
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_SEAR_UM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                           scale = FALSE,
                                           lambda_grid = lamb_grid,
                                           lambda_grid_reduced = lamb_grid0)
output_SEAR_UM_POEM$value[output_SEAR_UM_POEM$term == "pval_hdmm"] # 0.8255
output_SEAR_UM_POEM$value[output_SEAR_UM_POEM$term == "pval_pe"]   # 0.8255 
colnames(M)[attr(output_SEAR_UM_POEM, "active_mediators")] # null

## region=='SEAR' & income=='High' -- no samples
X <- matrix(data_gdpr_mat_SEAR_High$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_SEAR_High$imr; length(Y)  

## region=='WPR' & income=='Low' -- no samples
X <- matrix(data_gdpr_mat_WPR_Low$gdp,ncol=1); dim(X) # 0
Y <- data_imr_mat_WPR_Low$imr; length(Y)  # 0

## region=='WPR' & income=='Lower-middle'
X <- matrix(data_gdpr_mat_WPR_LM$gdp,ncol=1); dim(X) # 44
Y <- data_imr_mat_WPR_LM$imr; length(Y)
M <- data_ghed_WPR_LM[,6:ncol(data_ghed_WPR_LM)]; dim(M)
S <- matrix(data_ghed_WPR_LM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 44 by 56
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_WPR_LM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_WPR_LM_POEM$value[output_WPR_LM_POEM$term == "pval_hdmm"]  # 0.1853
output_WPR_LM_POEM$value[output_WPR_LM_POEM$term == "pval_pe"]    # 0.1853 
colnames(M)[attr(output_WPR_LM_POEM, "active_mediators")]  # null
 

## region=='WPR' & income=='Upper-middle'
X <- matrix(data_gdpr_mat_WPR_UM$gdp,ncol=1); dim(X) # 88
Y <- data_imr_mat_WPR_UM$imr; length(Y) 
M <- data_ghed_WPR_UM[,6:ncol(data_ghed_WPR_UM)]; dim(M)
S <- matrix(data_ghed_WPR_UM[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 88 by 56
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid) 
output_WPR_UM_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                          scale = FALSE,
                                          lambda_grid = lamb_grid,
                                          lambda_grid_reduced = lamb_grid0)
output_WPR_UM_POEM$value[output_WPR_UM_POEM$term == "pval_hdmm"] # 0.7309
output_WPR_UM_POEM$value[output_WPR_UM_POEM$term == "pval_pe"]   # 0.7309 
colnames(M)[attr(output_WPR_UM_POEM, "active_mediators")] # null

## region=='WPR' & income=='High'
X <- matrix(data_gdpr_mat_WPR_High$gdp,ncol=1); dim(X) # 66
Y <- data_imr_mat_WPR_High$imr; length(Y) 
M <- data_ghed_WPR_High[,6:ncol(data_ghed_WPR_High)]; dim(M)
S <- matrix(data_ghed_WPR_High[,c("year")], ncol=1); dim(S)
sum(is.na(scale(M))) == 0  # FALSE
## There are some constant columns in M, which becomes NA after scaling 
M <- M[,-which(is.na(scale(M)[1,]))]
dim(M) # 66 by 45
# Error in solve.default(t(MS) %*% MS) : system is computationally singular
scale(M)[,c("cfa_che", "vfa_che")]
colnames(M)[20]
M <- M[,-20]
dim(M) # 66 by 44
ngrid = 100
lamb_grid = seq(0.1,10,length.out = ngrid) 
lamb_grid0 = seq(0.1,10,length.out = ngrid)
output_WPR_High_POEM <- pe_mediation_linear(scale(X), Y-mean(Y), scale(M), S,
                                            scale = FALSE,
                                            lambda_grid = lamb_grid,
                                            lambda_grid_reduced = lamb_grid0)
output_WPR_High_POEM$value[output_WPR_High_POEM$term == "pval_hdmm"] # 0.0004
output_WPR_High_POEM$value[output_WPR_High_POEM$term == "pval_pe"] # 0
colnames(M)[attr(output_WPR_High_POEM, "active_mediators")]   # "cfa_che"  "pvtd_gdp"

###############################################################################
