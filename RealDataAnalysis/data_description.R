###################################################
# Infant Mortality Rate (IMR)
rm(list=ls())
load("data/data_GDPR_GHED_IMR.RData")
dim(data_gdpr) # 91 by 23
dim(data_imr)  # 91 by 23
dim(data_ghed) # 2002 by 62
## columns: country, code, region, income, year, 57 indicators 
## rows: 91 countries * 22 years 
###################################################

###################################################
## Under-Five Mortality Rate (U5MR)
rm(list=ls())
load("data/data_GDPR_GHED_U5MR.RData")
dim(data_gdpr) # 91 by 23
dim(data_u5mr) # 91 by 23
dim(data_ghed) # 2002 by 62
## columns: country, code, region, income, year, 57 indicators 
## rows: 91 countries * 22 years 
###################################################

###################################################
# Life Expectancy at Birth (LEB)
rm(list=ls())
load("data/data_GDPR_GHED_LEB.RData")
dim(data_gdpr) # 91 by 23
dim(data_leb)  # 91 by 23
dim(data_ghed) # 2002 by 62
## columns: country, code, region, income, year, 57 indicators 
## rows: 91 countries * 22 years 
###################################################

###################################################
#  the Prevalence of Low Birthweight Infants (LBW)
rm(list=ls())
load("data/data_GDPR_GHED_LBW.RData")
dim(data_gdpr) # 73 by 22
dim(data_lbw)  # 73 by 22
dim(data_ghed) # 1533 by 62
## columns: country, code, region, income, year, 57 indicators 
## rows: 73 countries * 21 years 
###################################################

###################################################
## the Prevalence of Undernourishment (PoU)
rm(list=ls())
load("data/data_GDPR_GHED_POU.RData")
dim(data_gdpr) # 79 by 22
dim(data_pou)  # 79 by 22
dim(data_ghed) # 1659 by 62
## columns: country, code, region, income, year, 57 indicators 
## rows: 79 countries * 21 years 
###################################################

