# README: Real Data Analysis

This repository contains the implementation details of the real data analysis for the paper "Power Enhancement in High-Dimensional Heterogeneous Mediation Analysis", submitted to the *Journal of the American Statistical Association*.

The raw data are publicly available at the World Bank Open Data [https://data.worldbank.org/] and the World Health Organization (WHO) Global Health Expenditure Database [https://apps.who.int/nha/database/]. In this repository, the `data` folder contains the processed datasets used in our paper. If you use these processed data directly from the repository, please cite the paper.


# Data Description

In this task, we study the question of how healthcare expenditures mediate the relationship between economic growth and public health outcomes in different countries, regions, and income groups.


The `data` folder contains five outcome-specific datasets:

- `data_GDPR_GHED_IMR.RData` Infant Mortality Rate (IMR)
- `data_GDPR_GHED_U5MR.RData` Under-Five Mortality Rate (U5MR)
- `data_GDPR_GHED_LEB.RData` Life Expectancy at Birth (LEB)
- `data_GDPR_GHED_LBW.RData` Prevalence of Low Birthweight Infants (LBW)
- `data_GDPR_GHED_POU.RData` Prevalence of Undernourishment (POU)

Each `.RData` file contains three objects:

- `data_gdpr`: the annual percentage growth rate of GDP per capita (exposure variable)
- `data_<outcome>`: the public health outcome of interest (outcome variable)
- `data_ghed`: various indicators of health expenditure (mediators)

# Details for Each Dataset

`data_GDPR_GHED_IMR.RData`

- `data_gdpr` 91 rows (countries), 23 columns (`Country.Code`, and yearly values from 2000 to 2021) 
- `data_imr` 91 rows (countries), 23 columns (`Country.Code`, and yearly values from 2000 to 2021)
- `data_ghed` 2002 rows (91 countries * 22 years) and 62 columns (`country`, `code`, `region`, `income`, `year`, and 57 health expenditure indicators). 

`data_GDPR_GHED_U5MR.RData`

- `data_gdpr` 91 rows (countries), 23 columns (`Country.Code`, and yearly values from 2000 to 2021)
- `data_u5mr` 91 rows (countries), 23 columns (`Country.Code`, and yearly values from 2000 to 2021)
- `data_ghed` 2002 rows (91 countries * 22 years) and 62 columns (`country`, `code`, `region`, `income`, `year`, and 57 health expenditure indicators). 

`data_GDPR_GHED_LEB.RData`

- `data_gdpr` 91 rows (countries), 23 columns (`Country.Code`, and yearly values from 2000 to 2021)
- `data_leb` 91 rows (countries), 23 columns (`Country.Code`, and yearly values from 2000 to 2021)
- `data_ghed` 2002 rows (91 countries * 22 years) and 62 columns (`country`, `code`, `region`, `income`, `year`, and 57 health expenditure indicators). 


`data_GDPR_GHED_LBW.RData`

- `data_gdpr` 73 rows (countries), 22 columns (`Country.Code`, and yearly values from 2000 to 2020)
- `data_lbw` 73 rows (countries), 22 columns (`Country.Code`, and yearly values from 2000 to 2020)
- `data_ghed` 1533 rows (73 countries * 21 years) and 62 columns (`country`, `code`, `region`, `income`, `year`, and 57 health expenditure indicators). 


`data_GDPR_GHED_POU.RData`

- `data_gdpr` 79 rows (countries), 22 columns (`Country.Code`, and yearly values from 2001 to 2021)
- `data_pou` 79 rows (countries), 22 columns (`Country.Code`, and yearly values from 2001 to 2021)
- `data_ghed` 1659 rows (79 countries * 21 years) and 62 columns (`country`, `code`, `region`, `income`, `year`, and 57 health expenditure indicators). 
