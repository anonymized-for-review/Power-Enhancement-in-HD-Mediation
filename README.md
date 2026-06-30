# Power Enhancement in High-Dimensional Heterogeneous Mediation Analysis

This repository contains the code and processed data used for the simulation studies and empirical analysis in the paper "Power Enhancement in High-Dimensional Heterogeneous Mediation Analysis", submitted to the *Journal of the American Statistical Association*.

The repository is organized into two main components:

- `Simulation/`: simulation scripts for continuous, dichotomous, and count outcomes
- `RealDataAnalysis/`: empirical analysis scripts and processed datasets used in the paper

An R package, **`POEM`**, is developed to facilitate practical implementations and is used throughout the codebase. 

- Package link: `[https://drive.google.com/drive/folders/1-hNpapZ51IE0m_flELAqohne1CPJrwgD?usp=sharing]`

[Note: To comply with JASA’s double-blind review policy, we have not yet submitted the package to R CRAN, but we plan to do so after the review process is complete.]

## Instructions to install the `POEM` package: 

`download.file( "https://drive.google.com/uc?export=download&id=1QHfmx55aJZJ8MKZSIKXCDCIpxJoLcSUH", "POEM_0.1.0.tar.gz", mode = "wb")`

`remotes::install_local("POEM_0.1.0.tar.gz", build_vignettes = TRUE, dependencies = TRUE)`

## Repository Structure

### `Simulation/`

This folder contains the main scripts for reproducing the simulation studies of the proposed power-enhanced heterogeneous mediation models.

- `code-main-simulation-linear.R`: simulation study for continuous outcomes
- `code-main-simulation-logistic.R`: simulation study for dichotomous outcomes
- `code-main-simulation-poisson.R`: simulation study for count outcomes
- `code-utils-comparison-methods-linear.R`: helper functions and comparison methods for continuous outcomes
- `code-utils-comparison-methods-logistic.R`: helper functions and comparison methods for dichotomous outcomes
- `code-utils-comparison-methods-poisson.R`: helper functions and comparison methods for count outcomes
- `README.md`: additional notes specific to the simulation studies

### `RealDataAnalysis/`

This folder contains the empirical example based on World Bank Open Data and the WHO Global Health Expenditure Database, along with processed datasets prepared for the paper.

- `code-PE-mediation-RealDataExample.R`: main real-data analysis script
- `code-utils-HDMM.R`: helper functions used in the real-data analysis
- `data/`: processed `.RData` files for five public-health outcomes
- `data_description.R`: quick data loading and dimension checks
- `README.md`: dataset descriptions and source notes

## Data Sources

The raw data used to construct the empirical example come from:

- World Bank Open Data: https://data.worldbank.org/
- WHO Global Health Expenditure Database: https://apps.who.int/nha/database/

This repository includes processed versions of the datasets used in the paper so the analysis can be reproduced directly. If you use these processed data directly from the repository, please cite the paper.


