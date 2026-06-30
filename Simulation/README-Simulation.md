# README: Simulation Studies

This repository contains the implementation details of simulation studies for the paper "Power Enhancement in High-Dimensional Heterogeneous Mediation Analysis", submitted to the *Journal of the American Statistical Association*.

This folder contains R scripts for simulation studies of the proposed `PE-HD(G)MM` method for:

- continuous outcomes
- dichotomous outcomes
- count outcomes

## Files

- `code-main-simulation-linear.R`: main script of HD mediation models for continuous outcomes
- `code-main-simulation-logistic.R`: main script of HD mediation models for dichotomous outcomes
- `code-main-simulation-poisson.R`: main script of HD mediation models for count outcomes
- `code-utils-comparison-methods-linear.R`: helper functions for comparisons of  HD linear mediation models 
- `code-utils-comparison-methods-logistic.R`: helper functions for comparisons of HD logistic mediation models 
- `code-utils-comparison-methods-poisson.R`: helper functions for comparisons of HD Poisson mediation models

## Requirements

Required R packages for the proposed `PE-HD(G)MM` implementation:

- `POEM`
- `stats`
- `glmnet`
- `ncvreg`

Required R packages for the comparison methods:

- `freebird`
- `globaltest`
