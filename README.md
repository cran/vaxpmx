
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vaxpmx

<!-- badges: start -->

<!-- badges: end -->

The goal of vaxpmx is to provide the functions for pharmacometric modeling in vaccines, specifically analyses related to correlates of protection (CoPs) and vaccine efficacy/effectiveness predictions using a CoP.

## Installation

You can install the released version of vaxpmx from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("vaxpmx")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(vaxpmx)
library(survival)

# Load an example dataset
data(data_temp)

# Fit logistic model relating neutralizing titer to disease status
logisticFit <- glm(disease_any ~ nAb1, data = data_temp, family = binomial())

# Fit Cox proportional hazards model relating neutralizing titer to time to disease or end of follow-up
# coxFit <- coxph(Surv(time_event, disease_any) ~ nAb1, data = data_temp)

# Estimate vaccine efficacy and 95\% confidence interval based on the fitted models
ve(logisticFit, data_temp, nboot = 500)
ve(coxFit, data_temp, nboot = 500)
```
