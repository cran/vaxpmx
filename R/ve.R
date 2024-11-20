#' @title Vaccine efficacy estimation
#'
#' @description Calculates vaccine efficacy and confidence interval as described in Dudasova et al., 2024, BMC Med Res Methodol and Dudasova et al., 2024, NPJ Vaccines
#' 
#' @param Fit an object of class inheriting from \code{"glm"} or \code{"coxph"} representing the fitted model
#' @param Data a data frame containing the variables in the fitted model; data must include a column called "vaccine" with binary indicator of vaccination status
#' @param nboot a numeric value for number of bootstrap samples for confidence interval construction
#'
#' @return a value of vaccine efficacy \code{VE} and lower and upper bound of confidence interval \code{CI}
#' 
#' @usage
#' ve(Fit, Data, nboot = 2000)
#'
#' @examples
#' #' # Load required packages
#' library(survival)
#' 
#' # Load an example dataset
#' data(data_temp)
#'
#' # Fit logistic model relating neutralizing titer to disease status
#' logisticFit <- glm(disease_any ~ nAb1, data = data_temp, family = binomial())
#' 
#' # Fit Cox proportional hazards model relating neutralizing titer 
#' # to time to disease or end of follow-up
#' coxFit <- coxph(Surv(time_event, disease_any) ~ nAb1, data = data_temp)
#'
#' # Estimate vaccine efficacy and 95\% confidence interval based on the fitted models
#' ve(logisticFit, data_temp, nboot = 500)
#' ve(coxFit, data_temp, nboot = 500)
#' 
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom stats predict
#' @export
ve <- function(Fit, Data, nboot = 2000){
  vaccine <- NULL
  Data.vaccinated <- Data %>% filter(vaccine == 1)
  Data.control <- Data %>% filter(vaccine == 0)
  if (as.character(Fit$call)[1] == "glm"){
    VE <- (1-(sum(predict(Fit,  Data.vaccinated, type="response"))/nrow(Data.vaccinated))/
             (sum(predict(Fit,  Data.control, type="response"))/nrow(Data.control))) * 100
    efficacySet <- glmParametricSampling(Fit, nboot, Data.vaccinated, Data.control)
    CI <- lapply(EfficacyCI(efficacySet),"*", 100)  
  } else if (as.character(Fit$call)[1] =="coxph"){
    VE <- (1-(sum(predict(Fit, Data.vaccinated, type = "risk"))/nrow(Data.vaccinated))/
             (sum(predict(Fit, Data.control, type = "risk"))/nrow(Data.control))) * 100
    efficacySet <- coxphParametricSampling(Fit, nboot, Data.vaccinated, Data.control)
    CI <- lapply(EfficacyCI(efficacySet),"*", 100)
  } else {
    return(NULL)
  }
  return (list(VE = VE,
               CI = list(LB = CI$CILow, UB = CI$CIHigh)))
}

#' @title Accounting for the uncertainty on the fitted \code{"glm"} model and observed data
#'
#' @description \code{glmParametricSampling} is used for vaccine efficacy confidence interval construction. 
#' It provides a vector of vaccine efficacy values, with length of \code{nboot}. 95\% confidence interval, defined by 2.5th and 97.5th percentile of this vector,
#' accounts for the uncertainty on the model fit (via parametric resampling of the posterior distribution of the model parameters) and observed data (via bootstrapping).
#' 
#' @param Fit an object of class inheriting from \code{"glm"} representing the fitted model
#' @param nboot a numeric value for number of bootstrap samples for confidence interval construction
#' @param Data.vaccinated a data frame for the vaccinated group, containing the variables in the fitted model; data must include a column called "vaccine" with binary indicator of vaccination status
#' @param Data.control a data frame for the control group, containing the variables in the fitted model; data must include a column called "vaccine" with binary indicator of vaccination status
#'
#' @return a vector of vaccine efficacy values \code{VE_set}, with length of \code{nboot}
#' 
#' @usage
#' glmParametricSampling(Fit, nboot = 2000, Data.vaccinated, Data.control)
#'
#' @examples
#' # Load required packages
#' library(dplyr)
#' 
#' # Load an example dataset
#' data(data_temp)
#' Data.vaccinated <- filter(data_temp, vaccine == 1)
#' Data.control <- filter(data_temp, vaccine == 0)
#'
#' # Fit logistic model relating neutralizing titer to disease status
#' logisticFit <- glm(disease_any ~ nAb1, data = data_temp, family = binomial())
#'
#' # Estimate 95\% confidence interval of vaccine efficacy based on the fitted model
#' efficacySet <- glmParametricSampling(logisticFit, nboot = 500, Data.vaccinated, Data.control)
#' CI <- lapply(EfficacyCI(efficacySet),"*", 100)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats vcov
#' @importFrom stats predict 
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @export
glmParametricSampling <- function (Fit, nboot = 2000, Data.vaccinated = NULL, Data.control = NULL){
  d <- data.frame(mvrnorm(n=nboot,
                          mu = Fit$coefficients,
                          Sigma = vcov(Fit)))
  Fit_util <- Fit
  VE_set<- c(0)
  for (i in 1:nboot){
    k <- ncol(d)
    for(j in 1:k){
      Fit_util$coefficients[[j]]<-d[i,j]
    }
    bootData.vaccinated <- Data.vaccinated[sample(seq_len(nrow(Data.vaccinated)), nrow(Data.vaccinated), replace=TRUE), ]
    bootData.control <- Data.control[sample(seq_len(nrow(Data.control)), nrow(Data.control), replace=TRUE), ]
    VE_set[i] <- 1-(sum(predict(Fit_util, bootData.vaccinated, type="response"))/nrow(bootData.vaccinated))/
      (sum(predict(Fit_util, bootData.control, type="response"))/nrow(bootData.control))
  }
  return(VE_set)
}

#' @title Accounting for the uncertainty on the fitted \code{"coxph"} model and observed data
#'
#' @description \code{coxphParametricSampling} is used for vaccine efficacy confidence interval construction. 
#' It provides a vector of vaccine efficacy values, with length of \code{nboot}. 95\% confidence interval, defined by 2.5th and 97.5th quantile of this vector,
#' accounts for the uncertainty on the model fit (via parametric resampling of the posterior distribution of the model parameters) and observed data (via bootstrapping).
#' 
#' @param Fit an object of class inheriting from \code{"coxph"} representing the fitted model
#' @param nboot a numeric value for number of bootstrap samples for confidence interval construction
#' @param Data.vaccinated a data frame for the vaccinated group, containing the variables in the fitted model
#' @param Data.control a data frame for the control group, containing the variables in the fitted model
#'
#' @return a vector of vaccine efficacy values \code{VE_set}, with length of \code{nboot}
#' 
#' @usage
#' coxphParametricSampling(Fit, nboot = 2000, Data.vaccinated, Data.control)
#'
#' @examples
#' # Load required packages
#' library(dplyr)
#' library(survival)
#' 
#' # Load an example dataset
#' data(data_temp)
#' Data.vaccinated <- filter(data_temp, vaccine == 1)
#' Data.control <- filter(data_temp, vaccine == 0)
#'
#' # Fit Cox proportional hazards model relating neutralizing titer 
#' # to time to disease or end of follow-up
#' coxFit <- coxph(Surv(time_event, disease_any) ~ nAb1, data = data_temp)
#'
#' # Estimate 95\% confidence interval of vaccine efficacy based on the fitted model
#' efficacySet <- coxphParametricSampling(coxFit, nboot = 500, Data.vaccinated, Data.control)
#' CI <- lapply(EfficacyCI(efficacySet),"*", 100)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats vcov
#' @importFrom stats predict 
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @export
coxphParametricSampling <- function (Fit, nboot = 2000, Data.vaccinated = NULL, Data.control = NULL){
  d <- data.frame(mvrnorm(n=nboot,
                          mu = Fit$coefficients,
                          Sigma = vcov(Fit)))
  Fit_util <- Fit
  VE_set<- c(0)
  for (i in 1:nboot){
    k <- ncol(d)
    for(j in 1:k){
      Fit_util$coefficients[[j]]<-d[i,j]
    }
    bootData.vaccinated <- Data.vaccinated[sample(seq_len(nrow(Data.vaccinated)), nrow(Data.vaccinated), replace=TRUE), ]
    bootData.control <- Data.control[sample(seq_len(nrow(Data.control)), nrow(Data.control), replace=TRUE), ]
    VE_set[i] <- (1-(sum(predict(Fit_util, bootData.vaccinated, type = "risk"))/nrow(bootData.vaccinated))/
                    (sum(predict(Fit_util,bootData.control, type = "risk"))/nrow(bootData.control)))
  }
  return(VE_set)
}

#' @title Efficacy summary (mean, median, confidence intervals)
#'
#' @description
#' Function summarizes efficacy statistics (mean, median, confidence intervals) based on the set of estimated efficacy values and chosen condfidence interval.
#'
#' @param efficacySet numeric vector - vector of estimated efficacy values
#' @param ci numeric - required confidence level
#'
#' @return
#' named list - mean, median, CILow, CIHigh
#'
#' @usage
#' EfficacyCI(efficacySet, ci = 0.95)
#'
#' @examples
#' # Load required packages
#' library(dplyr)
#' 
#' # Load an example dataset
#' data(data_temp)
#' Data.vaccinated <- filter(data_temp, vaccine == 1)
#' Data.control <- filter(data_temp, vaccine == 0)
#'
#' # Fit logistic model relating neutralizing titer to disease status
#' logisticFit <- glm(disease_any ~ nAb1, data = data_temp, family = binomial())
#'
#' # Estimate 95\% confidence interval of vaccine efficacy based on the fitted model
#' efficacySet <- glmParametricSampling(logisticFit, nboot = 500, Data.vaccinated, Data.control)
#' EfficacyCI(efficacySet)
#'
#' @details
#' Confidence intervals are calculated using quantiles of estimated efficacy values.
#'
#' @importFrom stats median
#' @importFrom stats quantile 
#' @export
EfficacyCI <- function(efficacySet, ci = 0.95) {
  CILow <- quantile(efficacySet, (1 - ci) / 2, names = F)
  CIHigh <- quantile(efficacySet, ci + (1 - ci) / 2, names = F)
  return(
    list(
      mean = mean(efficacySet),
      median = median(efficacySet),
      CILow = CILow,
      CIHigh = CIHigh
    )
  )
}