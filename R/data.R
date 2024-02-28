#' 	Example of a hypothetical vaccine clinical trial data set 
#'
#' A dataset containing immunogenicity data, and clinical outcome data in the vaccinated and control groups. The dataset is provided in the form of a data frame.
#'
#' @format Data frame:
#' \describe{
#'   \item{ID}{identification of subjects}
#'   \item{nAb1}{value of neutralizing titer for serotype 1}
#'   \item{nAb2}{value of neutralizing titer for serotype 2}
#'   \item{group}{binary indicator of a baseline demographic characteristics of interest}
#'   \item{vaccine}{binary indicator of treatment arm, with value 1 in vaccinated and 0 in control subjects}
#'   \item{type_disease}{serotype of disease}
#'   \item{disease_any}{binary indicator of disease caused by any serotype}
#' }
"data_temp"

