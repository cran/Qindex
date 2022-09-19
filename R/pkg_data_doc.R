

#' @title Subset of Ki67 and PR Data
#' 
#' @description A subset of Ki67 cell data containing 623 patients and PR cell data containing 631 patients.
#' A total of 572 patients are in common.
#' 
#' @format
#' 
#' \describe{
#'   \item{PATIENT_ID}{\link[base]{factor}, ..}
#'   \item{tissueID}{\link[base]{character}, ..}
#'   \item{RECURRENCE}{\link[survival]{Surv}, ...}
#'   \item{Marker}{\link[base]{numeric}, ...}
#'   \item{all other columns}{}
#' }
#' 
#' @name celldata
'Ki67'
#' @rdname celldata
'PR'