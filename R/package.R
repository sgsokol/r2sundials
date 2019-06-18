#' r2sundials
#'
#' Rcpp wrapper for SUNDIALS CVODE C code solving ordinary differential
#'
#' @docType package
#' @author Serguei Sokol
#' @import Rcpp RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @useDynLib r2sundials
#' @name r2sundials
NULL
# set useful constants
cnsts=c(
   "CV_SUCCESS",
   "CV_BDF",
   "CV_ADAMS",
   "R2SUNDIALS_EVENT_IGNORE",
   "R2SUNDIALS_EVENT_HOLD",
   "R2SUNDIALS_EVENT_STOP",
   "CV_SIMULTANEOUS",
   "CV_STAGGERED",
   "CV_STAGGERED1"
)
