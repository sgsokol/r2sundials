#' rsundials
#'
#' Rcpp wrapper for SUNDIALS CVODE C code solving ordinary differential
#'
#' @docType package
#' @author Serguei Sokol
#' @import Rcpp RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @useDynLib rsundials
#' @name rsundials
NULL
# set useful constants
cnsts=c(
   "CV_SUCCESS",
   "CV_BDF",
   "CV_ADAMS",
   "RSUNDIALS_EVENT_IGNORE",
   "RSUNDIALS_EVENT_HOLD",
   "RSUNDIALS_EVENT_STOP",
   "RMUMPS_PERM_AMD",
   "RMUMPS_PERM_AMF",
   "RMUMPS_PERM_SCOTCH",
   "RMUMPS_PERM_PORD",
   "RMUMPS_PERM_METIS",
   "RMUMPS_PERM_QAMD",
   "RMUMPS_PERM_AUTO",
   "CV_SIMULTANEOUS",
   "CV_STAGGERED",
   "CV_STAGGERED1"
)
