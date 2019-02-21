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
RSUNDIALS_EVENT_IGNORE=0
RSUNDIALS_EVENT_HOLD=1
RSUNDIALS_EVENT_STOP=-1
