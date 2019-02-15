// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/rsundials.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cvode
NumericMatrix cvode(const NumericVector& yv, const NumericVector& times, const RObject& fderiv, RObject& param, const double abstol, const double reltol, RObject fjac, RObject jacmat_, std::string rmumps_perm);
RcppExport SEXP _rsundials_cvode(SEXP yvSEXP, SEXP timesSEXP, SEXP fderivSEXP, SEXP paramSEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP fjacSEXP, SEXP jacmat_SEXP, SEXP rmumps_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type yv(yvSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type times(timesSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type fderiv(fderivSEXP);
    Rcpp::traits::input_parameter< RObject& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< const double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< RObject >::type fjac(fjacSEXP);
    Rcpp::traits::input_parameter< RObject >::type jacmat_(jacmat_SEXP);
    Rcpp::traits::input_parameter< std::string >::type rmumps_perm(rmumps_permSEXP);
    rcpp_result_gen = Rcpp::wrap(cvode(yv, times, fderiv, param, abstol, reltol, fjac, jacmat_, rmumps_perm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rsundials_cvode", (DL_FUNC) &_rsundials_cvode, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_rsundials(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
