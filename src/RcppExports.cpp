// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/rsundials.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cvodes
NumericMatrix cvodes(const NumericVector& yv, const vec& times, const RObject& frhs, RObject param, const double abstol, const double reltol, const Nullable<int> integrator_, const int maxord, const int maxsteps, const vec& constraints, const RObject fjac, const int nz, const Nullable<int> rmumps_perm_, const int nroot, const RObject froot, const RObject fevent, const int Ns, NumericVector psens, NumericVector psens_bar, const IntegerVector psens_list, const RObject fsens, const RObject fsens1, const Nullable<int> sens_method_, const bool errconS);
RcppExport SEXP _rsundials_cvodes(SEXP yvSEXP, SEXP timesSEXP, SEXP frhsSEXP, SEXP paramSEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP integrator_SEXP, SEXP maxordSEXP, SEXP maxstepsSEXP, SEXP constraintsSEXP, SEXP fjacSEXP, SEXP nzSEXP, SEXP rmumps_perm_SEXP, SEXP nrootSEXP, SEXP frootSEXP, SEXP feventSEXP, SEXP NsSEXP, SEXP psensSEXP, SEXP psens_barSEXP, SEXP psens_listSEXP, SEXP fsensSEXP, SEXP fsens1SEXP, SEXP sens_method_SEXP, SEXP errconSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type yv(yvSEXP);
    Rcpp::traits::input_parameter< const vec& >::type times(timesSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type frhs(frhsSEXP);
    Rcpp::traits::input_parameter< RObject >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< const double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< const Nullable<int> >::type integrator_(integrator_SEXP);
    Rcpp::traits::input_parameter< const int >::type maxord(maxordSEXP);
    Rcpp::traits::input_parameter< const int >::type maxsteps(maxstepsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< const RObject >::type fjac(fjacSEXP);
    Rcpp::traits::input_parameter< const int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< const Nullable<int> >::type rmumps_perm_(rmumps_perm_SEXP);
    Rcpp::traits::input_parameter< const int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< const RObject >::type froot(frootSEXP);
    Rcpp::traits::input_parameter< const RObject >::type fevent(feventSEXP);
    Rcpp::traits::input_parameter< const int >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psens(psensSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psens_bar(psens_barSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type psens_list(psens_listSEXP);
    Rcpp::traits::input_parameter< const RObject >::type fsens(fsensSEXP);
    Rcpp::traits::input_parameter< const RObject >::type fsens1(fsens1SEXP);
    Rcpp::traits::input_parameter< const Nullable<int> >::type sens_method_(sens_method_SEXP);
    Rcpp::traits::input_parameter< const bool >::type errconS(errconSSEXP);
    rcpp_result_gen = Rcpp::wrap(cvodes(yv, times, frhs, param, abstol, reltol, integrator_, maxord, maxsteps, constraints, fjac, nz, rmumps_perm_, nroot, froot, fevent, Ns, psens, psens_bar, psens_list, fsens, fsens1, sens_method_, errconS));
    return rcpp_result_gen;
END_RCPP
}
// get_cnst
int get_cnst(std::string s);
RcppExport SEXP _rsundials_get_cnst(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(get_cnst(s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rsundials_cvodes", (DL_FUNC) &_rsundials_cvodes, 24},
    {"_rsundials_get_cnst", (DL_FUNC) &_rsundials_get_cnst, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rsundials(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
