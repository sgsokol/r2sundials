#ifndef RSUNDIALS_H
#define RSUNDIALS_H

#include <cvode/cvode.h> // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h> // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_sparse.h> // access to sparse SUNMatrix
#include <sundials/sundials_types.h> // defs. of realtype, sunindextype

#include <sunlinsol_rmumps.h>

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define getmem(y, f) {(y)=(f); if ((y) == NULL) stop("no memory for " #y);}
#define check_retval(expr) {int retval=(expr); if (retval != CV_SUCCESS) stop("bad returned value: %d\ncall: %s", retval, #expr);}

// define a type for user supplied function rhs
typedef int (*rsunDerivFn)(double t, const vec &y, vec &ydot, RObject &param);
typedef int (*rsunJacFn)(double t, const vec &y, vec &ydot, mat &J, RObject &param);
typedef int (*rsunSpJacFn)(double t, vec &y, vec &ydot, uvec &i, uvec &p, vec &v, int n, int nz, RObject &param);

int derivwrap(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int jacsps(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int jacwrap(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int spjacwrap(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// error handler
void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data);
#endif
