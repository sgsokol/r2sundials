#ifndef RSUNDIALS_H
#define RSUNDIALS_H

#include <cvode/cvode.h> // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h> // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_sparse.h> // access to sparse SUNMatrix
#include <sundials/sundials_types.h> // defs. of realtype, sunindextype

#define RSUNDIALS_EVENT_IGNORE 0
#define RSUNDIALS_EVENT_HOLD 1
#define RSUNDIALS_EVENT_STOP -1

#include <sunlinsol_rmumps.h>

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define getmem(y, f) {(y)=(f); if ((y) == NULL) stop("no memory for " #y);}
#define check_retval(expr) {int retval=(expr); if (retval != CV_SUCCESS) stop("rsundials: call: %s\nreturned flag: %s", #expr, CVodeGetReturnFlagName(retval));}

// define a type for user supplied function rhs
typedef int (*rsunRhsFn)(double t, const vec &y, vec &ydot, RObject &param);
typedef int (*rsunJacFn)(double t, const vec &y, vec &ydot, mat &J, RObject &param);
typedef int (*rsunSpJacFn)(double t, vec &y, vec &ydot, uvec &i, uvec &p, vec &v, int n, int nz, RObject &param);
typedef int (*rsunRootFn)(double t, const vec &y, vec &vroot, RObject &param);
typedef int (*rsunEventFn)(double t, const vec &y, vec &ynew, const ivec &rootsfound, RObject &param);

int rhswrap(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int jacwrap(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int spjacwrap(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int rootwrap(realtype t, N_Vector y, realtype *rootout, void *user_data);

// error handler
void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data);
#endif
