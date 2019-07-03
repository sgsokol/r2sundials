## Introduction

Package `r2sundials` is an RcppArmadillo wrapper for a well known and
wide spread library
[SUNDIALS/CVODES](https://computation.llnl.gov/sites/default/files/public/cvs_guide.pdf)
from LLNL written in C. It provides an access from R to some basic
features of `cvodes` routine from this library which include:

  - solving user defined ODE via user provided functions calculating
    ODE’s right hand side (rhs);
  - calculating forward sensitivities to parameters on which ODE
    solution depends ;
  - setting key parameters for ODE solving method like explicit or
    implicit time scheme, minimal or maximal steps, error order etc.

Note that user has to install
[SUNDIALS/CVODES](https://computation.llnl.gov/projects/sundials/cvodes)
by his own means (cf. [Install](#install) section)

## Why another ODE solver for R?

The question is legitimate as there is already very furnished package
[deSolve](https://cran.r-project.org/package=deSolve). There were also
[Rsundials](https://cran.r-project.org/web/packages/Rsundials/index.html)
package which is archived since 2017-03-26. There is also a fresh
wrapper to the same library
[sundialr](https://cran.r-project.org/package=sundialr). So why a new
wrapper? Let see what are the novelties brought by `r2sundials` compared
to still active alternatives: `deSolve` and `sundialr`.

### deSolve

Compared to `deSolve`, `r2sundials` provides possibilities:

  - to do forward sensitivity calculations for all or selected
    parameters;
  - to write rhs in Rcpp where calculated values are stored “in-place”
    thus avoiding frequent memory reallocation;
  - to write dense or sparse Jacobian Rcpp functions based on the same
    principle, i.e. “in-place” storage;
  - to use [rmumps](https://cran.r-project.org/package=rmumps) package
    for solving underlying sparse linear systems;
  - to have much more flexible root finding and handling.
  - to see statistics of ODE methods (call number for rhs routines,
    number of Jacobian calculations and so on).

### sundialr

Compared to `sundialr`, `r2sundials` provides:

  - more complete access to fine tuning of `cvodes` methods;
  - more complete parameter infrastructure. In `r2sundials` parameters
    passed to callback functions can be of any R type (vector, list,
    environment, …), not only numeric vector as in `sundialr`;
  - sensitivity calculations possibly done with the help of user
    provided functions;
  - if automatic sensitivity calculations is used, it can be done on a
    selection of parameters, not necessarily on the totality of
    parameters.
  - Jacobian (dense or sparse) calculated with possibly user provided
    functions;
  - sparse calculations made with `rmumps` package;
  - root finding and handling;
  - some statistics of ODE methods (call number for rhs routines, number
    of Jacobian calculations and so on).

## Install

The wrapper itself can be installed as any other CRAN or github package.

However, prior to this, SUNDIALS CVODES library has to be installed on
the user’s system and corresponding locations of include and library
directories passed to installation routines. After downloading
[sundials/cvodes](https://computation.llnl.gov/projects/sundials/download/cvodes-4.1.0.tar.gz),
user has to configure and build it using cmake system. It is highly
advised to use the same compiler for compiling cvodes libraries as was
used for compiling R itself. It can be detected with a command `R CMD
config CC`. A mandatory configuration requirement for cvodes is to set
the index size to 32 bits. At the end of the configuration step,
corresponding line in `CMakeCache.txt` file in the build directory must
look like:

``` 
 SUNDIALS_INDEX_SIZE:STRING=32
```

Failing to do so will result in no compilation of `r2sundials`. User can
also to configure `cvodes` to use local BLAS/LAPACK libraries but this
point is out of our scope.

The rest of the installation procedure will be considered on two system
types: Linux-like and Windows. Let note installation directory of
`cvodes` by `<cvodes_prefix>`.

### Linux-like system

User can install `r2sundials` with

``` r
install.packages("r2sundials", configure.args="--with-cvodes_include=<cvodes_prefix>/include --with-cvodes_libs=<cvodes_prefix>/lib")
```

here `<cvodes_prefix>` has to be replaced by its actual value. Note that
if `<cvodes_prefix>=/usr/local` then installation incantation can be as
simple as

``` r
install.packages("r2sundials")
```

as this location is used as default.

### Windows

After installing `cvodes`, user has to permanently define 3 environment
variables:

    cvodes_include=<cvodes_prefix>\include
    cvodes_libs=<cvodes_prefix>\lib
    PATH=%PATH%;%cvodes_libs%

For example, on my windows machine, `cvodes` was installed in
`d:\\local_soft\cvodes` directory so that a file `cvodes.h` can be found
in `d:\\local_soft\cvodes\include\cvodes\cvodes.h`. In this case, I
could define (using a command `setx`):

    > setx cvodes_include d:\\local_soft\cvodes\include
    > setx cvodes_libs d:\\local_soft\cvodes\lib
    > setx PATH %PATH%;d:\\local_soft\cvodes\lib

This is to be done only once, at the first installation of `r2sundials`.

User has to restart `CMD.EXE` to see these new variables to be taken
into account and to use them in the rest of the installation procedure.
After that in R, do:

``` r
install.packages("r2sundials")
```

### Version note

`r2sundials` was developed and tested with CVODES version 4.1.0 released
on 2019-02-13. It can happen that `r2sundials` works with other versions
of CVODES but it was not tested. The versioning scheme of `r2sundials`
is based on the version of CVODES used during developments extended with
one number proper to `r2sundials`. For example, the first `r2sundials`
release has the version 4.1.0-1.

## Example

Let solve a very simple ODE \(y'(t)=-ν(y(t)-a)\), with \(ν=2\), \(a=1\)
and \(y(0)=0\) on a time interval \([0, 3]\). This equation describes an
exponential transition between two states 0 and \(a\) with a rate \(ν\).

With rhs written in R, we can do:

``` r
library(r2sundials)
ti=seq(0, 3, length.out=101)
p=c(nu=2, a=1)
y0=0
frhs=function(t, y, p, psens) -p["nu"]*(y-p["a"])
res=r2sundials::cvodes(y0, ti, frhs, param=p)
# compare with analytical solution
stopifnot(diff(range(p["a"]-exp(-p["nu"]*ti) - res)) < 1.e-6)
# see stats
print(attr(res, "stats"))
```

    ##                   NumSteps                NumRhsEvals 
    ##                        112                        136 
    ##           NumLinSolvSetups            NumErrTestFails 
    ##                         23                          3 
    ##         NumNonlinSolvIters     NumNonlinSolvConvFails 
    ##                        133                          0 
    ##                NumJacEvals             NumLinRhsEvals 
    ##                          2                          2 
    ##                  NumGEvals            SensNumRhsEvals 
    ##                          0                          0 
    ##            NumRhsEvalsSens       SensNumLinSolvSetups 
    ##                          0                          0 
    ##        SensNumErrTestFails     SensNumNonlinSolvIters 
    ##                          0                          0 
    ## SensNumNonlinSolvConvFails 
    ##                          0

The same problem solved with RcppArmadillo rhs can look like:

``` r
library(RcppXPtrUtils)
Sys.setenv(PKG_CXXFLAGS=paste0("-I ", gsub("\\", "/", readLines(system.file("cvodes.txt", package="r2sundials"))[1L], fixed=TRUE)))
ti=seq(0, 3, length.out=101)
p=c(nu=2, a=1)
y0=0
# next step can take few seconds at the first execution as it will compile C++ code.
ptr_exp=cppXPtr(code='
int rhs_exp(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
  NumericVector p(param);
  ydot[0] = -p["nu"]*(y[0]-p["a"]);
  return(CV_SUCCESS);
}
', depends=c("RcppArmadillo","r2sundials","rmumps"),
 includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)

res2=r2sundials::cvodes(y0, ti, ptr_exp, param=p)
stopifnot(diff(range(res2 - res)) < 1.e-14)
```

For more examples, see `?r2sundials::cvodes`.
