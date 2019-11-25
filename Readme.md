---
title: "ODE and forward sensitivity solving with `r2sundials`"
author: Serguei SOKOL
institute:  INRA / TBI / Mathematics cell
date: "November 25, 2019"
output:
  md_document:
    variant: gfm
---

## Introduction
Package `r2sundials` is an RcppArmadillo wrapper for a well known and wide spread library [SUNDIALS/CVODES](https://computation.llnl.gov/sites/default/files/public/cvs_guide.pdf) from LLNL written in C. It provides an access from R to some basic features of `cvodes` module from this library which include:

 - solving real valued, user defined ODE via user provided functions calculating ODE's right hand side (rhs);
 - calculating first order forward sensitivities to parameters on which ODE solution depends ;
 - setting key parameters for ODE solving method like explicit or implicit time scheme, minimal or maximal steps, error order etc.
 
Note that user has to install [SUNDIALS/CVODES](https://computation.llnl.gov/projects/sundials/cvodes) by his own means (cf. [Install] section)

## Why another ODE solver for R?

The question is legitimate as there is already very furnished package [deSolve](https://cran.r-project.org/package=deSolve). It was also published [Rsundials](https://cran.r-project.org/web/packages/Rsundials/index.html) package but it is archived now since 2017-03-26. There is a fresh wrapper to the same library [sundialr](https://cran.r-project.org/package=sundialr). You can find even more packages dedicated to ODE solving in this [task view](https://cran.r-project.org/web/views/DifferentialEquations.html). So why a new wrapper?
Let see what are the novelties brought by `r2sundials` compared to most close alternatives: `deSolve` and `sundialr`.

### deSolve
Compared to `deSolve`, `r2sundials` provides possibilities:

 - to do forward sensitivity calculations for all or selected parameters;
 - to write users callback functions (rhs, Jacobian, ...) in Rcpp(Armadillo) where calculated values are stored "in-place" thus avoiding frequent memory reallocation;
 - to pass *directly* to Rcpp functions parameter of any R type (vector, list, environment, ...);
 - to use [rmumps](https://cran.r-project.org/package=rmumps) package for solving underlying sparse linear systems;
 - to have much more flexible root finding and handling based on user callback functions.
 - to get statistics on ODE method used (call number for rhs routines, number of Jacobian calculations and so on).
 
### sundialr
Compared to `sundialr`, `r2sundials` provides:

 - more complete access to fine tuning of `cvodes` methods;
 - more complete parameter infrastructure. In `r2sundials`, parameters passed to callback functions can be of any R type (vector, list, environment, ...), not only numeric vector as in `sundialr`;
 - sensitivity calculations possibly done with the help of user provided functions (and not only with internal sundials procedure);
 - sensitivity calculations can be done on a selection of parameters, not necessarily on the totality of parameters.
 - Jacobian (dense or sparse) calculated with possibly user provided functions;
 - sparse system solving is made with `rmumps` package;
 - root finding and handling;
 - some statistics of ODE methods (call number for rhs routines, number of Jacobian calculations and so on).
 - (as of time of this writing, 2019-11-25) more thorough memory management which ensures that sundials' allocated memory is freed in due way, no matter what C++ exception and in what moment could happen. This avoids potential memory leaking problem in case of frequent package use during the same session.
 
## Install
The package itself can be installed as any other CRAN or github package.

However, prior to this, SUNDIALS CVODES library has to be installed on the user's system and corresponding locations of include and library directories passed to installation routines.
After downloading [sundials/cvodes](https://computation.llnl.gov/projects/sundials/download/cvodes-5.0.0.tar.gz), user has to configure and build it using cmake system. It is highly advised to use the same compiler for compiling cvodes libraries as was used for compiling R itself. It can be detected with a command `R CMD config CC`. A suggested configuration option for `cvodes` is to set the index size to 32 bits. If this option is retained, at the end of the configuration step, corresponding line in `CMakeCache.txt` file in the build directory will look like:

```
 SUNDIALS_INDEX_SIZE:STRING=32
```

Making index size set to 32 bits will prevent additional conversion between `int64` and `int32` types each time that a sparse Jacobian is calculated by a user supplied function.
So this option can provide some savings in calculation time but is not mandatory. User can leave the cvodes' default index size set to int64.
User can also configure `cvodes` to use local BLAS/LAPACK libraries but this point is out of our scope.

The rest of the installation procedure will be considered on two system types: Linux-like and Windows.
Let note  installation directory of `cvodes` by `<cvodes_prefix>`.


### Linux-like system

User can install `r2sundials` with a command similar to:


```r
install.packages("r2sundials", configure.args="--with-cvodes_include=<cvodes_prefix>/include --with-cvodes_libs=<cvodes_prefix>/lib")
```

here `<cvodes_prefix>` has to be replaced by its actual value. Note that if `<cvodes_prefix>=/usr/local` then installation incantation can be as simple as 


```r
install.packages("r2sundials")
```

as this location is used as default. In this case header cvodes files are loocked up in `/usr/local/include` and corresponding libraries are supposed to be in `/usr/local/lib64`.

### Windows
After installing `cvodes`, user has to permanently define 3 environment variables:
```
cvodes_include=<cvodes_prefix>\include
cvodes_libs=<cvodes_prefix>\lib
PATH=%PATH%;%cvodes_libs%
```

For example, on a windows machine, `cvodes` was installed in `d:\\local_soft\sundials` directory so that a file `cvodes.h` can be found in `d:\\local_soft\sundiales\include\cvodes\cvodes.h`. In this case, we could define (using a command `setx`):

```
> setx cvodes_include d:\\local_soft\sundials\include
> setx cvodes_libs d:\\local_soft\sundials\lib
> setx PATH %PATH%;d:\\local_soft\sundials\lib
```
This is to be done only once, at the first installation of `sundials`.

User has to restart `CMD.EXE` to see these new variables to be taken into account and to use them in the rest of the installation procedure.
After that, do in R:

```r
install.packages("r2sundials")
```

### Version note
`r2sundials` was developed and tested with CVODES version 5.0.0 released in October 2019. It can happen that `r2sundials` works with other versions of CVODES but it was not tested. The versioning scheme of `r2sundials` is based on the version of CVODES used during developments extended with one number proper to `r2sundials`. For example, `r2sundials` can have a version 5.0.0-2.

## Example
Let solve a very simple ODE $y'(t)=-ν(y(t)-a)$, with $ν=2$, $a=1$ and $y(0)=0$ on a time interval $[0, 3]$. This equation describes an exponential transition between two states 0 and $a$ with a rate $ν$.

With rhs written in R, we can do:

```r
library(r2sundials)
```

```
## Loading required package: rmumps
```

```
## 
## Attaching package: 'r2sundials'
```

```
## The following object is masked from 'package:rmumps':
## 
##     get_cnst
```

```r
ti=seq(0, 3, length.out=101) # set time grid
p=c(nu=2, a=1) # set parameter vector
y0=0 # set initial condition
frhs=function(t, y, p, psens) -p["nu"]*(y-p["a"]) # set rhs functions
res=r2sundials::r2cvodes(y0, ti, frhs, param=p) # solve ODE
# compare with analytical solution
stopifnot(diff(range(p["a"]-exp(-p["nu"]*ti) - res)) < 1.e-6)
# see stats
print(attr(res, "stats"))
```

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
```

The same problem solved with RcppArmadillo rhs can look like:

```r
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

res2=r2sundials::r2cvodes(y0, ti, ptr_exp, param=p)
stopifnot(diff(range(res2 - res)) < 1.e-14)
```

For more examples, see `?r2sundials::r2cvodes`.
