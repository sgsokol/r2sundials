context("Robertson")
yini <- c(y1=1, y2=0, y3=0)
neq=length(yini)
# parameters
parms <- c(k1 = 0.04, k2 = 3e7, k3 = 1e4)
# derivative functions
vRober <- function(y, parms) {
  dy1 <- -parms["k1"]*y[1]                + parms["k3"]*y[2]*y[3]
  dy2 <-  parms["k1"]*y[1] - parms["k2"]*y[2]*y[2] - parms["k3"]*y[2]*y[3]
  dy3 <-  parms["k2"]*y[2]*y[2]
  c(dy1, dy2, dy3)
}
Rober <- function(t, y, parms) {
  
#print(t)
#print(y)
#print(c(dy1, dy2, dy3))
  list(vRober(y, parms))
}
r_rober <- function(t, y, parms, psens) vRober(y, parms)

# -------------------------------------------------------
# run at high resolution 
# -------------------------------------------------------

times <- 10^(seq(from = -5, to = 11, by = 0.1))

# lsoda!
print (system.time(
out <- deSolve::ode(func = Rober, parms = parms, y = yini,
           times = times, atol = 1e-10, rtol = 1e-10,
           maxsteps = 1e5)
))

# now sundials cvode
# pointer to rhs function
includes="using namespace Rcpp;\nusing namespace arma;"
pfnd=cppXPtr(code='
int d_robertson(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
  NumericVector p(param);
  ydot[0] = -p["k1"]*y[0] + p["k3"]*y[1]*y[2];
  ydot[2] = p["k2"]*y[1]*y[1];
  ydot[1] = -ydot[0] - ydot[2];
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)
# pointer to dense jacobian function
pfnj=cppXPtr(code='
int jac_robertson(double t, const vec &y, vec &ydot, mat &J, RObject &param, NumericVector &psens) {
  NumericVector p(param);
  J.at(0, 0) = -p["k1"];
  J.at(1, 0) = p["k1"]; 
  J.at(2, 0) = 0.;

  J.at(0, 1) = p["k3"]*y[2];
  J.at(1, 1) = -p["k3"]*y[2]-2*p["k2"]*y[1];
  J.at(2, 1) = 2*p["k2"]*y[1];

  J.at(0, 2) = p["k3"]*y[1];
  J.at(1, 2) = -p["k3"]*y[1];
  J.at(2, 2) = 0.;
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)
# pointer to sparse jacobian function
# illustrates usage of named components of param vector
pfnspj=cppXPtr(code='
int spjac_robertson(double t, const vec &y, vec &ydot, uvec &ir, uvec &pj, vec &v, int n, int nz, RObject &param, NumericVector &psens) {
  if (nz < 7)
    stop("spjac_robertson: not enough room for non zeros, must have at least 7, instead got %d", nz);
  NumericVector prm(param);
  int i=0;
  pj[0] = 0; // init pj
  // first column
  ir.at(i) = 0;
  v.at(i++) = -prm["k1"];
  ir.at(i) = 1;
  v.at(i++) = prm["k1"]; 
  pj[1] = i;
  // second column
  ir.at(i) = 0;
  v.at(i++) = prm["k3"]*y[2];
  ir.at(i) = 1;
  v.at(i++) = -prm["k3"]*y[2]-2*prm["k2"]*y[1];
  ir.at(i) = 2;
  v.at(i++) = 2*prm["k2"]*y[1];
  pj[2] = i;
  // third column
  ir.at(i) = 0;
  v.at(i++) = prm["k3"]*y[1];
  ir.at(i) = 1;
  v.at(i++) = -prm["k3"]*y[1];
  pj[3] = i;
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)
# pointer to sensitivity rhs function
pfnsens1=cppXPtr(code='
int sens_robertson1(int Ns, double t, const vec &y, vec &ydot, int iS, vec &yS, vec &ySdot, RObject &param, NumericVector &p, vec &tmp1, vec &tmp2) {
  // calculate (∂f /∂y)s_i(t) + (∂f /∂p_i) for i = iS
  // (∂f /∂y)s_i(t)
//print(p);
//stop("print");
  ySdot[0] = -p["k1"]*yS[0] + p["k3"]*y[2]*yS[1] + p["k3"]*y[1]*yS[2];
  ySdot[1] = p["k1"]*yS[0] - (p["k3"]*y[2]+2*p["k2"]*y[1])*yS[1] - p["k3"]*y[1]*yS[2]; 
  ySdot[2] = 2*p["k2"]*y[1]*yS[1];
  // + (∂f /∂p_i)
  switch(iS) {
    case 0:
      ySdot[0] -= y[0];
      ySdot[1] += y[0];
      break;
    case 1:
      ySdot[1] -= y[1]*y[1];
      ySdot[2] += y[1]*y[1];
      break;
    case 2:
      ySdot[0] += y[1]*y[2];
      ySdot[1] -= y[1]*y[2];
  }
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)

# sparse Jacobian
#print (system.time(
#out1 <- rsundials::cvodes(yini, times, pfnd, param=parms, fjac=pfnspj, nz=8)
#))
# dense Jacobian + forward sensitivity 1 by 1
#print (system.time(
#out2 <- rsundials::cvodes(yini, times, pfnd, param=parms, fjac=pfnj, Ns=3, psens=parms, fsens1=pfnsens1)
#))

# bouncing ball example (to illustrate discontinuties handling)
# A ball falls from some height. At this moment, it has 0 vertical speed and some non zero horizontal speed.
# When it hits the ground, it bounce and looses ky part of its vertical speed and kx part of the horizontal one.
# Simulation should stop after 5th hit of the ground (unknown time point before running simulation).

yinib=c(x=0, y=1, vx=0.5, vy=0) # initial values and their names
paramb=c(g=9.81, kx=0.1, ky=0.3, nbounce=5) # usefull parameters
timesb=seq(0, 3, length.out=101) # time points
# pointer to rhs function
pball=cppXPtr(code='
int d_ball(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
  NumericVector p(param);
  ydot[0] = y[2];
  ydot[1] = y[3];
  ydot[2] = y[1] > 0 ? 0. : -y[2]; // falling till y=0 then damping
  ydot[3] = y[1] > 0 ? -p["g"] : -y[3]; // falling till y=0 then damping
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)
# pointer to root function
proot=cppXPtr(code='
int root_ball(double t, const vec &y, vec &vroot, RObject &param) {
  vroot[0] = y[1]; // y==0
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)
# pointer to event handler function
Sys.setenv(PKG_CXXFLAGS="-I/home/sokol/dev/R/rcpp-pkgs/rsundials/inst/include")
pevt=cppXPtr(code='
int event_ball(double t, const vec &y, vec &ynew, ivec &rootsfound, RObject &param, NumericVector &psens) {
  NumericVector p(param);
  static int nbounce=0;
  if (y[3] > 0) // we cross 0 in ascending trajectory, it can happen when y < 0 in limits of abstol
    return(RSUNDIALS_EVENT_IGNORE);
  ynew=y;
  if (++nbounce < p["nbounce"]) {
    // here nbounce=1:4
    ynew[2] *= 1.-p["kx"]; // horizontal speed is lowered
    ynew[3] *= -(1.-p["ky"]); // vertical speed is lowered and reflected
    return(RSUNDIALS_EVENT_HOLD);
  } else {
    // here nbounce=5
    nbounce=0; // reinit counter for possible next calls to cvode
    return(RSUNDIALS_EVENT_STOP);
  }
}
', depends=c("RcppArmadillo", "rmumps"), includes=paste0(includes, "\n#include <rsundials.h>"), cacheDir="lib", verbose=FALSE)
#outb <- rsundials::cvodes(yinib, timesb, pball, paramb, nroot=1, froot=proot, fevent=pevt)
#class(outb)=class(out); plot(outb)

rhs_ball_r=function(t, y, p, psens) {
  ydot=y # same length as y
  ydot[1] = y[3];
  ydot[2] = y[4];
  ydot[3] = if (y[2] > 0) 0. else -y[3] # falling till y=0 then damping
  ydot[4] = if (y[2] > 0) -p["g"] else -y[4] # falling till y=0 then damping
  return(ydot)
}
root_ball_r=function(t, y, p, psens) y[2]
event_ball_r=local({
  nbounce=0 # workaround for static variable
  function(t, y, rootsfound, p, psens) {
    if (y[4] > 0) # we cross 0 in ascending trajectory, it can happen when y < 0 in limits of abstol
      return(list(flag=RSUNDIALS_EVENT_IGNORE, ynew=y))
    nbounce <<- nbounce + 1
    ynew=y;
    if (nbounce < p["nbounce"]) {
      # here nbounce=1:4
      ynew[3] = ynew[3]*(1.-p["kx"]) # horizontal speed is lowered
      ynew[4] = -ynew[4]*(1.-p["ky"]) # vertical speed is lowered and reflected
      return(list(flag=RSUNDIALS_EVENT_HOLD, ynew=ynew))
    } else {
      # here nbounce=5
      nbounce <<- 0 # reinit counter for possible next calls to cvode
      return(list(flag=RSUNDIALS_EVENT_STOP, ynew=ynew))
    }
  }
})
#system.time(outbr <- rsundials::cvodes(yinib, timesb, rhs_ball_r, paramb, nroot=1, froot=root_ball_r, fevent=event_ball_r))
#class(outbr)=class(out); plot(outbr)

# decaying exp example
# y'=-nu*(y-ylim)
# pointer to rhs function
pexp=cppXPtr(code='
int d_exp(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
  ydot[0] = -psens["nu"]*(y[0]-psens["lim"]);
  return(0);
}
', depends="RcppArmadillo", includes=includes, cacheDir="lib", verbose=FALSE)
par_exp=c("nu"=1, "lim"=1)
#system.time(oute <- rsundials::cvodes(0., seq(0, 5, length.out=11), pexp, Ns=2, psens=par_exp))
#class(oute)=class(out); plot(oute)
