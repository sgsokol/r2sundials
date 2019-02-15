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
# pointer to deriv function
pfnd=cppXPtr(code='
int d_robertson(double t, const arma::vec &y, arma::vec &ydot, Rcpp::RObject &param) {
  arma::vec p(as<Rcpp::NumericVector>(param));
//printf("d_robertson: t=%g\\n", t);
//Rcpp::print(param);
//y.print("y");
  ydot[0] = -p[0]*y[0] + p[2]*y[1]*y[2];
  ydot[2] = p[1]*y[1]*y[1];
  ydot[1] = -ydot[0] - ydot[2];
//ydot.print("ydot");
  return(0);
}
', depends="RcppArmadillo", cacheDir="lib", verbose=FALSE)
# pointer to dense jacobian function
pfnj=cppXPtr(code='
int jac_robertson(double t, const arma::vec &y, arma::vec &ydot, arma::mat &J, Rcpp::RObject &param) {
  arma::vec p(Rcpp::as<Rcpp::NumericVector>(param));
//printf("jac_robertson: t=%g\\n", t);
//Rcpp::print(param);
//y.print("y");
//ydot.print("ydot");
  J.at(0, 0) = -p[0];
  J.at(1, 0) = p[0]; 
  J.at(2, 0) = 0.;

  J.at(0, 1) = p[2]*y[2];
  J.at(1, 1) = -p[2]*y[2]-2*p[1]*y[1];
  J.at(2, 1) = 2*p[1]*y[1];

  J.at(0, 2) = p[2]*y[1];
  J.at(1, 2) = -p[2]*y[1];
  J.at(2, 2) = 0.;
//J.print("J");
  return(0);
}
', depends="RcppArmadillo", cacheDir="lib", verbose=FALSE)
# pointer to sparse jacobian function
pfnspj=cppXPtr(code='
int spjac_robertson(double t, const arma::vec &y, arma::vec &ydot, arma::uvec &ir, arma::uvec &pj, arma::vec &v, int n, int nz, Rcpp::RObject &param) {
  arma::vec prm(Rcpp::as<Rcpp::NumericVector>(param));
  int i=0;
  pj[0] = 0;
  // first column
  ir.at(i) = 0;
  v.at(i++) = -prm[0];
  ir.at(i) = 1;
  v.at(i++) = prm[0]; 
  pj[1] = i;
  // second column
  ir.at(i) = 0;
  v.at(i++) = prm[2]*y[2];
  ir.at(i) = 1;
  v.at(i++) = -prm[2]*y[2]-2*prm[1]*y[1];
  ir.at(i) = 2;
  v.at(i++) = 2*prm[1]*y[1];
  pj[2] = i;
  // third column
  ir.at(i) = 0;
  v.at(i++) = prm[2]*y[1];
  ir.at(i) = 1;
  v.at(i++) = -prm[2]*y[1];
  pj[3] = i;
  if (i != nz-1)
    stop("disconcordant i (%d) and nz-1 (%d)", i, nz-1);
//Rcout << "t=" << t << "\\n";
//v.print("jac v");
  return(0);
}
', depends="RcppArmadillo", cacheDir="lib", verbose=FALSE)

# sparse Jacobian
jproto=matrix(1., nrow=neq, ncol=neq)
jproto[3,1]=jproto[3,3]=0
jproto=jproto+diag(neq)
jproto=as.simple_triplet_matrix(jproto)
#print (system.time(
#out1 <- rsundials::cvode(yini, times, pfnd, parms, abstol=1.e-10, reltol=1.e-10, fjac=pfnspj, jacmat_=jproto)
#))
# dense Jacobian
jd=matrix(1., nrow=neq, ncol=neq)
#print (system.time(
#out2 <- rsundials::cvode(yini, times, pfnd, parms, abstol=1.e-10, reltol=1.e-10, fjac=pfnj, jacmat_=jd)
#))
