
R2cpp <- function(ode) {
  stopifnot("ode is not of type function"=is.function(ode))
  odefcpp <- ast2ast::translate(ode, output = "XPtr",
                                types_of_args = c("double", rep("sexp", 4)),
                                reference = TRUE,
                                return_type = "void",
                                getsource = TRUE)
  
  odefcpp <- strsplit(odefcpp, "\n")[[1]]
  odetmp <- character(length(odefcpp))
  counter <- 1
  for(i in odefcpp) {
    if("SEXP getXPtr() {" == i) break
    odetmp[counter] <- i
    counter <- counter + 1
  }
  odetmp <- paste(odetmp, collapse = "\n")
  
  fctptr <- '
  int wrapper_ode(double t, const arma::vec& y, arma::vec& ydot,
                Rcpp::RObject& param, Rcpp::NumericVector& psens) {
    sexp y_ = etr::vector(0.0, y.size());
    for(int i = 0; i < y.size(); i++) {
      y_[i] = y[i];
    }
    sexp ydot_(ydot.size(), ydot.memptr(), 2);
    Rcpp::NumericVector p(param);
    sexp param_(p.size(), p.begin(), 2);
    Rcpp::NumericVector ps(psens);
    sexp psens_(ps.size(), ps.begin(), 2);
    ode(t, y_, ydot_, param_, psens_);
    return(0);
  }
  SEXP getXPtr() {                                                                                                   
  typedef int (*fct_ptr) (
                double t, const arma::vec& y, arma::vec& ydot,
                Rcpp::RObject& param, Rcpp::NumericVector& psens);                        
  return Rcpp::XPtr<fct_ptr>(new fct_ptr(&  wrapper_ode ));                                                                  
  }    
  '

  ode <- paste(odetmp, fctptr)

  env <- new.env()
  Rcpp::sourceCpp(code = ode, verbose = TRUE, env = env) 
  fct_ret <- env$getXPtr()
  attributes(fct_ret) <- list(class = "XPtr")

  return(fct_ret)
}




