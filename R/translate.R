R2cpp <- function(ode) {
  stopifnot("ode is not of type function" = is.function(ode))
  odefcpp <- ast2ast::translate(ode,
    references = rep(TRUE, 5),
    types_of_args = c("double", "const double", rep("double", 3)),
    data_structures = c("scalar", rep("borrow", 4)),
    handle_inputs = rep("", 5),
    output = "XPtr", getsource = TRUE
  )

  odefcpp <- strsplit(odefcpp, "\n")[[1]]
  odetmp <- character(length(odefcpp))
  counter <- 1
  for (i in odefcpp) {
    if ("SEXP getXPtr ()  {" == i) break
    odetmp[counter] <- i
    counter <- counter + 1
  }
  odetmp <- paste(odetmp, collapse = "\n")

  fctptr <- "
  int wrapper_ode(double t, const arma::vec& y, arma::vec& ydot,
                Rcpp::RObject& param, Rcpp::NumericVector& psens) {
    etr::Vec<const double, etr::Borrow<const double>> y_(y.memptr(), y.size());
    etr::Vec<double, etr::Borrow<double>> ydot_(ydot.memptr(), ydot.size());
    Rcpp::NumericVector p(param);
    etr::Vec<double, etr::Borrow<double>> param_(p.begin(),p.size());
    Rcpp::NumericVector ps(psens);
    etr::Vec<double, etr::Borrow<double>> psens_(ps.begin(),ps.size());
    ode(t, y_, ydot_, param_, psens_);
    return(0);
  }
  SEXP getXPtr() {
  typedef int (*fct_ptr) (
                double t, const arma::vec& y, arma::vec& ydot,
                Rcpp::RObject& param, Rcpp::NumericVector& psens);
  return Rcpp::XPtr<fct_ptr>(new fct_ptr(&  wrapper_ode ));
  }
  "
  ode <- paste(odetmp, fctptr)
  env <- new.env()
  Rcpp::sourceCpp(code = ode, verbose = TRUE, env = env)
  fct_ret <- env$getXPtr()
  attributes(fct_ret) <- list(class = "XPtr")
  return(fct_ret)
}
