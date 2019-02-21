// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <rsundials.h>
// [[Rcpp::export]]
NumericMatrix cvode(const NumericVector &yv, const vec &times, const RObject &frhs, RObject param=R_NilValue, const double abstol=1.e-8, const double reltol=1.e-8, const std::string integrator="bdf", const int maxord=0, const int maxsteps=0, const vec &constraints=NumericVector::create(), const RObject fjac=R_NilValue, const int nz=0, const std::string rmumps_perm="auto", const int nroot=0, const RObject froot=R_NilValue, const RObject fevent=R_NilValue) {
  //long clk_tck = CLOCKS_PER_SEC;
  //clock_t t1, t2;
  //t1 = clock();
  realtype t, tout;
  N_Vector y, nv_constraints;
  SUNMatrix A;
  SUNLinearSolver LS;
  IntegerVector vpermu=IntegerVector::create(_["amd"]=perm_amd, _["amf"]=perm_amf, _["scotch"]=perm_scotch, _["pord"]=perm_pord, _["metis"]=perm_metis, _["qamd"]=perm_qamd,  _["auto"]=perm_auto);
  int retval, retevent, i;
  ivec rootsfound;
  mat mroots;
  Function rf_event(Environment::global_env()["ls"]); // just a placeholder

  if (frhs.sexp_type() != CLOSXP && !frhs.inherits("XPtr"))
    stop("frhs must be either R function or of class 'XPtr'");
  if (!fjac.isNULL() && fjac.sexp_type() != CLOSXP && !fjac.inherits("XPtr"))
    stop("if not NULL, fjac must be either R function or of class 'XPtr'");
  if (vpermu.findName(rmumps_perm) == -1)
    stop("cvode: invalide rmumps_perm '%s'", rmumps_perm);
  if (integrator != "adams" && integrator != "bdf")
    stop("cvode: invalide integrator '%s' (must be either 'adams' or 'bdf')", integrator);

  List lp=List::create(_("param")=param, _("frhs")=frhs, _("fjac")=fjac, _("froot")=froot, _("fevent")=fevent, _("nroot")=nroot);
//print(wrap(lp));
  void *cvode_mem, *user_data=static_cast<void*>(&lp);
  int iout, nnz;
  
  int neq=yv.size(), nti=times.size();
  mat res(neq+1, nti); // will be transposed at return, can be resized if roots are added and/or event handling causes premature stop
  // copy times in the first row of result matrix
  res.row(0)=times.t();
  
  /* Create serial vector of length neq from yv (initial conditions)*/
  getmem(y, N_VNew_Serial(neq));
  vec yvec(NV_DATA_S(y), neq, false), ynew; // vec proxy for y and y after event treatment by user function
  // copy init values
  yvec.subvec(0, neq-1)=vec(yv);
  
  // create the solver memory and specify the Backward Differentiation Formula
  getmem(cvode_mem, CVodeCreate(integrator == "bdf" ? CV_BDF : CV_ADAMS));
  check_retval(CVodeSetErrHandlerFn(cvode_mem, rsunerr, NULL));
  // Set cvode_mem and put different solver components
  if (maxord > 0)
    check_retval(CVodeSetMaxOrd(cvode_mem, maxsteps));
  check_retval(CVodeSetMaxNumSteps(cvode_mem, maxord));
  check_retval(CVodeInit(cvode_mem, rhswrap, times[0], y));
  check_retval(CVodeSStolerances(cvode_mem, reltol, abstol));
  check_retval(CVodeSetUserData(cvode_mem, user_data));
  if (constraints.size() != 0 && constraints.size() != neq)
    stop("cvode: if not empty, constraints vector must be of the same size as yv (%d), actually its size is %d", neq, constraints.size());
  if (any(constraints)) {
    // check that only 0, ±1 and ±2 are present
    for (double it: constraints) {
      double ita=fabs(it);
      if (ita != 0. && ita != 1. && ita != 2.)
        stop("cvode: values in constraints must be 0, ±1 or ±2, instead %g found", it);
    }
    getmem(nv_constraints, N_VNewEmpty_Serial(neq));
    NV_DATA_S(nv_constraints) = (realtype *) constraints.begin();
    check_retval(CVodeSetConstraints(cvode_mem, nv_constraints));
  }

  // which jacobian: dense or sparse
  if (nz > 0) {
    if (fjac.isNULL())
      stop("cvode: fjac must not be NULL if nz > 0 (%d)", nz);
    getmem(A, SUNSparseMatrix(neq, neq, nz, CSC_MAT));
    getmem(LS, SUNLinSol_RMUMPS(y, A, (rmumps_perm_t) as<int>(vpermu[rmumps_perm])));
    check_retval(CVodeSetLinearSolver(cvode_mem, LS, A));
    check_retval(CVodeSetJacFn(cvode_mem, spjacwrap)); //jacsps)); //
  } else {
    // dense SUNMatrix for use in linear solves
    getmem(A, SUNDenseMatrix(neq, neq));
    // dense SUNLinearSolver object
    getmem(LS, SUNLinSol_Dense(y, A));
    check_retval(CVodeSetLinearSolver(cvode_mem, LS, A));
    if (!fjac.isNULL())
      check_retval(CVodeSetJacFn(cvode_mem, jacwrap)); // Jac)); //
  }
  // setup root treatment
  if (nroot > 0) {
    if (froot.isNULL())
      stop("rsundials: froot cannot be NULL when nroot > 0 (%d)", nroot);
    if (froot.sexp_type() != CLOSXP && froot.sexp_type() != EXTPTRSXP)
      stop("rsundials: froot must be an R function or an XPtr returned by RcppXPtrUtils::cppXPtr()");
    if (fevent.isNULL())
      stop("rsundials: fevent cannot be NULL when nroot > 0 (%d)", nroot);
    if (fevent.sexp_type() != CLOSXP && fevent.sexp_type() != EXTPTRSXP)
      stop("rsundials: fevent must be an R function or an XPtr returned by RcppXPtrUtils::cppXPtr()");
    rootsfound.resize(nroot);
    ynew.set_size(neq);
    mroots.set_size(nroot+1, 0);
    check_retval(CVodeRootInit(cvode_mem, nroot, rootwrap));
    if (fevent.sexp_type() == CLOSXP)
      rf_event=as<Function>(fevent);
  }

  // copy t0 y -> res(_, 0)
  res.col(0).subvec(1, neq)=yvec;
  // main loop
  for (int iout=1, insroot=0; iout < nti; iout++) {
    retval=CVode(cvode_mem, times[iout], y, &t, CV_NORMAL);
    if (retval == CV_ROOT_RETURN) {
      // root found proceed it
      check_retval(CVodeGetRootInfo(cvode_mem, rootsfound.begin()));
      if (fevent.sexp_type() == CLOSXP) {
        // call R fevent
//Rcout << "rf_event=" << rf_event << "\n";
//print(fevent);
        List lres=rf_event(t, yvec, rootsfound, param);
//Rcout << "lres=\n";
//print(lres);
        ynew.subvec(0, neq-1)=as<vec>(lres["ynew"]);
        retevent=as<int>(lres["flag"]);
      } else {
        // call XPtr fevent
        XPtr< rsunEventFn > xuser_fn=as<XPtr< rsunEventFn >>(fevent);
        rsunEventFn user_fn=*xuser_fn;
        retevent=user_fn(t, yvec, ynew, rootsfound, param);
      }
      // treat retevent
      if (retevent == RSUNDIALS_EVENT_IGNORE) {
        ; // do nothing special
      } else {
        mroots.insert_cols(mroots.n_cols, 1, false);
        mroots.col(mroots.n_cols-1)[0]=t;
        mroots.col(mroots.n_cols-1).subvec(1, nroot)=conv_to<vec>::from(rootsfound);
        // insert this time point and current state
        if (t < times[iout]) {
          res.insert_cols(iout+insroot, 1, false);
          res.col(iout+insroot)[0]=t;
          res.col(iout+insroot).subvec(1, neq)=yvec;
          insroot++;
        }
        if (any(ynew != yvec)) {
          // insert new state too
          res.insert_cols(iout+insroot, 1, false);
          res.col(iout+insroot)[0]=t;
          res.col(iout+insroot).subvec(1, neq)=ynew;
          insroot++;
          // reinit cvode
          yvec.subvec(0, neq-1)=ynew;
          check_retval(CVodeReInit(cvode_mem, t, y));
        }
        if (retevent == RSUNDIALS_EVENT_STOP) {
          // free extra cols in res
          res.resize(neq+1, iout+insroot);
          break;
        } else if (retevent != RSUNDIALS_EVENT_HOLD) {
          stop("cvode: fevent() returned unknown flag %d", retevent);
        }
      }
      if (t < times[iout])
        iout--; // stay at the same time point for the next iteration
    } else {
      // simple storage of current solution
      res.col(iout+insroot).subvec(1, neq)=yvec;
    }
    checkUserInterrupt();
  }
  // final stats as attribute to output
  IntegerVector stats(9);
  stats.attr("names")=CharacterVector::create("NumSteps", "NumRhsEvals", "NumLinSolvSetups", "NumErrTestFails", "NumNonlinSolvIters", "NumNonlinSolvConvFails", "NumJacEvals", "NumLinRhsEvals", "NumGEvals");
  i=0;
  check_retval(CVodeGetNumSteps(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumRhsEvals(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumLinSolvSetups(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumErrTestFails(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumNonlinSolvIters(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumNonlinSolvConvFails(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumJacEvals(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumLinRhsEvals(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumGEvals(cvode_mem, (long int*) &stats[i++]));

  // Free allocated memory
  N_VDestroy(y);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  
  // set colnames in resout
  StringVector colnm(neq);
  if (yv.hasAttribute("names"))
    colnm=yv.attr("names");
  else
    for (int i=0; i < neq; i++)
      colnm[i]="V"+std::to_string(i+1);
  colnm.insert(0, "time");
  NumericMatrix resout=wrap(res.t());
  resout.attr("dimnames")=List::create(R_NilValue, colnm);
  resout.attr("stats")=stats;
  if (nroot > 0)
    resout.attr("roots")=wrap(mroots.t());
  //t2 = clock();
  //(void)printf("Temps consomme (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);
  return(resout);
}

int rhswrap(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false), ydotv(NV_DATA_S(ydot), neq, false);
  if (TYPEOF(lp["frhs"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["frhs"]);
    ydotv.subvec(0, neq-1)=as<vec>(user_fr(wrap(t), wrap(yv), param));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunRhsFn > xuser_fn=as<XPtr< rsunRhsFn >>(lp["frhs"]);
    rsunRhsFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, param));
  }
}
int jacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false),
      ydotv(NV_DATA_S(ydot), neq, false);
  mat ja(SM_DATA_D(J), neq, neq, false);
  if (TYPEOF(lp["fjac"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fjac"]);
    ja.submat(0, neq-1, 0, neq-1)=as<mat>(user_fr(wrap(t), wrap(yv), wrap(ydotv), param));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunJacFn > xuser_fn=as<XPtr< rsunJacFn >>(lp["fjac"]);
    rsunJacFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, ja, param));
  }
}
int spjacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y), nz=SM_NNZ_S(J);
  vec yv(NV_DATA_S(y), neq, false), ydotv(NV_DATA_S(ydot), neq, false);
  uvec ir((unsigned int *) SM_INDEXVALS_S(J), nz, false);
  uvec p((unsigned int *) SM_INDEXPTRS_S(J), neq+1, false);
  vec v((double *) SM_DATA_S(J), nz, false);
  if (TYPEOF(lp["fjac"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fjac"]);
    List res=user_fr(wrap(t), wrap(yv), wrap(ydotv), List::create(_["i"]=ir, _["p"]=p, _["v"]=v), param);
    uvec rir=res["i"];
    uvec rp=res["p"];
    vec rv=res["v"];
    ir.subvec(0, rir.n_elem-1)=rir;
    p.subvec(0, neq)=rp;
    v.subvec(0, rv.n_elem-1)=rv;
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunSpJacFn > xuser_fn=as<XPtr< rsunSpJacFn >>(lp["fjac"]);
    rsunSpJacFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, ir, p, v, neq, nz, param));
  }
}
int rootwrap(realtype t, N_Vector y, realtype *rootout, void *user_data) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  int nroot=lp["nroot"];
  vec yv(NV_DATA_S(y), neq, false);
  vec vroot((double *) rootout, nroot, false);
  if (TYPEOF(lp["froot"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["froot"]);
    vroot.subvec(0, nroot-1)=as<vec>(user_fr(wrap(t), wrap(yv), param));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunRootFn > xuser_fn=as<XPtr< rsunRootFn >>(lp["froot"]);
    rsunRootFn user_fn=*xuser_fn;
    return(user_fn(t, yv, vroot, param));
  }
}
void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
  stop("%s: %s", function, msg);
}
