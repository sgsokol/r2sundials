// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <rsundials.h>
// [[Rcpp::export]]
NumericMatrix cvodes(const NumericVector &yv, const vec &times, const RObject &frhs, RObject param=R_NilValue, const double abstol=1.e-8, const double reltol=1.e-8, const Nullable<int> integrator_=R_NilValue, const int maxord=0, const int maxsteps=0, const vec &constraints=NumericVector::create(), const RObject fjac=R_NilValue, const int nz=0, const Nullable<int> rmumps_perm_=R_NilValue, const int nroot=0, const RObject froot=R_NilValue, const RObject fevent=R_NilValue,
const int Ns=0, NumericVector psens=NumericVector::create(), NumericVector psens_bar=NumericVector::create(), const IntegerVector psens_list=IntegerVector::create(), const RObject fsens=R_NilValue, const RObject fsens1=R_NilValue, const Nullable<int> sens_method_=R_NilValue, const bool errconS=true) {
  //long clk_tck = CLOCKS_PER_SEC;
  //clock_t t1, t2;
  //t1 = clock();
  UserData udata;
  realtype t, tout;
  N_Vector nv_y, nv_constraints, *yS;
  SUNMatrix A;
  SUNLinearSolver LS;
  int retval, retevent, i;
  ivec rootsfound, plist;
  mat mroots;
  Function rf_event(Environment::global_env()["ls"]); // just a placeholder
  rsunEventFn user_event_fn;
  const int integrator=(integrator_.isNotNull() ? as<int>(integrator_) : CV_BDF),
    rmumps_perm=(rmumps_perm_.isNotNull() ? as<int>(rmumps_perm_) : RMUMPS_PERM_AUTO),
    sens_method=(sens_method_.isNotNull() ? as<int>(sens_method_) : CV_SIMULTANEOUS);

  if (frhs.sexp_type() != CLOSXP && !frhs.inherits("XPtr"))
    stop("rsundials: frhs must be either R function or of class 'XPtr'");
  if (!fjac.isNULL() && fjac.sexp_type() != CLOSXP && !fjac.inherits("XPtr"))
    stop("rsundials: if not NULL, fjac must be either R function or of class 'XPtr'");
  if (!fevent.isNULL() && fevent.sexp_type() != CLOSXP && !fevent.inherits("XPtr"))
    stop("rsundials: if not NULL, fevent must be either R function or of class 'XPtr'");
  if (!fevent.isNULL() && fevent.inherits("XPtr")) {
    user_event_fn=*(as<XPtr< rsunEventFn >>(fevent));
  }

  if (integrator != CV_ADAMS && integrator != CV_BDF)
    stop("cvode: invalide integrator_ '%d' (must be either NULL, %d (CV_ADAMS) or %d (CV_BDF default))", integrator, CV_ADAMS, CV_BDF);

  udata.lp=List::create(_("param")=param, _("frhs")=frhs, _("fjac")=fjac, _("froot")=froot, _("fevent")=fevent, _("nroot")=nroot, _("fsens")=fsens, _("fsens1")=fsens1);
  udata.psens=psens;
//print(wrap(udata.lp));
  void *cvode_mem, *user_data=static_cast<void*>(&udata);
  int iout, nnz;
  
  int neq=yv.size(), nti=times.size();
  mat res(neq+1, nti); // will be transposed at return, can be resized if roots are added and/or event handling causes premature stop
  // copy times in the first row of result matrix
  res.row(0)=times.t();
  cube asens(neq, Ns, nti);
  std::vector<vec> ySv(Ns);
  
  /* Create serial vector of length neq from yv (initial conditions)*/
  getmem(nv_y, N_VNew_Serial(neq));
  vec yvec(NV_DATA_S(nv_y), neq, false), ynew; // vec proxy for y and y after event treatment by user function
  // copy init values
  yvec.subvec(0, neq-1)=vec(yv);
  
  // create the solver memory and specify the Backward Differentiation Formula
  getmem(cvode_mem, CVodeCreate(integrator));
  check_retval(CVodeSetErrHandlerFn(cvode_mem, rsunerr, NULL));
  // Set cvode_mem and put different solver components
  if (maxord > 0)
    check_retval(CVodeSetMaxOrd(cvode_mem, maxsteps));
  check_retval(CVodeSetMaxNumSteps(cvode_mem, maxord));
  check_retval(CVodeInit(cvode_mem, rhswrap, times[0], nv_y));
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
    getmem(LS, SUNLinSol_RMUMPS(nv_y, A, rmumps_perm));
    check_retval(CVodeSetLinearSolver(cvode_mem, LS, A));
    check_retval(CVodeSetJacFn(cvode_mem, spjacwrap)); //jacsps)); //
  } else {
    // dense SUNMatrix for use in linear solves
    getmem(A, SUNDenseMatrix(neq, neq));
    // dense SUNLinearSolver object
    getmem(LS, SUNLinSol_Dense(nv_y, A));
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
  // setup sensitivity
  bool sens_user=!fsens.isNULL() || !fsens1.isNULL();
  if (Ns > 0) {
    if (psens.size() != Ns) {
      stop("rsundials: if Ns > 0 (%d) then length(ps) (%d) must be equal to Ns", Ns, psens.size());
    }
    if(!fsens.isNULL() && !fsens1.isNULL()) {
      stop("rsundials: fsens and fsens1 cannot both be to not NULL. At least one of them or both must NULL");
    }
    if (!sens_user) {
      // prepare param vectors
      if (psens.size() == 0) {
        stop("rsundials: psens must not be empty if no fsens neither fsens1 are set while Ns > 0");
      }
      if (psens.size() < Ns) {
        stop("rsundials: if not empty, psens must not be shorter than Ns (%d < %d)", psens.size(), Ns);
      }
      if (psens.size() > Ns && psens_list.size() == 0)
        stop("rsundials: if length(psens) > Ns, psens_list must not be empty");
      if (psens_list.size() == 0) {
        plist=regspace<ivec>(0, Ns-1);
      } else {
        if (psens_list.size() != Ns)
          stop("rsundials: if psens_list is not empty, its size (%s) must be equal to Ns (%d)", psens_list.size(), Ns);
        if (psens_list.size() > psens.size())
          stop("rsundials: length(psens) must not be shorter than length(psens_list) (%d < %d)", psens.size(), psens_list.size());
        ivec u=unique(psens_list);
        if (u.size() < psens_list.size())
          stop("rsundials: psens_list has non unique entries");
        if (min(psens_list) < 1)
          stop("rsundials: minimal value in psens_list is less than 1 (%d)", min(psens_list));
        if (max(psens_list) > psens.size())
          stop("rsundials: maximal value in psens_list is greater than length(psens) (%d > %d)", max(psens_list), psens.size());
        plist=as<ivec>(psens_list);
        plist -= 1;
      }
    }
    // prepare yS=0
    yS = N_VCloneVectorArray_Serial(Ns, nv_y);
    for (int is=0; is < Ns; is++) {
      N_VConst(0., yS[is]);
      ySv[is]=vec(NV_DATA_S(yS[is]), neq, false); // vec proxies of yS "matrix"
    }
    // set sens problem
    if (!fsens.isNULL()) {
      check_retval(CVodeSensInit(cvode_mem, Ns, sens_method, sensrhswrap, yS));
    } else if (!fsens1.isNULL()) {
      check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method, sensrhs1wrap, yS));
    } else {
      check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method, NULL, yS));
    }
    check_retval(CVodeSetSensParams(cvode_mem, psens.begin(), psens_bar.size() == 0 ? NULL : psens_bar.begin(), plist.begin()));
    check_retval(CVodeSensEEtolerances(cvode_mem));
    check_retval(CVodeSetSensErrCon(cvode_mem, errconS));
  }

  // time course
  // copy t0 nv_y -> res(_, 0)
  res.col(0).subvec(1, neq)=yvec;
  // main loop
  for (int iout=1, insroot=0; iout < nti; iout++) {
    retval=CVode(cvode_mem, times[iout], nv_y, &t, CV_NORMAL);
    if (retval == CV_ROOT_RETURN) {
      // root found proceed it
      check_retval(CVodeGetRootInfo(cvode_mem, rootsfound.begin()));
      if (fevent.sexp_type() == CLOSXP) {
        // call R fevent
//Rcout << "rf_event=" << rf_event << "\n";
//print(fevent);
        List lres=rf_event(t, yvec, rootsfound, param, psens);
//Rcout << "lres=\n";
//print(lres);
        ynew.subvec(0, neq-1)=as<vec>(lres["ynew"]);
        retevent=as<int>(lres["flag"]);
      } else {
        // call XPtr fevent
        retevent=user_event_fn(t, yvec, ynew, rootsfound, param, psens);
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
          if (Ns > 0) {
            CVodeGetSens(cvode_mem, &t, yS);
            asens.insert_slices(iout+insroot, 1, false);
            for (int is=0; is < Ns; is++)
              asens.slice(iout+insroot).col(is)=ySv[is];
          }
          insroot++;
        }
        if (any(ynew != yvec)) {
          // insert new state too
          res.insert_cols(iout+insroot, 1, false);
          res.col(iout+insroot)[0]=t;
          res.col(iout+insroot).subvec(1, neq)=ynew;
          if (Ns > 0) {
            // reinit yS to 0
            asens.insert_slices(iout+insroot, 1, false);
            for (int is=0; is < Ns; is++) {
              N_VConst(0., yS[is]); // reinit sensitivities
              asens.slice(iout+insroot).col(is)=ySv[is];
            }
          }
          insroot++;
          // reinit cvode
          yvec.subvec(0, neq-1)=ynew;
          check_retval(CVodeReInit(cvode_mem, t, nv_y));
          if (Ns > 0) {
            // reinit sensitivities
            CVodeSensFree(cvode_mem);
            if (!fsens.isNULL()) {
              check_retval(CVodeSensInit(cvode_mem, Ns, sens_method, sensrhswrap, yS));
            } else if (!fsens1.isNULL()) {
              check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method, sensrhs1wrap, yS));
            } else {
              check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method, NULL, yS));
            }
          }
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
      if (Ns > 0) {
        check_retval(CVodeGetSens(cvode_mem, &t, yS));
        for (int is=0; is < Ns; is++)
          asens.slice(iout+insroot).col(is)=ySv[is];
      }
    }
    checkUserInterrupt();
  }
  // final stats as attribute to output
  IntegerVector stats(15);
  stats.attr("names")=CharacterVector::create("NumSteps", "NumRhsEvals", "NumLinSolvSetups", "NumErrTestFails", "NumNonlinSolvIters", "NumNonlinSolvConvFails", "NumJacEvals", "NumLinRhsEvals", "NumGEvals",
    "SensNumRhsEvals", "NumRhsEvalsSens", "SensNumLinSolvSetups", "SensNumErrTestFails", "SensNumNonlinSolvIters", "SensNumNonlinSolvConvFails"
  );
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

  // stats and Free allocated memory
  if (Ns > 0) {
    check_retval(CVodeGetSensNumRhsEvals(cvode_mem, (long int*) &stats[i++]));
    check_retval(CVodeGetNumRhsEvalsSens(cvode_mem, (long int*) &stats[i++]));
    check_retval(CVodeGetSensNumLinSolvSetups(cvode_mem, (long int*) &stats[i++]));
    if (errconS) {
      check_retval(CVodeGetSensNumErrTestFails(cvode_mem, (long int*) &stats[i++]));
    }
    if (sens_method == CV_STAGGERED) {
      check_retval(CVodeGetSensNumNonlinSolvIters(cvode_mem, (long int*) &stats[i++]));
      check_retval(CVodeGetSensNumNonlinSolvConvFails(cvode_mem, (long int*) &stats[i++]));
    }
    CVodeSensFree(cvode_mem);
    N_VDestroyVectorArray(yS, Ns);
  }
  N_VDestroy(nv_y);
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
  if (Ns > 0)
    resout.attr("sens")=wrap(asens);
  //t2 = clock();
  //(void)printf("Temps consomme (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);
  return(resout);
}

int rhswrap(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false), ydotv(NV_DATA_S(ydot), neq, false);
  if (TYPEOF(lp["frhs"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["frhs"]);
    ydotv.subvec(0, neq-1)=as<vec>(user_fr(wrap(t), wrap(yv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunRhsFn > xuser_fn=as<XPtr< rsunRhsFn >>(lp["frhs"]);
    rsunRhsFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, param, udata->psens));
  }
}
int jacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false),
    ydotv(NV_DATA_S(ydot), neq, false),
    tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false),
    tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false),
    tmp3v(NV_DATA_S(tmp3), NV_LENGTH_S(tmp3), false);
  mat ja(SM_DATA_D(J), neq, neq, false);
  if (TYPEOF(lp["fjac"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fjac"]);
    ja.submat(0, neq-1, 0, neq-1)=as<mat>(user_fr(wrap(t), wrap(yv), wrap(ydotv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunJacFn > xuser_fn=as<XPtr< rsunJacFn >>(lp["fjac"]);
    rsunJacFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, ja, param, udata->psens, tmp1v, tmp2v, tmp3v));
  }
}
int spjacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y), nz=SM_NNZ_S(J);
  vec yv(NV_DATA_S(y), neq, false), ydotv(NV_DATA_S(ydot), neq, false),
    tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false),
    tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false),
    tmp3v(NV_DATA_S(tmp3), NV_LENGTH_S(tmp3), false);
  uvec ir((unsigned int *) SM_INDEXVALS_S(J), nz, false);
  uvec p((unsigned int *) SM_INDEXPTRS_S(J), neq+1, false);
  vec v((double *) SM_DATA_S(J), nz, false);
  if (TYPEOF(lp["fjac"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fjac"]);
    List res=user_fr(wrap(t), wrap(yv), wrap(ydotv), List::create(_["i"]=ir, _["p"]=p, _["v"]=v), param, udata->psens, tmp1v, tmp2v, tmp3v);
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
    return(user_fn(t, yv, ydotv, ir, p, v, neq, nz, param, udata->psens, tmp1v, tmp2v, tmp3v));
  }
}
int rootwrap(realtype t, N_Vector y, realtype *rootout, void *user_data) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  int nroot=lp["nroot"];
  vec yv(NV_DATA_S(y), neq, false);
  vec vroot((double *) rootout, nroot, false);
  if (TYPEOF(lp["froot"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["froot"]);
    vroot.subvec(0, nroot-1)=as<vec>(user_fr(wrap(t), wrap(yv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunRootFn > xuser_fn=as<XPtr< rsunRootFn >>(lp["froot"]);
    rsunRootFn user_fn=*xuser_fn;
    return(user_fn(t, yv, vroot, param, udata->psens));
  }
}
int sensrhswrap(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false);
  vec ydotv(NV_DATA_S(ydot), neq, false);
  vec tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false);
  vec tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false);
  std::vector<vec> ySv(Ns), ySdotv(Ns);
  // init stdvec
  for (int i=0; i < Ns; i++) {
    ySv[i]=vec(NV_DATA_S(yS[i]), neq, false);
    ySdotv[i]=vec(NV_DATA_S(ySdot[i]), neq, false);
  }
  if (TYPEOF(lp["fsens"]) == CLOSXP) {
    // plain R function
    mat ySm(neq, Ns);
    for (int is=0; is < Ns; is++)
      ySm.col(is)=ySv[is];
    Function user_fr=as<Function>(lp["fsens"]);
    DataFrame res=user_fr(wrap(Ns), wrap(t), wrap(yv), wrap(ydotv), wrap(ySm), param, udata->psens);
    for (int i=0; i < Ns; i++)
      ySdotv[i].subvec(0, neq-1)=as<vec>(res[i]);
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunSensFn > xuser_fn=as<XPtr< rsunSensFn >>(lp["fsens"]);
    rsunSensFn user_fn=*xuser_fn;
    return(user_fn(Ns, t, yv, ydotv, ySv, ySdotv, param, udata->psens, tmp1v, tmp2v));
  }
}
int sensrhs1wrap(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false);
  vec ydotv(NV_DATA_S(ydot), neq, false);
  vec ySv(NV_DATA_S(yS), neq, false);
  vec ySdotv(NV_DATA_S(ySdot), neq, false);
  vec tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false);
  vec tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false);
  
  if (TYPEOF(lp["fsens1"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fsens1"]);
    ySv.subvec(0, neq-1)=as<vec>(user_fr(wrap(Ns), wrap(t), wrap(yv), wrap(iS), wrap(ydotv), wrap(ySv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunSens1Fn > xuser_fn=as<XPtr< rsunSens1Fn >>(lp["fsens1"]);
    rsunSens1Fn user_fn=*xuser_fn;
    return(user_fn(Ns, t, yv, ydotv, iS, ySv, ySdotv, param, udata->psens, tmp1v, tmp2v));
  }
}

// helpers
void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
  stop("%s: %s", function, msg);
}

#define cnst_pair(c) {#c, c}

// [[Rcpp::export]]
int get_cnst(std::string s) {
  static std::map<std::string, int> dict={
    cnst_pair(CV_SUCCESS),
    
    cnst_pair(CV_BDF),
    cnst_pair(CV_ADAMS),
    
    cnst_pair(RSUNDIALS_EVENT_IGNORE),
    cnst_pair(RSUNDIALS_EVENT_HOLD),
    cnst_pair(RSUNDIALS_EVENT_STOP),
    
    cnst_pair(RMUMPS_PERM_AMD),
    cnst_pair(RMUMPS_PERM_AMF),
    cnst_pair(RMUMPS_PERM_SCOTCH),
    cnst_pair(RMUMPS_PERM_PORD),
    cnst_pair(RMUMPS_PERM_METIS),
    cnst_pair(RMUMPS_PERM_QAMD),
    cnst_pair(RMUMPS_PERM_AUTO),
    
    cnst_pair(CV_SIMULTANEOUS),
    cnst_pair(CV_STAGGERED),
    cnst_pair(CV_STAGGERED1)
  };
  if (dict.count(s))
    return(dict[s]);
  else
    stop("get_cnst: constant '%s' is not in dictionary", s);
  return(NA_INTEGER);
}