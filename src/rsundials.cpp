// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <rsundials.h>
// [[Rcpp::export]]
NumericMatrix cvode(const NumericVector &yv, const NumericVector &times, const RObject &fderiv, RObject &param, const double abstol=1.e-8, const double reltol=1.e-8, RObject fjac=R_NilValue, RObject jacmat_=R_NilValue, std::string rmumps_perm="auto") {
  long clk_tck = CLOCKS_PER_SEC;
  clock_t t1, t2;
  t1 = clock();
  realtype t, tout;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  rsunDerivFn user_fn=NULL;
  IntegerVector vpermu=IntegerVector::create(_["amd"]=perm_amd, _["amf"]=perm_amf, _["scotch"]=perm_scotch, _["pord"]=perm_pord, _["metis"]=perm_metis, _["qamd"]=perm_qamd,  _["auto"]=perm_auto);
  int i;

  if (!fderiv.inherits("XPtr"))
    stop("fderiv must be of class 'XPtr'");
  if (!fjac.isNULL() && !fjac.inherits("XPtr"))
    stop("if not NULL, fjac must be of class 'XPtr'");
  if (vpermu.findName(rmumps_perm) == -1)
    stop("cvode: invalide rmumps_perm '%s'", rmumps_perm);

  List lp=List::create(_("param")=param, _("fderiv")=fderiv, _("fjac")=fjac);
//print(wrap(lp));
  void *cvode_mem, *user_data=static_cast<void*>(&lp);
  int iout, nnz;
  
  int neq=yv.size(), nti=times.size();
  NumericMatrix res(nti, neq+1);
  // copy times in the first column of result matrix
  res(_, 0)=times;
  
  /* Create serial vector of length neq from yv (initial conditions)*/
  getmem(y, N_VNew_Serial(neq));
  // copy init values
  std::copy(yv.begin(), yv.begin()+neq, NV_DATA_S(y));
  
  // create the solver memory and specify the Backward Differentiation Formula
  getmem(cvode_mem, CVodeCreate(CV_BDF));
  check_retval(CVodeSetErrHandlerFn(cvode_mem, rsunerr, NULL));

  // which jacobian: dense or sparse
  bool set_sparse;
  sp_mat asp;
  if (!jacmat_.isNULL() && jacmat_.inherits("simple_triplet_matrix")) {
    set_sparse=true;
    List jacmat(jacmat_);
    if (as<int>(jacmat["nrow"]) != neq)
      stop("nrow(jacmat) (%d) must be equal to nrow(yv) (%d)", as<int>(jacmat["nrow"]), neq);
    if (as<int>(jacmat["ncol"]) != neq)
      stop("ncol(jacmat) (%d) must be equal to nrow(yv) (%d)", as<int>(jacmat["ncol"]), neq);
  } else if (!jacmat_.isNULL() && jacmat_.inherits("dgCMatrix")) {
    set_sparse=true;
    S4 jacmat(jacmat_);
    IntegerVector d=jacmat.slot("Dim");
    if (d[0] != neq)
      stop("nrow(jacmat) (%d) must be equal to nrow(yv) (%d)", d[0], neq);
    if (d[1] != neq)
      stop("ncol(jacmat) (%d) must be equal to nrow(yv) (%d)", d[1], neq);
  } else {
    set_sparse=false;
  }
  if (set_sparse) {
    //nnz=as<NumericVector>(jacmat["v"]).size();
    asp=as<sp_mat>(jacmat_);
//asp.print("asp");
    nnz=asp.n_nonzero;
//Rcout << "nz=" << nnz << "\n";
    // check that main diagonal is in asp
    if (!all(vec(asp.diag())))
      stop("rsundials: all main diagonal elements must be present in sparse Jacobian");
    getmem(A, SUNSparseMatrix(neq, neq, nnz, CSC_MAT));
//printf("irp=%x; jcp=%x; vp=%x\n", SM_INDEXPTRS_S(A), SM_INDEXVALS_S(A), SM_DATA_S(A));
    // fill A by copying asp content
    std::copy(asp.col_ptrs, asp.col_ptrs+neq+1, SM_INDEXPTRS_S(A));
    std::copy(asp.row_indices, asp.row_indices+nnz, SM_INDEXVALS_S(A));
    //std::copy(asp.values, asp.values + nnz, SM_DATA_S(A));
//vec(SM_DATA_S(A), nnz, false).print("vasp");
    getmem(LS, SUNLinSol_RMUMPS(y, A, (rmumps_perm_t) as<int>(vpermu[rmumps_perm])));
  } else {
    // dense SUNMatrix for use in linear solves
    getmem(A, SUNDenseMatrix(neq, neq));
    // dense SUNLinearSolver object
    getmem(LS, SUNLinSol_Dense(y, A));
  }
  // Set cvode_mem and put different solver components
  check_retval(CVodeInit(cvode_mem, derivwrap, times[0], y));
  check_retval(CVodeSStolerances(cvode_mem, reltol, abstol));
  check_retval(CVodeSetUserData(cvode_mem, user_data));
  // attach the matrix and linear solver
  check_retval(CVodeSetLinearSolver(cvode_mem, LS, A));
  //check_retval(CVodeSetJacFn(cvode_mem, Jac));
//Rcout << "set_sparse=" << set_sparse << "\n";
  if (set_sparse) {
    check_retval(CVodeSetJacFn(cvode_mem, spjacwrap)); //jacsps)); //
  } else {
    check_retval(CVodeSetJacFn(cvode_mem, jacwrap)); // Jac)); //
  }

  // copy t0 y -> res(0,_)
  std::copy(NV_DATA_S(y), NV_DATA_S(y)+neq, res(0,_).begin()+1);
  // main loop
  for (int iout = 1; iout < nti; iout++) {
    check_retval(CVode(cvode_mem, times[iout], y, &t, CV_NORMAL));
    std::copy(NV_DATA_S(y), NV_DATA_S(y)+neq, res(iout,_).begin()+1);
//printf("ti=%g; tret=%g\n", times[iout], t);
//res.col(iout).print("yres");
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

/*
long int nje;
check_retval(CVodeGetNumJacEvals(cvode_mem, &nje));
Rcout << "n jac eval=" << nje << "\n";
*/
  // Free allocated memory
  N_VDestroy(y);
//Rcout << "cvode y mem is freed\n";
  CVodeFree(&cvode_mem);
//Rcout << "cvode_mem mem is freed\n";
  SUNLinSolFree(LS);
//Rcout << "cvode LS mem is freed\n";
  SUNMatDestroy(A);
//Rcout << "cvode A mem is freed\n";
  
  // set colnames in res
  StringVector colnm(neq);
  if (yv.hasAttribute("names"))
    colnm=yv.attr("names");
  else
    for (int i=0; i < neq; i++)
      colnm[i]="V"+std::to_string(i+1);
  colnm.insert(0, "time");
  res.attr("dimnames")=List::create(R_NilValue, colnm);
  res.attr("stats")=stats;
  t2 = clock();
  (void)printf("Temps consomme (s) : %lf \n",
                (double)(t2-t1)/(double)clk_tck);
  return(res);
}
//#define Ith(v,i)    NV_Ith_S(v,i-1) // Ith numbers components 1..NEQ
int derivwrap(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp[0];
  // get user function pointer
  XPtr< rsunDerivFn > xuser_fn=as<XPtr< rsunDerivFn >>(lp[1]);
  rsunDerivFn user_fn=*xuser_fn;
  vec yv(NV_DATA_S(y), NV_LENGTH_S(y), false),
      ydotv(NV_DATA_S(ydot), NV_LENGTH_S(ydot), false);
  return(user_fn(t, yv, ydotv, param));
}
int jacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp[0];
  // get user function pointer
  XPtr< rsunJacFn > xuser_fn=as<XPtr< rsunJacFn >>(lp[2]);
  rsunJacFn user_fn=*xuser_fn;
//printf("xuser_fn=%x; user_fn=%x\n", xuser_fn, user_fn);
  vec yv(NV_DATA_S(y), NV_LENGTH_S(y), false),
      ydotv(NV_DATA_S(ydot), NV_LENGTH_S(ydot), false);
  mat ja(SM_DATA_D(J), SM_ROWS_D(J), SM_COLUMNS_D(J), false);
  return(user_fn(t, yv, ydotv, ja, param));
}
int spjacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  List lp=*(static_cast<List *>(user_data));
  RObject param=lp[0];
  // get user function pointer
  XPtr< rsunSpJacFn > xuser_fn=as<XPtr< rsunSpJacFn >>(lp[2]);
  rsunSpJacFn user_fn=*xuser_fn;
//printf("xuser_fn=%x; user_fn=%x\n", xuser_fn, user_fn);
  vec yv(NV_DATA_S(y), NV_LENGTH_S(y), false),
      ydotv(NV_DATA_S(ydot), NV_LENGTH_S(ydot), false);
  int n=SM_COLUMNS_S(J), nz=SM_NNZ_S(J);
  uvec ir((unsigned int *) SM_INDEXVALS_S(J), nz, false);
  uvec p((unsigned int *) SM_INDEXPTRS_S(J), n+1, false);
  vec v((double *) SM_DATA_S(J), nz, false);
  return(user_fn(t, yv, ydotv, ir, p, v, n, nz, param));
}
/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */
#define Ith(v,i)    NV_Ith_S(v,i-1)         // Ith numbers components 1..NEQ
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) // IJth numbers rows,cols 1..NEQ
#define ZERO  RCONST(0.0)

int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y2, y3;
//printf("Jac t=%g\n", t);
  y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = RCONST(-0.04);
  IJth(J,1,2) = RCONST(1.0e4)*y3;
  IJth(J,1,3) = RCONST(1.0e4)*y2;

  IJth(J,2,1) = RCONST(0.04); 
  IJth(J,2,2) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
  IJth(J,2,3) = RCONST(-1.0e4)*y2;

  IJth(J,3,1) = ZERO;
  IJth(J,3,2) = RCONST(6.0e7)*y2;
  IJth(J,3,3) = ZERO;
  
/*mat jd(SUNDenseMatrix_Data(J), 3, 3, false);
jd.print("jd");*/

  return(0);
}

void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
  stop("cvode: '%s' failed: %s", function, msg);
}

#define ZERO  RCONST(0.0)
int jacsps(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  if (SUNMatGetID(J) != SUNMATRIX_SPARSE)
    stop("jacsps: wrong matrix type (expected sparse)");
  realtype *yval;
  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(J);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(J);
  realtype *data = SUNSparseMatrix_Data(J);
  
  yval = N_VGetArrayPointer(y);

  SUNMatZero(J);

  colptrs[0] = 0;
  colptrs[1] = 3;
  colptrs[2] = 6;
  colptrs[3] = 9;

  data[0] = RCONST(-0.04);
  rowvals[0] = 0;
  data[1] = RCONST(0.04);
  rowvals[1] = 1;
  data[2] = ZERO;
  rowvals[2] = 2;

  data[3] = RCONST(1.0e4)*yval[2];
  rowvals[3] = 0;
  data[4] = (RCONST(-1.0e4)*yval[2]) - (RCONST(6.0e7)*yval[1]);
  rowvals[4] = 1;
  data[5] = RCONST(6.0e7)*yval[1];
  rowvals[5] = 2;

  data[6] = RCONST(1.0e4)*yval[1];
  rowvals[6] = 0;
  data[7] = RCONST(-1.0e4)*yval[1];
  rowvals[7] = 1;
  data[8] = ZERO;
  rowvals[8] = 2;

Rcout << "t=" << t << "\n";
vec js(data, SM_NNZ_S(J), false);
js.print("js");

  return(0);
}
