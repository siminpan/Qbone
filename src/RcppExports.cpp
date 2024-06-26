// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// eigenmm
SEXP eigenmm(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _Qbone_eigenmm(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenmm(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenmt
SEXP eigenmt(Eigen::MatrixXd A);
RcppExport SEXP _Qbone_eigenmt(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(eigenmt(A));
    return rcpp_result_gen;
END_RCPP
}
// eigenmapmt
SEXP eigenmapmt(Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _Qbone_eigenmapmt(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(eigenmapmt(A));
    return rcpp_result_gen;
END_RCPP
}
// eigenmapmm
SEXP eigenmapmm(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _Qbone_eigenmapmm(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenmapmm(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenmapmmt
SEXP eigenmapmmt(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _Qbone_eigenmapmmt(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenmapmmt(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenmapmtm
SEXP eigenmapmtm(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _Qbone_eigenmapmtm(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenmapmtm(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Qbone_eigenmm", (DL_FUNC) &_Qbone_eigenmm, 2},
    {"_Qbone_eigenmt", (DL_FUNC) &_Qbone_eigenmt, 1},
    {"_Qbone_eigenmapmt", (DL_FUNC) &_Qbone_eigenmapmt, 1},
    {"_Qbone_eigenmapmm", (DL_FUNC) &_Qbone_eigenmapmm, 2},
    {"_Qbone_eigenmapmmt", (DL_FUNC) &_Qbone_eigenmapmmt, 2},
    {"_Qbone_eigenmapmtm", (DL_FUNC) &_Qbone_eigenmapmtm, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Qbone(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
