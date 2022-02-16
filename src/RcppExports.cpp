// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// eigenMatMult
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _Qbone_eigenMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMatMulttrans
SEXP eigenMatMulttrans(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _Qbone_eigenMatMulttrans(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMulttrans(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMapMatMult
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _Qbone_eigenMapMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMapMatMulttrans
SEXP eigenMapMatMulttrans(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _Qbone_eigenMapMatMulttrans(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapMatMulttrans(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigentransMapMatMult
SEXP eigentransMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _Qbone_eigentransMapMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigentransMapMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Qbone_eigenMatMult", (DL_FUNC) &_Qbone_eigenMatMult, 2},
    {"_Qbone_eigenMatMulttrans", (DL_FUNC) &_Qbone_eigenMatMulttrans, 2},
    {"_Qbone_eigenMapMatMult", (DL_FUNC) &_Qbone_eigenMapMatMult, 2},
    {"_Qbone_eigenMapMatMulttrans", (DL_FUNC) &_Qbone_eigenMapMatMulttrans, 2},
    {"_Qbone_eigentransMapMatMult", (DL_FUNC) &_Qbone_eigentransMapMatMult, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Qbone(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}