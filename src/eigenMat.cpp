#include <RcppEigen.h>

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// references  https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// penppml 	MIT + file LICENSE

// [[Rcpp::export]]
SEXP eigenmm(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenmt(Eigen::MatrixXd A){

    return Rcpp::wrap(A.transpose());
}

// [[Rcpp::export]]
SEXP eigenmapmt(Eigen::Map<Eigen::MatrixXd> A){

    return Rcpp::wrap(A.transpose());
}

// [[Rcpp::export]]
SEXP eigenmapmm(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenmapmmt(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B.transpose();

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenmapmtm(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A.transpose() * B;

    return Rcpp::wrap(C);
}
