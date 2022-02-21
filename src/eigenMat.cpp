#include <RcppEigen.h>
#include <Rcpp.h>
#include "qbone_types.h"

using namespace Rcpp;
// using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// references  https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// penppml 	MIT + file LICENSE

//' Faster Matrix Multiplication
//'
//' Faster matrix multiplication using C++ Eigen.
//' A * B
//'
//' @param A,B Matrices.
//' @rdname eigenmm
//' @export
// [[Rcpp::export]]
SEXP eigenmm(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

//' Faster Matrix Transpose
//'
//' Faster matrix multiplication using C++ Eigen.
//' t(A)
//'
//' @param A Matrices.
//' @rdname eigenmt
// [[Rcpp::export]]
SEXP eigenmt(Eigen::MatrixXd A){

    return Rcpp::wrap(A.transpose());
}

//' Faster Matrix Transpose (Eigen::Map)
//'
//' Faster matrix multiplication using C++ Eigen::Map
//' t(A)
//'
//' @param A Matrices.
//' @rdname eigenmapmt
// [[Rcpp::export]]
SEXP eigenmapmt(Eigen::Map<Eigen::MatrixXd> A){

    return Rcpp::wrap(A.transpose());
}

//' Faster Matrix Multiplication (Eigen::Map)
//'
//' Faster matrix multiplication using C++ Eigen::Map.
//' A * B
//'
//' @param A,B Matrices.
//' @rdname eigenmapmm
//' @export
// [[Rcpp::export]]
SEXP eigenmapmm(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

//' Faster Matrix Multiplication (Eigen::Map)
//'
//' Faster matrix multiplication using C++ Eigen::Map.
//' A * t(B)
//'
//' @param A,B Matrices.
//' @rdname eigenmapmmt
// [[Rcpp::export]]
SEXP eigenmapmmt(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B.transpose();

    return Rcpp::wrap(C);
}

//' Faster Matrix Multiplication (Eigen::Map)
//'
//' Faster matrix multiplication using C++ Eigen::Map.
//' t(A) * B
//'
//' @param A,B Matrices.
//' @rdname eigenmapmtm
// [[Rcpp::export]]
SEXP eigenmapmtm(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A.transpose() * B;

    return Rcpp::wrap(C);
}
