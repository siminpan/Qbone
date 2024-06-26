% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisfunction.R
\name{lassoList}
\alias{lassoList}
\title{Use penalized regression (lasso) to find a sparse subset of dictionary
elements}
\usage{
lassoList(
  object,
  verbose = TRUE,
  new.assay.name = "Lasso.list",
  data.assay = defaultAssay(object),
  alpha = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  beta = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  assay.seed = .Random.seed,
  parallel = T,
  ...
)
}
\arguments{
\item{object}{A Qboneobject}

\item{verbose}{Print a progress bar}

\item{new.assay.name}{New assay name assigned to the lassolist data}

\item{data.assay}{It is the name of the assay whose data will be used to
compute the lasso list. Default is the data from the defaultAssay(object).}

\item{alpha}{Vector containing sequence of beta parameter for internal
function \code{generateBetaCDF()}}

\item{beta}{Vector containing sequence of beta parameter for internal
function \code{generateBetaCDF()}}

\item{assay.seed}{assay information to add into the QboneData object
scale.data. The default of \code{lassolist()} will save the random seed for
the run.
Use \code{.Random.seed <-  object@assays[["Lasso.list"]]@scale.data[["lassolist"]] }
before run \code{lassolist()} for the same results.}

\item{parallel}{If TRUE, use parallel foreach to fit each fold in
\code{cv.glmnet()}. Default use \code{registerDoMC()} to register parallel.
There is another function \code{lassolist_parallel()} for overall parallel
computing for all sample.}

\item{...}{Arguments passed to other methods}
}
\description{
First construct overcomplete dictionary (Beta CDF). Then uses penalized
regression (lasso) to find a sparse subset of dictionary elements.
}
