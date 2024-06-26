% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionalregression.R
\name{qfrModel}
\alias{qfrModel}
\title{Fit the quantile functional regression model}
\usage{
qfrModel(
  object,
  new.assay.name = "Q.F.Regression",
  data.assay = defaultAssay(object),
  X = NULL,
  X1 = NULL,
  delta2 = 0.95,
  H = NULL,
  pct.range = c(0.05, 0.95),
  assay.seed2 = .Random.seed,
  ...
)
}
\arguments{
\item{object}{A Qboneobject}

\item{new.assay.name}{New assay name assigned to the quantlets data}

\item{data.assay}{It is the name of the assay whose data will be used
to compute the lasso list. Default is the data from the
defaultAssay(object).}

\item{X}{Covariates (N by A matrix, N is the number of observations)}

\item{X1}{New covariates for inferencing MCMC fit result.}

\item{delta2}{Cutoff percentage of the energy for the functional coefficients}

\item{H}{Number of clustering groups. Cluster basis indices based on their eigen-values.}

\item{pct.range}{Percentage range from original data set. Use as the range of the domain for the density estimate. Default is c(0.05, 0.95). You can override it by assign \code{xranges = c(min, max)}}

\item{...}{Arguments passed to other methods}
}
\description{
Fit the quantlet space model in using Markov chain Monte Carlo (MCMC)
}
