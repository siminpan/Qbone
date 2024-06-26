% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisfunction.R
\name{ecQuantlets}
\alias{ecQuantlets}
\title{Get the Compute Empirical Coefficients for quantlets}
\usage{
ecQuantlets(
  object,
  new.assay.name = "Empirical.Coefficients",
  data.assay = defaultAssay(object),
  k = NULL,
  sparsity = 0.001,
  ...
)
}
\arguments{
\item{object}{A Qboneobject}

\item{new.assay.name}{New assay name assigned to the quantlets data}

\item{data.assay}{It is the name of the assay whose data will be used
to compute the lasso list. Default is the data from the
defaultAssay(object).}

\item{k}{number of basis coefficients to keep. This is based on the number
of dictionary elements (C) to keep. Default will pick the one chosen from
\code{dxPlot()}.}

\item{sparsity}{Sparsity regularization parameter.}

\item{...}{Arguments passed to other methods}
}
\description{
Compute empirical coefficients for quantlets
}
