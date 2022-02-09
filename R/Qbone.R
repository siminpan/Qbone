#' @include generics.R
#' @include assay.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Qbone Class
#'
#' The Qbone object is a representation of data analysis using Quantile Functional Regression using Quantlets (doi: 10.1080/01621459.2019.1609969) for R.
#'

Qbone <- setClass(
  Class = 'Qbone',
  slots = c(
    assays = 'list',
    # raw.data = 'list',
    # thin.data = 'list',
    meta.data = 'data.frame',
    # thin.meta = 'list',
    active.assay = 'character',
    active.ident = 'factor',
    # lasso.list = 'list',
    graphs = 'list',
    images = 'list',
    project.name = 'character',
    misc = 'list',
    # version = 'package_version',
    commands = 'list',
    tools = 'list'
  ),
  prototype = c(
    assays = list(),
    # raw.data = list(),
    # thin.data = list(),
    # meta.data = data.frame(id = NULL),
    # thin.meta = list(), # list(T, "0.1")
    active.assay = NA_character_,
    active.ident = factor(),
    # lasso.list = list(),
    graphs = list(),
    images = list(),
    project.name = "Qbone",
    misc = list(),
    # version = 'package_version',
    commands = list(),
    tools = list()
  )
)
