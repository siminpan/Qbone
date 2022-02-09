#' @include generics.R
#' @include qboneassay.R
#' @importFrom methods setClass
#' @importFrom data.table fread
#' @useDynLib Qbone
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Qbone Class
#'
#' The Qbone object is a representation of data analysis using Quantile
#' Functional Regression using Quantlets (doi: 10.1080/01621459.2019.1609969)
#' for R.
#'

Qbone <- setClass(
  Class = 'Qbone',
  slots = c(
    assays = 'list',
    meta.data = 'data.frame',
    active.assay = 'character',
    active.ident = 'factor',
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
    active.assay = NA_character_,
    active.ident = factor(),
    graphs = list(),
    images = list(),
    project.name = "Qbone",
    misc = list(),
    # version = 'package_version',
    commands = list(),
    tools = list()
  )
)
