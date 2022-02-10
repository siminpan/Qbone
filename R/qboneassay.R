#' @include generics.R
#' @include utils.R
#' @importFrom methods setClass new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## The QboneAssay Class ----
#' The QboneAssay Class
#' The QboneAssay object is the basic unit of Qbone.


QboneAssay <- setClass(
  Class = 'QboneAssay',
  slots = c(
    data = 'list',
    scale.data = 'list',
    key = 'character',
    assay.orig = 'character',
    meta.assays = 'data.frame',
    misc = 'list'
  ),
  prototype = c(
    data = list(),
    scale.data = list(),
    key = character(),
    assay.orig = character(),
    meta.assays = data.frame(),
    misc = list()
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

