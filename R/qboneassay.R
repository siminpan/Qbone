#' @include zzz.R
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
    assay.orig = 'optionalCharacter',
    meta.assays = 'data.frame',
    misc = 'optionalList'
  ),
  prototype = c(
    data = list(),
    scale.data = list(),
    key = character(),
    assay.orig = NULL,
    meta.assays = data.frame(),
    misc = NULL
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

