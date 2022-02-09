#' @include generics.R
#' @importFrom methods new slot slot<-
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
