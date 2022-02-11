#' @include zzz.R
#' @include generics.R
#' @include utils.R
#' @importFrom methods setClass new slot slot<- slotNames
#' @importFrom stats na.omit
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## The QboneData Class ----
#' The QboneData Class
#' The QboneData object is the basic unit of Qbone.
#'
#' @slot data Normalized expression data
#' @slot scale.data Scaled expression data
#' @slot key Key for the Assay
#' @slot assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @slot var.features Vector of features exhibiting high variance across
#' single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay
#'
#' @name QboneData-class
#' @rdname QboneData-class
#' @exportClass QboneData
#'
#'
#' @seealso \code{\link{QboneData-methods}}
#'
QboneData <- setClass(
  Class = 'QboneData',
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

#' Create an QboneData object
#'
#' @param data Prenormalized data; if provided, do not pass \code{counts}
#' @param min.cells Include features detected in at least this many cells. Will
#' subset the counts matrix as well. To reintroduce excluded features, create a
#' new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#' detected.
#' @param check.matrix Check counts matrix for NA, NaN, Inf, and non-integer values
#' @param ... Arguments passed to \code{\link{as.sparse}}
#'
#' @return A \code{\link{Assay}} object
#'
#' @importFrom methods as
#'
#' @export
#'
CreateAssayObject <- function(
  data,
  min.cells = 0,
  min.features = 0,
  check.matrix = FALSE,
  ...
) {

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Qbone-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# R-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
