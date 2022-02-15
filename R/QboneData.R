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
#' @slot key Key for the Assay                                                || double check
#' @slot assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @slot var.features Vector of features exhibiting high variance across      || double check
#' single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay   || double check
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
    assay.name = 'optionalCharacter',
    assay.orig = 'optionalCharacter',
    meta.assays = 'data.frame',
    misc = 'optionalList'
  ),
  prototype = c(
    data = list(),
    scale.data = list(),
    key = character(),
    assay.name = NULL,
    assay.orig = NULL,
    meta.assays = data.frame(),
    misc = NULL
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' create an QboneData object
#'
#' @param data data
#' @param meta.assays Include cells where at least this many features are
#' detected.
#' @param sampleid column number of sample name in meta.assay
#' @param assay.name assay name for the data.
#' @param assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @param ... || double check
#'
#' @return A \code{\link{QboneData}} object
#'
#' @importFrom methods as
#'
#' @export
#'
createQboneData <- function(
  data,
  meta.assays = NULL,
  sampleid = 1,
  assay.name = "Bone",
  assay.orig = NULL,
  ...
) {
  if (missing(x = data)) {
    stop("Must provide either 'data'")
  }
  if (anyDuplicated(x = meta.assays[,sampleid])) {
    warning(
      "Non-unique sample names (meta.assays[,sampleid]) present in the input meta.assays, making unique",
      call. = FALSE,
      immediate. = TRUE
    )
    meta.assays[,sampleid] <- make.unique(names = meta.assays[,sampleid])
  }
  if (is.list(data)){
    if (!is.null(x = meta.assays)) {
      if (nrow(x = meta.assays) != length(data)) {
        stop("There is a mismatch between the number of Metadata and the number of input data.")
      }
    } else {
      warning("meta.assays is not provided.")
    }
    for (i in 1:length(data)){
      if (is.atomic(data[[i]])){
        data[[i]] <- sort(data[[i]])
      } else {
        stop(paste0("Input data data[[",i,"]] is not atomic, please confirm data structure."))
      }
    }
  } else {
    stop("Input data expected to be a list. Other class is not supported yet.")
  }
  rownames(meta.assays) = meta.assays[,sampleid]
  names(data) <- rownames(meta.assays)
  qbonedata <- new(
    Class = 'QboneData',
    data = data,
    assay.name = assay.name,
    assay.orig = NULL,
    meta.assays = meta.assays,
    misc = NULL
  )
  return(qbonedata)
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
