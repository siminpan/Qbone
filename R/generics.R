#' @include zzz.R
#' @importFrom methods setGeneric
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Samples ----
#' Get samples present in an object
#' @param x An object
#'
#' @return A vector of sample names
#'
#' @rdname samples
#' @export samples
#'
samples <- function(x) UseMethod(generic = 'samples', object = x)
# as S3 method
# samples <- function(x) UseMethod(generic = 'samples', object = x)
# as S4 method
# setGeneric("samples", function(x) standardGeneric("samples"))

## createQboneObject ----

#' Create a \code{Qbone} object
#'
#' @inheritParams createQboneData
#'
#' @param data data
#' @param project \link{Project} name for the \code{Qbone} object
#' @param assay Name of the initial assay
#' @param names.field || double check
#' @param names.delim || double check
#' @param meta.data Include cells where at least this many features are
#' detected.
#' @param sampleid Column number of sample name in meta.data
#'
#' @rdname createQboneObject
#' @export
#'
createQboneObject <- function(
  data,
  project = 'createQboneObject',
  assay = 'Bone',
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  sampleid = 1,
  ...
) {
  UseMethod(generic = 'createQboneObject', object = data)
}

# setGeneric("createQboneObject", function(x) standardGeneric("createQboneObject"))

## DefaultAssay ----

#' Get and set the default assay
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{DefaultAssay}: The name of the default assay
#'
#' @rdname defaultAssay
#' @export defaultAssay
#'
defaultAssay <- function(object, ...) {
  UseMethod(generic = 'defaultAssay', object = object)
}

# setGeneric("defaultAssay", function(x) standardGeneric("defaultAssay"))

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value Name of assay to set as default
#'
#' @return \code{defaultAssay<-}: An object with the default assay updated
#'
#' @rdname defaultAssay
#' @export defaultAssay<-
#'
"defaultAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'sefaultAssay<-', object = object)
}

# setGeneric("defaultAssay<-", function(x) standardGeneric("defaultAssay<-"))

## addMetaData ----
#' Add in metadata associated
#' @param object An object
#' @param metadata A vector, list, or data.frame with metadata to add
#' @param col.name A name for meta data if not a named list or data.frame
#'
#' @return \code{object} with metadata added
#'
#' @rdname addMetaData
#' @export addMetaData
#'
#'
addMetaData <- function(object, metadata, col.name = NULL) {
  UseMethod(generic = 'addMetaData', object = object)
}

# setGeneric("addMetaData", function(x) standardGeneric("addMetaData"))


## Idents ----
#' Get, set, and manipulate an object's identity classes
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname idents
#' @export idents
#'
idents <- function(object, ... ) {
  UseMethod(generic = 'idents', object = object)
}

# setGeneric("idents", function(x) standardGeneric("idents"))

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value The name of the identities to pull from object metadata or the
#' identities themselves
#'
#' @return \code{idents<-}: \code{object} with the cell identities changed
#'
#' @rdname idents
#' @export idents<-
#'
"idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'idents<-', object = object)
}

# setGeneric("idents<-", function(x) standardGeneric("idents<-"))
