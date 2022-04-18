#' @include zzz.R
#' @importFrom methods setGeneric
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 samples ----
#' Get samples present in an object
#' @param x An object
#'
#' @return A vector of sample names
#'
#' @rdname samples
#' @export samples
#'
#' @examples
#' samples(x = qbone1)
#'
#'
samples <- function(x) UseMethod(generic = 'samples', object = x)
# as S3 method
# samples <- function(x) UseMethod(generic = 'samples', object = x)
# as S4 method
# setGeneric("samples", function(x) standardGeneric("samples"))

## 1.2 createQboneObject ----

#' Create a \code{Qbone} object
#'
#' @inheritParams createQboneData
#'
#' @param data data
#' @param project \link{project} name for the \code{Qbone} object
#' @param assay Name of the initial assay
#' @param names.field number of field of names
#' @param names.delim delimiter of names
#' @param meta.data Include cells where at least this many features are
#' detected.
#' @param sampleid column number of sample name in mate.data
#'
#' @rdname createQboneObject
#' @export
#'
#' @examples
#' \dontrun{
#' n = 10000
#' list1 = list(rnorm(n, mean = 0, sd = 1),
#'              rnorm(n, mean = 0, sd = 1),
#'              rnorm(n, mean = 0.5, sd = 1),
#'              rnorm(n, mean = 0.5, sd = 1))
#' meta.data = data.frame(name = c("a_1", "a_2", "b_1", "b_2"),
#'                        group = c("a", "a", "b", "b"))
#' rownames(meta.data) = meta.data[,1]
#' qbonedata = createQboneData(list1, meta.data, sampleid = 1)
#' qbone1 = createQboneObject(qbonedata, meta.data = meta.data)
#' }
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
){
  UseMethod(generic = 'createQboneObject', object = data)
}

# setGeneric("createQboneObject", function(x) standardGeneric("createQboneObject"))

## 1.3 defaultAssay ----

#' Get and set the default assay
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{DefaultAssay}: The name of the default assay
#'
#' @rdname defaultAssay
#' @export defaultAssay
#'
defaultAssay <- function(object, ...){
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
"defaultAssay<-" <- function(object, ..., value){
  UseMethod(generic = 'defaultAssay<-', object = object)
}

# setGeneric("defaultAssay<-", function(x) standardGeneric("defaultAssay<-"))

## 1.4 addMetaData ----
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
#' @examples
#' meta.data2 = data.frame(group2 = c("a2", "a2", "b2", "b2"))
#' qbone1 <- addMetaData(object = qbone1, metadata = meta.data2[,1], col.name = "group2")
#' # or
#' qbone1[["group2"]] <- meta.data2[,1]
#'
addMetaData <- function(object, metadata, col.name = NULL){
  UseMethod(generic = 'addMetaData', object = object)
}

# setGeneric("addMetaData", function(x) standardGeneric("addMetaData"))


## 1.5 idents ----
#' Get, set, and manipulate an object's identity classes
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname idents
#' @export idents
#'
#' @examples
#' # Get sample identity classes
#' idents(qbone1)
#'
idents <- function(object, ... ){
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
#' @examples
#' # Set sample identity classes
#' # Can be used to set identities for specific samples to a new level
#' idents(qbone1, samples = 3:4) <- 'c'
#' idents(qbone1)
#'
#' # Can also set idents from a value in object metadata
#' colnames(qbone1[[]])
#' idents(qbone1) <- colnames(qbone1[[]])[3]
#' # or
#' idents(qbone1) <- 'group2'
#' idents(qbone1)
#'
"idents<-" <- function(object, ..., value){
  UseMethod(generic = 'idents<-', object = object)
}

# setGeneric("idents<-", function(x) standardGeneric("idents<-"))

## 1.6 getQboneData ----
#' Get and Set Assay Data
#'
#' General accessor and setter functions for \code{\link{QboneData}} objects.
#' \code{getQboneData} can be used to pull information from any of the
#' expression matrices (eg. \dQuote{counts}, \dQuote{data}, or
#' \dQuote{scale.data}). \code{SetAssayData} can be used to replace one of these
#' expression matrices                                              || double check
#'
#' @param object An object
#' @param slot Specific assay data to get or set
#' @param ... Arguments passed to other methods
#'
#' @return \code{getQboneData}: returns the specified assay data
#'
#' @name assayData
#' @rdname assayData
#' @export getQboneData
#'
#' @order 1
#'
#' @concept data-access
#'
getQboneData <- function(object, slot, ...){
  UseMethod(generic = 'getQboneData', object = object)
}
