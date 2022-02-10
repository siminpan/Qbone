#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Samples ----
# Samples <- function(x) {
#   UseMethod(generic = 'Samples', object = x)
# }

setGeneric("samples", function(x) standardGeneric("samples"))

## CreateQboneObject ----
# CreateQboneObject <- function(
#   counts,
#   project = 'CreateQboneObject',
#   assay = 'Bone',
#   names.field = 1,
#   names.delim = '_',
#   meta.data = NULL,
#   ...
# ) {
#   UseMethod(generic = 'CreateQboneObject', object = counts)
# }

setGeneric("createQboneObject", function(x) standardGeneric("createQboneObject"))

## DefaultAssay ----
# DefaultAssay <- function(object, ...) {
#   UseMethod(generic = 'defaultAssay', object = object)
# }

setGeneric("defaultAssay", function(x) standardGeneric("defaultAssay"))


## "DefaultAssay<-" ----
# "DefaultAssay<-" <- function(object, ..., value) {
#   UseMethod(generic = 'DefaultAssay<-', object = object)
# }

setGeneric("defaultAssay<-", function(x) standardGeneric("defaultAssay<-"))

## AddMetaData ----
# AddMetaData <- function(object, metadata, col.name = NULL) {
#   UseMethod(generic = 'AddMetaData', object = object)
# }

setGeneric("addMetaData", function(x) standardGeneric("addMetaData"))


## Idents ----
# Idents <- function(object, ... ) {
#   UseMethod(generic = 'Idents', object = object)
# }

setGeneric("idents", function(x) standardGeneric("idents"))


## "Idents<-" ----
# "Idents<-" <- function(object, ..., value) {
#   UseMethod(generic = 'Idents<-', object = object)
# }

setGeneric("idents<-", function(x) standardGeneric("idents<-"))
