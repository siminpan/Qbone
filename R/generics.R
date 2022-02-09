#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Samples <- function(x) {
  UseMethod(generic = 'Samples', object = x)
}

CreateQboneObject <- function(
  counts,
  project = 'CreateQboneObject',
  assay = 'Bone',
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  ...
) {
  UseMethod(generic = 'CreateQboneObject', object = counts)
}

DefaultAssay <- function(object, ...) {
  UseMethod(generic = 'DefaultAssay', object = object)
}

"DefaultAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultAssay<-', object = object)
}

AddMetaData <- function(object, metadata, col.name = NULL) {
  UseMethod(generic = 'AddMetaData', object = object)
}

Idents <- function(object, ... ) {
  UseMethod(generic = 'Idents', object = object)
}

"Idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'Idents<-', object = object)
}
