#' @include utils.R
#' @importFrom methods setOldClass setClassUnion slot slot<-
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(name = 'optionalCharacter', members = c('NULL', 'character'))
setClassUnion(name = 'optionalList', members = c('NULL', 'list'))

setOldClass(Classes = 'package_version')
