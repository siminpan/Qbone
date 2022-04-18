#' @include zzz.R
#' @include generics.R
#' @include utils.R
#' @importFrom methods setClass new slot slot<- slotNames
#' @importFrom stats na.omit
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Class definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 The QboneData Class ----
#' The QboneData Class
#' The QboneData object is the basic unit of Qbone.
#'
#' @slot data Sample data
#' @slot scale.data Parameter of data processing
#' @slot assay.name name of assay
#' @slot assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @slot meta.assays Metadata for the data processing
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
    assay.name = 'optionalCharacter',
    assay.orig = 'optionalCharacter',
    meta.assays = 'data.frame'
  ),
  prototype = c(
    data = list(),
    scale.data = list(),
    assay.name = NULL,
    assay.orig = NULL,
    meta.assays = data.frame()
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 create an QboneData object ----
#' create an QboneData object
#'
#' @param data Sample data
#' @param meta.assays metadata for the assay with sample names
#' @param sampleid.assays column number of sample name in meta.assay
#' @param assay.name assay name for the data.
#' @param assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @param sort sort the data before put into data slot if T, default is F.
#' @param ...
#'
#' @return A \code{\link{QboneData}} object
#'
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' \dontrun{
#' n = 10000
#' list1 = list(rnorm(n, mean = 0, sd = 1),
#'              rnorm(n, mean = 0, sd = 1),
#'             rnorm(n, mean = 0.5, sd = 1),
#'             rnorm(n, mean = 0.5, sd = 1))
#' meta.data = data.frame(name = c("a_1", "a_2", "b_1", "b_2"), group = c("a", "a", "b", "b"))
#' rownames(meta.data) = meta.data[,1]
#' qbonedata = createQboneData(list1, meta.data, sampleid = 1)
#' }
#'
createQboneData <- function(
  data,
  meta.assays = NULL,
  sampleid.assays = 1,
  assay.name = "Bone",
  assay.orig = NULL,
  sort = F,
  ...
){
  if (missing(x = data)){
    stop("Must provide either 'data'")
  }
  if (anyDuplicated(x = meta.assays[,sampleid.assays])){
    warning(
      "Non-unique sample names (meta.assays[,sampleid.assays]) present in the input meta.assays, making unique",
      call. = FALSE,
      immediate. = TRUE
    )
    meta.assays[,sampleid.assays] <- make.unique(names = meta.assays[,sampleid.assays])
  }
  if (is.list(data)){
    if (!is.null(x = meta.assays)){
      if (nrow(x = meta.assays) != length(data)){
        stop("There is a mismatch between the number of Metadata and the number of input data.")
      }
    } else {
      warning("meta.assays is not provided.")
    }
    for (i in 1:length(data)){
      if (is.atomic(data[[i]])){
        if (sort == T){
          data[[i]] <- sort(data[[i]], method = "quick")
        } else {
          data[[i]] <- data[[i]]
          }
      } else {
        stop(paste0("Input data data[[",i,"]] is not atomic, please confirm data structure."))
      }
    }
  } else {
    stop("Input data expected to be a list. Other class is not supported yet.")
  }
  rownames(meta.assays) = meta.assays[,sampleid.assays]
  names(data) <- rownames(meta.assays)
  qbonedata <- new(
    Class = 'QboneData',
    data = data,
    assay.name = assay.name,
    assay.orig = assay.orig,
    meta.assays = data.frame(id = meta.assays[,sampleid.assays],
                             row.names = meta.assays[,sampleid.assays])
  )
  # qbonedata@scale.data <- append(qbonedata@scale.data,
  #                                list(id = meta.assays[,sampleid.assays]))
  if (sort == T){
    qbonedata@scale.data <- append(qbonedata@scale.data, list(sort = T))
  } else {
    qbonedata@scale.data <- append(qbonedata@scale.data, list(sort = F))
  }
  return(qbonedata)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Qbone-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## 3.1 getQboneData.QboneData ----
#' @rdname assayData
#' @export
#' @method getQboneData QboneData
#'
#' @examples
#' # Get the data directly from an QboneData object
#' getQboneData(qbone1[["Bone"]], slot = "data")
#'
getQboneData.QboneData <- function(
  object,
  slot = c('data', 'scale.data'), # , 'scale.data', 'counts'            || Double check
  ...
){
  CheckDots(...)
  slot <- slot[1]
  slot <- match.arg(arg = slot)
  return(slot(object = object, name = slot))
}

## 3.4 samples.QboneData ----
#' @rdname samples
#' @export
#'
samples.QboneData <- function(x){
  return(names(x@data))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. R-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 QboneData Methods ----
#'  QboneData Methods
#' \code{QboneData} Methods
#'
#' Methods for \code{\link{QboneData}} objects for generics defined in
#' other packages
#'
#' @param x,object An \code{\link{QboneData}} object
#' @param i,features For \code{[[}: metadata names; for all other methods,
#' feature names or indices
#' @param j,samples sample names or indices
#' @param ... Arguments passed to other methods
#'
#' @name QboneData-methods
#' @rdname QboneData-methods
#'
#' @concept QboneData
#'
NULL

## 4.1.1 [[.QboneData ----
#' @describeIn QboneData-methods Metadata and associated object accessor
#'
#' @param drop See \code{\link[base]{drop}}
#'
#' @return \code{[[}: If \code{i} is missing, the metadata data frame; if
#' \code{i} is a vector of metadata names, a data frame with the requested
#' metadata, otherwise, the requested associated object
#'
#' @export
#' @method [[ QboneData
#'
"[[.QboneData" <- function(x, i, ..., drop = FALSE){
  if (missing(x = i)){
    i <- colnames(x = slot(object = x, name = 'meta.assay'))
  }
  data.return <- slot(object = x, name = 'meta.assay')[, i, drop = FALSE, ...]
  if (drop){
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
  }
  return(data.return)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. S4 methods ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Internal ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
