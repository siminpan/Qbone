#' @include zzz.R
#' @include generics.R
#' @include utils.R
#' @include QboneData.R
#' @importFrom methods setClass new slot slot<- slotNames setMethod
#' @importFrom stats na.omit
#' @importFrom utils argsAnywhere isS3method isS3stdGeneric methods
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Class definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 The Qbone Class ----
#' The Qbone Class
#'
#' The Qbone object is a representation of data analysis using Quantile
#' Functional Regression using Quantlets (doi: 10.1080/01621459.2019.1609969)
#' for R.
#'
#' @slot assays A list of QboneData Object for this project.
#' @slot meta.data Meta information regarding each sample.
#' @slot active.assay Name of the active, or default, assay; settable using
#' \code{\link{defaultAssay}}
#' @slot active.ident The active identity for this Qbone object;
#' settable using \code{\link{idents}}
#' @slot project.name Name of the project
#' @slot version Version of Qbone this object was built under
#'
#' @name Qbone-class
#' @rdname Qbone-class
#' @exportClass Qbone
#'
Qbone <- setClass(
  Class = 'Qbone',
  slots = c(
    assays = 'list',
    meta.data = 'data.frame',
    active.assay = 'character',
    active.ident = 'factor',
    project.name = 'character',
    version = 'package_version' # function OK but infinite list under this slot.
    # list(packageVersion()) ended up the same.
  ),
  prototype = c(
    assays = list(),
    active.assay = NA_character_,
    active.ident = factor(),
    project.name = "Qbone"
    # version = packageVersion(pkg = "Qbone"),
  )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1.1 createQboneObject.default ----
#'
#' @rdname createQboneObject
#' @method createQboneObject default
#' @export
#'
#' @examples
#' \dontrun{
#' n = 10000
#' list1 = list(rnorm(n, mean = 0, sd = 1),
#'              rnorm(n, mean = 0, sd = 1),
#'              rnorm(n, mean = 0.5, sd = 1),
#'              rnorm(n, mean = 0.5, sd = 1))
#' meta.data = data.frame(name = c("a_1", "a_2", "b_1", "b_2"), group = c("a", "a", "b", "b"))
#' rownames(meta.data) = meta.data[,1]
#' qbonedata = createQboneData(list1, meta.data, sampleid = 1)
#' qbone1 = createQboneObject(qbonedata, meta.data = meta.data)
#' }
#'
createQboneObject.default <- function(
  data = data,
  project = 'QboneProject',
  assay = "Bone",
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  sampleid = 1,
  ...
){
  if (!is.null(x = meta.data)){
    if (nrow(x = meta.data) != length(data)){
      stop("There is a mismatch between the number of Metadata and the number of input data.")
    }
  }
  data.Qbone <- createQboneData(
    data,
    meta.assays = meta.data,
    sampleid.assays = sampleid,
    assay.name = assay,
    assay.orig = NULL,
    sort = T
    )
    return(
      createQboneObject(
        data = data.Qbone,
        project = 'QboneProject',
        assay = assay,
        names.field = names.field,
        names.delim = names.delim,
        meta.data = meta.data,
        sampleid = sampleid,
        ...
        )
      )
}

## 2.1.2 createQboneObject.QboneData ----
#' @rdname createQboneObject
#' @method createQboneObject QboneData
#' @export
#'
createQboneObject.QboneData <- function(
  data = data,
  project = 'QboneProject',
  assay = "Bone",
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  sampleid = 1,
  ...
){
  if (!is.null(x = meta.data)){
    # if (is.null(x = rownames(x = meta.data))){
    #   stop("Row names not set in metadata. Please ensure that rownames of metadata match sample names of data")
    # }
    # if (length(x = setdiff(x = rownames(x = meta.data), y = rownames(data@meta.assays)))){
    #   warning("Some samples in meta.data not present in provided data.")
    #   meta.data <- meta.data[intersect(x = rownames(x = meta.data), y = rownames(data@meta.assays)), , drop = FALSE]
    # }
    if (length(x = setdiff(x = meta.data[,sampleid], y = rownames(data@meta.assays
)))){
      warning("Some samples in meta.data not present in provided data.")
      meta.data <- meta.data[intersect(x = meta.data[,sampleid], y = rownames(data@meta.assays
      )), , drop = FALSE]
    }
    if (is.data.frame(x = meta.data)){
      new.meta.data <- data.frame(row.names = rownames(data@meta.assays
      ))
      for (ii in 1:ncol(x = meta.data)){
        new.meta.data[meta.data[,sampleid], colnames(x = meta.data)[ii]] <- meta.data[, ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  } else {
    meta.data <- data@meta.assays
  }
  # samples <- list(data@data)
  # names(x = samples) <- "samples"
  # assay.list <- list(samples)
  assay.list <- list(data)
  names(x = assay.list) <- assay
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    # X = colnames(x = meta.data),
    X = names(data@data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  if (any(is.na(x = idents))){
    warning(
      "Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)){
    idents <- rep.int(x = factor(x = project), times = length(data@data))
  }
  names(x = idents) <- names(data@data)

  object <- new(
    Class = 'Qbone',
    assays = assay.list,
    meta.data =  data.frame(row.names = names(data@data)), # meta.data,
    active.assay = assay,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'Qbone'),
    ...
  )
  if (!is.null(x = meta.data)){
    object <- addMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

## 2.2 assays ----
#' Query Specific Object Types
#'
#' Adopted from Seurateobject package
#' List the names of \code{\link{assay}}
#'
#' @param object A \code{\link{Qbone}} object
#' @param slot Name of component object to return
#'
#' @return If \code{slot} is \code{NULL}, the names of all component objects
#' in this \code{Qbone} object. Otherwise, the specific object specified
#'
#' @rdname ObjectAccess
#'
#' @export
#'
#' @concept data-access
#'
assays <- function(object, slot = NULL){
  assays <- FilterObjects(object = object, classes.keep = 'QboneData')
  if (is.null(x = slot)){
    return(assays)
  }
  if (!slot %in% assays){
    warning(
      "Cannot find an assay of name ",
      slot,
      " in this Seurat object",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(slot(object = object, name = 'assays')[[slot]])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Qbone-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 addMetaData ----
#' @rdname addMetaData
#' @export
#' @method addMetaData Qbone
#'
addMetaData.Qbone <- .addMetaData

## 3.2.1 defaultAssay.Qbone ----
#' @rdname defaultAssay
#' @export
#' @method defaultAssay Qbone
#'
defaultAssay.Qbone <- function(object, ...){
  CheckDots(...)
  object <- updateSlots(object = object)
  return(slot(object = object, name = 'active.assay'))
}

## 3.2.2 defaultAssay<-.Qbone" ----
#' @rdname defaultAssay
#' @export
#' @method defaultAssay<- Qbone
#'
"defaultAssay<-.Qbone" <- function(object, ..., value){
  CheckDots(...)
  object <- updateSlots(object = object)
  if (!value %in% names(x = slot(object = object, name = 'assays'))){
    stop("Cannot find assay ", value)
  }
  slot(object = object, name = 'active.assay') <- value
  return(object)
}

## 3.3.1 idents.Qbone ----
#' @param object An Qbone object
#' @rdname idents
#' @export
#' @method idents Qbone
#'
idents.Qbone <- function(object, ...){
  CheckDots(...)
  object <- updateSlots(object = object)
  return(slot(object = object, name = 'active.ident'))
}


## 3.3.2 idents<-.Qbone ----
#' @param object An Qbone object
#' @param samples Set cell identities for specific samples
#' @param drop Drop unused levels
#'
#' @rdname idents
#' @export
#' @method idents<- Qbone
#'
"idents<-.Qbone" <- function(object, samples = NULL, drop = FALSE, ..., value){
  CheckDots(...)
  object <- updateSlots(object = object)
  samples <- samples %||% samples(object)
  if (is.numeric(x = samples)){
    samples <- samples(object)[samples]
  }
  samples <- intersect(x = samples, y = samples(object))
  samples <- match(x = samples, table = samples(object))
  if (length(x = samples) == 0){
    warning("Cannot find samples provided")
    return(object)
  }
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[[]])){
    unlist(x = object[[value]], use.names = FALSE)[samples]
  } else {
    if (is.list(x = value)){
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = samples))
  }
  new.levels <- if (is.factor(x = idents.new)){
    levels(x = idents.new)
  } else {
    unique(x = idents.new)
  }
  old.levels <- levels(x = object)
  levels <- c(new.levels, old.levels)
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = idents(object = object))
  idents[samples] <- idents.new
  idents[is.na(x = idents)] <- 'NA'
  levels <- intersect(x = levels, y = unique(x = idents))
  names(x = idents) <- samples(object)
  missing.samples <- which(x = is.na(x = names(x = idents)))
  if (length(x = missing.samples) > 0){
    idents <- idents[-missing.samples]
  }
  idents <- factor(x = idents, levels = levels)
  slot(object = object, name = 'active.ident') <- idents
  if (drop){
    object <- droplevels(x = object)
  }
  return(object)
}


## 3.4 samples.Qbone ----
#' @rdname samples
#' @export
#'
samples.Qbone <- function(x){
  return(names(x@assays[[defaultAssay(x)]]@data))
}
# as S3 method
# samples.Qbone <- function(x){
#   return(names(x@assays[[defaultAssay(x)]][["samples"]]))
# }
# as S4 method
# setMethod("samples", "Qbone", function(x){
#   return(names(x@assays[[defaultAssay(x)]][["samples"]]))
# })

## 3.5 getQboneData.Qbone ----
#' @param assay Specific assay to get data from or set data for; defaults to
#' the \link[QboneObject:defaultAssay]{default assay}
#'
#' @rdname assayData
#' @export
#' @method getQboneData Qbone
#'
#' @order 3
#'
#' @concept data-access
#'
getQboneData.Qbone <- function(object, slot = 'data', assay = NULL, ...){
  CheckDots(...)
  object <- updateSlots(object = object)
  assay <- assay %||% defaultAssay(object = object)
  if (!assay %in% assays(object = object)){
    stop("'", assay, "' is not an assay", call. = FALSE)
  }
  return(getQboneData(
    object = object[[assay]],
    slot = slot
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. R-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 Qbone Methods ----
#' Qbone Methods
#'
#' Methods for \code{\link{Qbone}} objects for generics defined in other
#' packages
#'
#' @param x,object A \code{\link{Qbone}} object
#' @param i Depends on the method
#' \describe{
#  \item{\code{[}, \code{subset}}{Feature names or indices}
#  \item{\code{$}, \code{$<-}}{Name of a single metadata column}
#'  \item{\code{[[}, \code{[[<-}}{
#'   Name of one or more metadata columns or an associated object; associated
#'   objects include \code{\link{QboneData}}                  || double check
#   \code{\link{Graph}}, \code{\link{SeuratCommand}}, or
#   \code{\link{SpatialImage}} objects
#'  }
#' }
#' @param j,samples Sample names or indices
#' @param n The number of rows of metadata to return
#' @param ... Arguments passed to other methods
#'
#' @name Qbone-methods
#' @rdname Qbone-methods
#'
#'
NULL


### 4.1.1 [[.Qbone ----
#' @describeIn Qbone-methods Metadata and associated object accessor
#'
#' @param drop See \code{\link[base]{drop}}
#'
#' @return \code{[[}: If \code{i} is missing, the metadata data frame; if
#' \code{i} is a vector of metadata names, a data frame with the requested
#' metadata, otherwise, the requested associated object
#'
#' @export
#' @method [[ Qbone
#'
"[[.Qbone" <- function(x, i, ..., drop = FALSE){
  x <- updateSlots(object = x)
  if (missing(x = i)){
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  if (length(x = i) == 0){
    return(data.frame(row.names = names(getQboneData(x))))
  } else if (length(x = i) > 1 || any(i %in% colnames(x = slot(object = x, name = 'meta.data')))){
    if (any(!i %in% colnames(x = slot(object = x, name = 'meta.data')))){
      warning(
        "Cannot find the following bits of meta data: ",
        paste0(
          i[!i %in% colnames(x = slot(object = x, name = 'meta.data'))],
          collapse = ', '
        )
      )
    }
    i <- i[i %in% colnames(x = slot(object = x, name = 'meta.data'))]
    data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
    if (drop){
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = colnames(x = x), times = length(x = i))
    }
  } else {
    slot.use <- unlist(x = lapply(
      X = c('assays'), #, 'graphs', 'images'),
      FUN = function(s){
        if (any(i %in% names(x = slot(object = x, name = s)))){
          return(s)
        }
        return(NULL)
      }
    ))
    if (is.null(x = slot.use)){
      stop("Cannot find '", i, "' in this Qbone object", call. = FALSE)
    }
    data.return <- slot(object = x, name = slot.use)[[i]]
  }
  return(data.return)
}

## 4.2 levels.Qbone ----
#' @rdname idents
#' @export
#' @method levels Qbone
#'
#' @examples
#' \dontrun{
#' # Get the levels of identity classes of a Qbone object
#' levels(x = qbone1)
#' }
#'
levels.Qbone <- function(x) {
  x <- updateSlots(object = x)
  return(levels(x = idents(object = x)))
}

## 4.3 levels<- .Qbone ----
#' @rdname idents
#' @export
#' @method levels<- Qbone
#'
#' @examples
#' \dontrun{
#' # Reorder identity classes
#' levels(x = qbone1)
#' levels(x = qbone1) <- c('A', 'B')
#' levels(x = qbone1)
#' }
#'
"levels<-.Qbone" <- function(x, value) {
  x <- updateSlots(object = x)
  idents <- idents(object = x)
  if (!all(levels(x = idents) %in% value)) {
    stop("NA's generated by missing levels", call. = FALSE)
  }
  idents <- factor(x = idents, levels = value)
  idents(object = x) <- idents
  return(x)
}

#' @rdname idents
#' @export
#' @method droplevels Qbone
#'
droplevels.Qbone <- function(x, ...) {
  x <- updateSlots(object = x)
  slot(object = x, name = 'active.ident') <- droplevels(x = idents(object = x), ...)
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. S4 methods ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 5.1 [[<- Qbone ----
#' @describeIn Qbone-methods Metadata and associated object accessor
#'
#' @param value Additional metadata or associated objects to add; \strong{note}:
#' can pass \code{NULL} to remove metadata or an associated object
#'
#' @return \code{[[<-}: \code{x} with the metadata or associated objects added
#' as \code{i}; if \code{value} is \code{NULL}, removes metadata or associated
#' object \code{i} from object \code{x}
#'
#' @export
#'
setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'Qbone'),
  definition = function(x, i, ..., value){
    x <- updateSlots(object = x)
    # Require names, no index setting
    if (!is.character(x = i)){
      stop("'i' must be a character", call. = FALSE)
    }
    # Allow removing of other object
    if (is.null(x = value)){
      slot.use <- if (i %in% colnames(x = x[[]])){
        'meta.data'
      } else {
        FindObject(object = x, name = i)
      }
      if (is.null(x = slot.use)){
        stop("Cannot find object ", i, call. = FALSE)
      }
      if (i == defaultAssay(object = x)){
        stop("Cannot delete the default assay", call. = FALSE)
      }
    }
    # remove disallowed characters from object name
    newi <- if (is.null(x = value)){
      i
    } else {
      make.names(names = i)
    }
    if (any(i != newi)){
      warning(
        "Invalid name supplied, making object name syntactically valid. New object name is ",
        newi,
        "; see ?make.names for more details on syntax validity",
        call. = FALSE,
        immediate. = TRUE
      )
      i <- newi
    }
    # Figure out where to store data
    slot.use <- if (inherits(x = value, what = 'QboneData')){
      # Ensure we have the same number of samples
      if (length(value@data) != nrow(x = x@meta.data)){
        stop(
          "Cannot add a different number of samples than already present",
          call. = FALSE
        )
      }
      # Ensure cell order stays the same
      if (all(samples(x = value) %in% samples(x = x)) && !all(samples(x = value) == samples(x = x))){
        for (slot in c('data', 'scale.data')){
          assay.data <- getQboneData(object = value, slot = slot)
          if (length(assay.data!=0)){
            assay.data <- assay.data[samples(x = x), drop = FALSE]
          }
          # Use slot because SetAssayData is being weird
          slot(object = value, name = slot) <- assay.data
        }
      }
      'assays'
    # } else if (inherits(x = value, what = 'SpatialImage')){
    #   # Ensure that all samples for this image are present
    #   if (!all(samples(x = value) %in% samples(x = x))){
    #     stop("", call. = FALSE)
    #   }
    #   # Ensure Assay that SpatialImage is associated with is present in Seurat object
    #   if (!defaultAssay(object = value) %in% Assays(object = x)){
    #     warning(
    #       "Adding image data that isn't associated with any assay present",
    #       call. = FALSE,
    #       immediate. = TRUE
    #     )
    #   }
    #   'images'
    # } else if (inherits(x = value, what = 'Graph')){
    #   # Ensure Assay that Graph is associated with is present in the Seurat object
    #   if (is.null(x = defaultAssay(object = value))){
    #     warning(
    #       "Adding a Graph without an assay associated with it",
    #       call. = FALSE,
    #       immediate. = TRUE
    #     )
    #   } else if (!any(defaultAssay(object = value) %in% Assays(object = x))){
    #     stop("Cannot find assay '", defaultAssay(object = value), "' in this Qbone object", call. = FALSE)
    #   }
    #   # Ensure Graph object is in order
    #   if (all(samples(x = value) %in% samples(x = x)) && !all(samples(x = value) == samples(x = x))){
    #     value <- value[samples(x = x), samples(x = x)]
    #   }
    #   'graphs'
    # } else if (inherits(x = value, what = 'DimReduc')){
    #   # All DimReducs must be associated with an Assay
    #   if (is.null(x = defaultAssay(object = value))){
    #     stop("Cannot add a DimReduc without an assay associated with it", call. = FALSE)
    #   }
    #   # Ensure Assay that DimReduc is associated with is present in the Qbone object
    #   if (!IsGlobal(object = value) && !any(defaultAssay(object = value) %in% Assays(object = x))){
    #     stop("Cannot find assay '", defaultAssay(object = value), "' in this Qbone object", call. = FALSE)
    #   }
    #   # Ensure DimReduc object is in order
    #   if (all(samples(x = value) %in% samples(x = x)) && !all(samples(x = value) == samples(x = x))){
    #     slot(object = value, name = 'cell.embeddings') <- value[[samples(x = x), ]]
    #   }
    #   'reductions'
    # } else if (inherits(x = value, what = "Neighbor")){
    #   # Ensure all samples are present in the Qbone object
    #   if (length(x = samples(x = value)) > length(x = samples(x = x))){
    #     stop(
    #       "Cannot have more samples in the Neighbor object than are present in the Qbone object.",
    #       call. = FALSE
    #     )
    #   }
    #   if (!all(samples(x = value) %in% samples(x = x))){
    #     stop(
    #       "Cannot add samples in the Neighbor object that aren't present in the Qbone object.",
    #       call. = FALSE
    #     )
    #   }
    #   'neighbors'
    # } else if (inherits(x = value, what = 'SeuratCommand')){
    #   # Ensure Assay that SeuratCommand is associated with is present in the Qbone object
    #   if (is.null(x = defaultAssay(object = value))){
    #     warning(
    #       "Adding a command log without an assay associated with it",
    #       call. = FALSE,
    #       immediate. = TRUE
    #     )
    #   } else if (!any(defaultAssay(object = value) %in% Assays(object = x))){
    #     stop("Cannot find assay '", defaultAssay(object = value), "' in this Qbone object", call. = FALSE)
    #   }
    #   'commands'
    } else if (is.null(x = value)){
      slot.use
    } else {
      'meta.data'
    }
    if (slot.use == 'meta.data'){
      # Add data to object metadata
      meta.data <- x[[]]
      cell.names <- rownames(x = meta.data)
      # If we have metadata with names, ensure they match our order
      if (is.data.frame(x = value) && !is.null(x = rownames(x = value))){
        meta.order <- match(x = rownames(x = meta.data), table = rownames(x = value))
        value <- value[meta.order, , drop = FALSE]
      }
      if (length(x = i) > 1){
        # Add multiple pieces of metadata
        value <- rep_len(x = value, length.out = length(x = i))
        for (index in 1:length(x = i)){
          meta.data[i[index]] <- value[index]
        }
      } else {
        # Add a single column to metadata
        if (length(x = intersect(x = names(x = value), y = cell.names)) > 0){
          meta.data[, i] <- value[cell.names]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)){
          meta.data[, i] <- value
        } else {
          stop("Cannot add more or fewer sample meta.data information without values being named with sample names", call. = FALSE)
        }
      }
      # Check to ensure that we aren't adding duplicate names
      if (any(colnames(x = meta.data) %in% FilterObjects(object = x))){
        bad.cols <- colnames(x = meta.data)[which(colnames(x = meta.data) %in% FilterObjects(object = x))]
        stop(
          paste0(
            "Cannot add a metadata column with the same name as an Assay or DimReduc - ",
            paste(bad.cols, collapse = ", ")),
          call. = FALSE
        )
      }
      # Store the revised metadata
      slot(object = x, name = 'meta.data') <- meta.data
    } else {
      # Add other object to Qbone object
      # Ensure samples match in value and order
      if (!inherits(x = value, what = c('SeuratCommand', 'NULL', 'SpatialImage', 'Neighbor')) && !all(samples(x = value) == samples(x = x))){
        stop("All samples in the object being added must match the samples in this object", call. = FALSE)
      }
      # Ensure we're not duplicating object names
      duplicate <- !is.null(x = FindObject(object = x, name = i)) &&
        !inherits(x = value, what = c(class(x = x[[i]]), 'NULL')) &&
        !inherits(x = x[[i]], what = class(x = value))
      if (isTRUE(x = duplicate)){
        stop(
          "This object already contains ",
          i,
          " as a",
          ifelse(
            test = tolower(x = substring(text = class(x = x[[i]]), first = 1, last = 1)) %in% c('a', 'e', 'i', 'o', 'u'),
            yes = 'n ',
            no = ' '
          ),
          class(x = x[[i]]),
          ", so ",
          i,
          " cannot be used for a ",
          class(x = value),
          call. = FALSE
        )
      }
      #   # Check keyed objects
      #   if (inherits(x = value, what = c('Assay', 'DimReduc', 'SpatialImage'))){
      #     if (length(x = Key(object = value)) == 0 || nchar(x = Key(object = value)) == 0){
      #       Key(object = value) <- paste0(tolower(x = i), '_')
      #     }
      #     Key(object = value) <- UpdateKey(key = Key(object = value))
      #     # Check for duplicate keys
      #     object.keys <- Key(object = x)
      #     vkey <- Key(object = value)
      #     if (vkey %in% object.keys && !isTRUE(x = object.keys[i] == vkey)){
      #       new.key <- if (is.na(x = object.keys[i])){
      #         # Attempt to create a duplicate key based off the name of the object being added
      #         new.keys <- paste0(
      #           paste0(tolower(x = i), c('', RandomName(length = 2L))),
      #           '_'
      #         )
      #         # Select new key to use
      #         key.use <- min(which(x = !new.keys %in% object.keys))
      #         new.key <- if (is.infinite(x = key.use)){
      #           RandomName(length = 17L)
      #         } else {
      #           new.keys[key.use]
      #         }
      #         warning(
      #           "Cannot add objects with duplicate keys (offending key: ",
      #           Key(object = value),
      #           "), setting key to '",
      #           new.key,
      #           "'",
      #           call. = FALSE
      #         )
      #         new.key
      #       } else {
      #         # Use existing key
      #         warning(
      #           "Cannot add objects with duplicate keys (offending key: ",
      #           Key(object = value),
      #           ") setting key to original value '",
      #           object.keys[i],
      #           "'",
      #           call. = FALSE
      #         )
      #         object.keys[i]
      #       }
      #       # Set new key
      #       Key(object = value) <- new.key
      #     }
      #   }

        # # For Assays, run CalcN
        # if (inherits(x = value, what = 'QboneData')){
        #   if ((!i %in% assays(object = x)) |
        #       (i %in% assays(object = x) && !identical(
        #         x = getQboneData(object = x, assay = i, slot = "data"),
        #         y = getQboneData(object = value, slot = "data"))
        #       )){
        #     n.calc <- CalcN(object = value)
        #     if (!is.null(x = n.calc)){
        #       names(x = n.calc) <- paste(names(x = n.calc), i, sep = '_')
        #       x[[names(x = n.calc)]] <- n.calc
        #     }
        #   }
        # }

        # When removing an Assay, clear out associated DimReducs, Graphs, and SeuratCommands
        if (is.null(x = value) && inherits(x = x[[i]], what = 'Assay')){
          objs.assay <- FilterObjects(
            object = x,
            classes.keep = c('DimReduc', 'SeuratCommand', 'Graph')
          )
          objs.assay <- Filter(
            f = function(o){
              return(all(defaultAssay(object = x[[o]]) == i) && !IsGlobal(object = x[[o]]))
            },
            x = objs.assay
          )
          for (o in objs.assay){
            x[[o]] <- NULL
          }
        }
        # If adding a command, ensure it gets put at the end of the command list
        if (inherits(x = value, what = 'SeuratCommand')){
          slot(object = x, name = slot.use)[[i]] <- NULL
          slot(object = x, name = slot.use) <- Filter(
            f = Negate(f = is.null),
            x = slot(object = x, name = slot.use)
          )
        }
        slot(object = x, name = slot.use)[[i]] <- value
        slot(object = x, name = slot.use) <- Filter(
          f = Negate(f = is.null),
          x = slot(object = x, name = slot.use)
        )
    }
    # CheckGC()
    gc(verbose = FALSE)
    return(x)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Internal ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 6.1 FindObject ----
#' Find the collection of an object within a Qbone object
#' Adopted from Seurateobject package
#'
#' @param object A \code{\link{Qbone}} object
#' @param name Name of object to find
#'
#' @return The collection (slot) of the object
#'
#' @keywords internal
#'
#' @noRd
#'
FindObject <- function(object, name){
  collections <- c(
    'assays'
    # 'graphs',
    # 'commands',
    # 'images'
  )
  object.names <- lapply(
    X = collections,
    FUN = function(x){
      return(names(x = slot(object = object, name = x)))
    }
  )
  names(x = object.names) <- collections
  object.names <- Filter(f = Negate(f = is.null), x = object.names)
  for (i in names(x = object.names)){
    if (name %in% names(x = slot(object = object, name = i))){
      return(i)
    }
  }
  return(NULL)
}

## 6.2 FilterObjects ----
#'
#' Get the names of objects within a \code{Qbone} object that are of a
#' certain class
#' Adopted from Seurateobject package
#'
#' @param object A \code{\link{Qbone}} object
#' @param classes.keep A vector of names of classes to get
#'
#' @return A vector with the names of objects within the \code{Qbone} object
#' that are of class \code{classes.keep}
#'
#' @keywords internal
#'
#' @noRd
#'
FilterObjects <- function(object, classes.keep = c('QboneData')){
  object <- updateSlots(object = object)
  slots <- na.omit(object = Filter(
    f = function(x){
      sobj <- slot(object = object, name = x)
      return(is.list(x = sobj) && !is.data.frame(x = sobj) && !is.package_version(x = sobj))
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x){
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i){
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
  return(names(x = object.classes))
}
