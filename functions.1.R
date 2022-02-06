library("Rcpp")

# Get Script dir----
# script.dir <- dirname(sys.frame(1)$ofile)
# getSrcDirectory(function(x) {x})
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (cmdArgs[1] == "RStudio") {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  } else if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

thisFile()

setwd(thisFile())


# test static_assert not really working----

sourceCpp("cmatrix_2.cpp")

A <- matrix(rnorm(10000), 100, 100)
B <- matrix(rnorm(10000), 100, 100)


o1 = A%*%B

e1 = eigenMatMult(A, B)
e2 = eigenMapMatMult(A, B)

o3 = A%*%t(B)
e3 = eigenMapMatMulttrans(A,B)
e4 = eigenMatMulttrans(A,B)

# test class example ----
studentBio <- list(studentName = "Harry Potter", studentAge = 19, studentContact="London")
class(studentBio) <- "StudentInfo"
studentBio

contact <- function(object) {
  UseMethod("contact")
}

contact.StudentInfo <- function(object) {
  cat("Your contact is", object$studentContact, "\n")
}

contact.StudentInfo <- function(object) {
  cat("Your contact is", object$studentName, "\n")
}
contact(studentBio)
# define class ----
# https://rdrr.io/cran/SeuratObject/src/R/seurat.R
Qbone <- setClass(
  Class = 'Qbone',
  slots = c(
    assays = 'list',
    # raw.data = 'list',
    # thin.data = 'list',
    meta.data = 'data.frame',
    # thin.meta = 'list',
    active.assay = 'character',
    active.ident = 'factor',
    # lasso.list = 'list',
    graphs = 'list',
    images = 'list',
    project.name = 'character',
    misc = 'list',
    # version = 'package_version',
    commands = 'list',
    tools = 'list'
  ), 
  prototype = c(
    assays = list(),
    # raw.data = list(),
    # thin.data = list(),
    # meta.data = data.frame(id = NULL),
    # thin.meta = list(), # list(T, "0.1")
    active.assay = NA_character_,
    active.ident = factor(),
    # lasso.list = list(),
    graphs = list(),
    images = list(),
    project.name = "Qbone",
    misc = list(),
    # version = 'package_version',
    commands = list(),
    tools = list()
  )
)

q0 = new("Qbone")

# https://adv-r.hadley.nz/s4.html#helper
# CreateSeuratObject.default <- 

# generic note----
# standardGeneric() is the S4 equivalent to UseMethod().

# It is bad practice to use {} in the generic as it triggers a special case that is more expensive, and generally best avoided.


# Read Qbone ----
library("methods")
# For this reason, it’s a good idea to include an explicit library(methods) whenever you’re using S4.
library("fpeek")
ReadQbone0 <- function(
  data.dir,
  groupbyfolder = F, 
  data.column = 1,
  # skip = 1,
  # header = F,
  project = "QboneProject"
  # gene.column = 2,
  # cell.column = 1,
  # unique.features = TRUE,
  # strip.suffix = FALSE
) {
  full.data <- list()
  meta.file <- c()
  meta.group <- c()
  # if (groupbyfolder == T){
  #   meta.group <- c()
  # }
  # check dir/file existence
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    file.list = list.files(path = paste0(data.dir), pattern = ".csv$", recursive = TRUE)
    if (length(file.list) ==0){
      stop("Directory provided does not contain any .csv file")
    }
  }
  # read file & sort
  for (i in 1:length(file.list)){
    # no need for V3
    # peek.n10 = textConnection(peek_head(paste0(data.dir,"/",file.list[i]), n = 10, intern = TRUE)[-1])
    # col.n = max(count.fields(peek.n10, sep = ","))
    
    file1 =  drop(as.matrix(drop(as.matrix(
      # V1
      # read.csv(paste0(data.dir,"/",file.list[i]),
      #          skip = skip , header = header,
      #          colClasses = c(rep("NULL", c(1:col.n)[data.column]-1),
      #                         rep("numeric", 1),
      #                         rep("NULL", (col.n-c(1:col.n)[data.column]))))
      # V2
      # readr::read_csv(paste0(data.dir,"/",file.list[i]),
      #                 skip = skip , col_names = header,
      #                 col_types = paste0(
      #                   stringr::str_dup("-", c(1:col.n)[data.column]-1),
      #                   stringr::str_dup("n", 1),
      #                   stringr::str_dup("-", (col.n-c(1:col.n)[data.column]))
      #                 )
      # V3
      data.table::fread(paste0(data.dir,"/",file.list[i]), 
                        # skip = skip , 
                        select=c(data.column)
      )
    ))))
    full.data[[i]] <- sort(file1)
    meta.file <- c(meta.file, gsub(".csv", "", basename(file.list[i])))
    if (groupbyfolder == T){
      meta.group <- c(meta.group, dirname(file.list[i]))
    }
  }
  if (groupbyfolder == T){
    meta.data = data.frame(id = meta.file,
                           group = meta.group)
  } else {
    meta.data = data.frame(id = meta.file)
  }
  object <- new(
    Class = 'Qbone',
    raw.data = full.data,
    meta.data = meta.data,
    project.name = project
    # version = packageVersion(pkg = 'Qbone'),
  )
  return(object)
}

# if groupbyfolder ? dirname(file.list)
# file name gsub(".csv", "", basename(file.list[i]))
# data raw.data thin.data
# Ident

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv"
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"
q1 = ReadQbone0(data.dir)
q1[["id"]]

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group2"
q2 = ReadQbone0(data.dir, groupbyfolder = F)
q3 = ReadQbone0(data.dir, groupbyfolder = T)

# need utils and mthods ----

## Subset method ----
# https://rdrr.io/cran/SeuratObject/man/Seurat-methods.html
# "[[.Seurat"  in seurat.R
# https://rdrr.io/cran/SeuratObject/src/R/seurat.R

# ../R/utils.R ----
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# SeuratObject  /R/seurat.R ----

UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(mode = class(x = xobj), length = 1L)
    }
  }
  return(object)
}

"[[.Qbone" <- function(x, i, ..., drop = FALSE) {
  x <- UpdateSlots(object = x)
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  if (length(x = i) == 0) {
    return(data.frame(row.names = colnames(x = x)))
  } else if (length(x = i) > 1 || any(i %in% colnames(x = slot(object = x, name = 'meta.data')))) {
    if (any(!i %in% colnames(x = slot(object = x, name = 'meta.data')))) {
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
    if (drop) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = colnames(x = x), times = length(x = i))
    }
  } else {
    slot.use <- unlist(x = lapply(
      # X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands', 'images'),
      X = c('assays', 'graphs', 'images'),
      FUN = function(s) {
        if (any(i %in% names(x = slot(object = x, name = s)))) {
          return(s)
        }
        return(NULL)
      }
    ))
    if (is.null(x = slot.use)) {
      stop("Cannot find '", i, "' in this Qbone object", call. = FALSE)
    }
    data.return <- slot(object = x, name = slot.use)[[i]]
  }
  return(data.return)
}

# slotNames(object)

FindObject <- function(object, name) {
  collections <- c(
    'assays',
    'graphs',
    'commands',
    'images'
  )
  object.names <- lapply(
    X = collections,
    FUN = function(x) {
      return(names(x = slot(object = object, name = x)))
    }
  )
  names(x = object.names) <- collections
  object.names <- Filter(f = Negate(f = is.null), x = object.names)
  for (i in names(x = object.names)) {
    if (name %in% names(x = slot(object = object, name = i))) {
      return(i)
    }
  }
  return(NULL)
}

FilterObjects <- function(object, classes.keep = c('Assay')) {
  object <- UpdateSlots(object = object)
  slots <- na.omit(object = Filter(
    f = function(x) {
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
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(inherits(x = object[[i]], what = classes.keep))
    }
  )
  object.classes <- which(x = object.classes, useNames = TRUE)
  return(names(x = object.classes))
}

DefaultAssay.Qbone <- function(object, ...) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'active.assay'))
}

DefaultAssay(cqo2)

"DefaultAssay<-.Qbone" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  if (!value %in% names(x = slot(object = object, name = 'assays'))) {
    stop("Cannot find assay ", value)
  }
  slot(object = object, name = 'active.assay') <- value
  return(object)
}

CheckGC <- function(option = 'Qbone.memsafe') {
  if (isTRUE(x = getOption(x = option, default = FALSE))) {
    gc(verbose = FALSE)
  }
  return(invisible(x = NULL))
}

CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(
    X = fxn.args,
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  )
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found", call. = FALSE)
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Qbone.checkdots", default = 'warn'),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Qbone.checkdots option. Please choose one of warn, stop, silent")
      )
      # unused.hints <- sapply(X = unused, FUN = OldParamHints)
      # names(x = unused.hints) <- unused
      # unused.hints <- na.omit(object = unused.hints)
      # if (length(x = unused.hints) > 0) {
      #   message(
      #     "Suggested parameter: ",
      #     paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
      #     "\n"
      #   )
      # }
    }
  }
  return(invisible(x = NULL))
}

setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'Qbone'),
  definition = function(x, i, ..., value) {
    x <- UpdateSlots(object = x)
    # Require names, no index setting
    if (!is.character(x = i)) {
      stop("'i' must be a character", call. = FALSE)
    }
    # Allow removing of other object
    if (is.null(x = value)) {
      slot.use <- if (i %in% colnames(x = x[[]])) {
        'meta.data'
      } else {
        FindObject(object = x, name = i)
      }
      if (is.null(x = slot.use)) {
        stop("Cannot find object ", i, call. = FALSE)
      }
      if (i == DefaultAssay(object = x)) {
        stop("Cannot delete the default assay", call. = FALSE)
      }
    }
    # remove disallowed characters from object name
    newi <- if (is.null(x = value)) {
      i
    } else {
      make.names(names = i)
    }
    if (any(i != newi)) {
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
    slot.use <- if (inherits(x = value, what = 'Assay')) {
      # Ensure we have the same number of cells
      if (ncol(x = value) != ncol(x = x)) {
        stop(
          "Cannot add a different number of cells than already present",
          call. = FALSE
        )
      }
      # Ensure cell order stays the same
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        for (slot in c('counts', 'data', 'scale.data')) {
          assay.data <- GetAssayData(object = value, slot = slot)
          if (!IsMatrixEmpty(x = assay.data)) {
            assay.data <- assay.data[, Cells(x = x), drop = FALSE]
          }
          # Use slot because SetAssayData is being weird
          slot(object = value, name = slot) <- assay.data
        }
      }
      'assays'
    } else if (inherits(x = value, what = 'SpatialImage')) {
      # Ensure that all cells for this image are present
      if (!all(Cells(x = value) %in% Cells(x = x))) {
        stop("", call. = FALSE)
      }
      # Ensure Assay that SpatialImage is associated with is present in Seurat object
      if (!DefaultAssay(object = value) %in% Assays(object = x)) {
        warning(
          "Adding image data that isn't associated with any assay present",
          call. = FALSE,
          immediate. = TRUE
        )
      }
      'images'
    } else if (inherits(x = value, what = 'Graph')) {
      # Ensure Assay that Graph is associated with is present in the Seurat object
      if (is.null(x = DefaultAssay(object = value))) {
        warning(
          "Adding a Graph without an assay associated with it",
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (!any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Qbone object", call. = FALSE)
      }
      # Ensure Graph object is in order
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        value <- value[Cells(x = x), Cells(x = x)]
      }
      'graphs'
    } else if (inherits(x = value, what = 'DimReduc')) {
      # All DimReducs must be associated with an Assay
      if (is.null(x = DefaultAssay(object = value))) {
        stop("Cannot add a DimReduc without an assay associated with it", call. = FALSE)
      }
      # Ensure Assay that DimReduc is associated with is present in the Qbone object
      if (!IsGlobal(object = value) && !any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Qbone object", call. = FALSE)
      }
      # Ensure DimReduc object is in order
      if (all(Cells(x = value) %in% Cells(x = x)) && !all(Cells(x = value) == Cells(x = x))) {
        slot(object = value, name = 'cell.embeddings') <- value[[Cells(x = x), ]]
      }
      'reductions'
    } else if (inherits(x = value, what = "Neighbor")) {
      # Ensure all cells are present in the Qbone object
      if (length(x = Cells(x = value)) > length(x = Cells(x = x))) {
        stop(
          "Cannot have more cells in the Neighbor object than are present in the Qbone object.",
          call. = FALSE
        )
      }
      if (!all(Cells(x = value) %in% Cells(x = x))) {
        stop(
          "Cannot add cells in the Neighbor object that aren't present in the Qbone object.",
          call. = FALSE
        )
      }
      'neighbors'
    } else if (inherits(x = value, what = 'SeuratCommand')) {
      # Ensure Assay that SeuratCommand is associated with is present in the Qbone object
      if (is.null(x = DefaultAssay(object = value))) {
        warning(
          "Adding a command log without an assay associated with it",
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (!any(DefaultAssay(object = value) %in% Assays(object = x))) {
        stop("Cannot find assay '", DefaultAssay(object = value), "' in this Qbone object", call. = FALSE)
      }
      'commands'
    } else if (is.null(x = value)) {
      slot.use
    } else {
      'meta.data'
    }
    if (slot.use == 'meta.data') {
      # Add data to object metadata
      meta.data <- x[[]]
      cell.names <- rownames(x = meta.data)
      # If we have metadata with names, ensure they match our order
      if (is.data.frame(x = value) && !is.null(x = rownames(x = value))) {
        meta.order <- match(x = rownames(x = meta.data), table = rownames(x = value))
        value <- value[meta.order, , drop = FALSE]
      }
      if (length(x = i) > 1) {
        # Add multiple pieces of metadata
        value <- rep_len(x = value, length.out = length(x = i))
        for (index in 1:length(x = i)) {
          meta.data[i[index]] <- value[index]
        }
      } else {
        # Add a single column to metadata
        if (length(x = intersect(x = names(x = value), y = cell.names)) > 0) {
          meta.data[, i] <- value[cell.names]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)) {
          meta.data[, i] <- value
        } else {
          stop("Cannot add more or fewer cell meta.data information without values being named with cell names", call. = FALSE)
        }
      }
      # Check to ensure that we aren't adding duplicate names
      if (any(colnames(x = meta.data) %in% FilterObjects(object = x))) {
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
      # Ensure cells match in value and order
      if (!inherits(x = value, what = c('SeuratCommand', 'NULL', 'SpatialImage', 'Neighbor')) && !all(Cells(x = value) == Cells(x = x))) {
        stop("All cells in the object being added must match the cells in this object", call. = FALSE)
      }
      # Ensure we're not duplicating object names
      duplicate <- !is.null(x = FindObject(object = x, name = i)) &&
        !inherits(x = value, what = c(class(x = x[[i]]), 'NULL')) &&
        !inherits(x = x[[i]], what = class(x = value))
      if (isTRUE(x = duplicate)) {
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
    #   if (inherits(x = value, what = c('Assay', 'DimReduc', 'SpatialImage'))) {
    #     if (length(x = Key(object = value)) == 0 || nchar(x = Key(object = value)) == 0) {
    #       Key(object = value) <- paste0(tolower(x = i), '_')
    #     }
    #     Key(object = value) <- UpdateKey(key = Key(object = value))
    #     # Check for duplicate keys
    #     object.keys <- Key(object = x)
    #     vkey <- Key(object = value)
    #     if (vkey %in% object.keys && !isTRUE(x = object.keys[i] == vkey)) {
    #       new.key <- if (is.na(x = object.keys[i])) {
    #         # Attempt to create a duplicate key based off the name of the object being added
    #         new.keys <- paste0(
    #           paste0(tolower(x = i), c('', RandomName(length = 2L))),
    #           '_'
    #         )
    #         # Select new key to use
    #         key.use <- min(which(x = !new.keys %in% object.keys))
    #         new.key <- if (is.infinite(x = key.use)) {
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
    #   # For Assays, run CalcN
    #   if (inherits(x = value, what = 'Assay')) {
    #     if ((!i %in% Assays(object = x)) |
    #         (i %in% Assays(object = x) && !identical(
    #           x = GetAssayData(object = x, assay = i, slot = "counts"),
    #           y = GetAssayData(object = value, slot = "counts"))
    #         )) {
    #       n.calc <- CalcN(object = value)
    #       if (!is.null(x = n.calc)) {
    #         names(x = n.calc) <- paste(names(x = n.calc), i, sep = '_')
    #         x[[names(x = n.calc)]] <- n.calc
    #       }
    #     }
    #   }
    #   # When removing an Assay, clear out associated DimReducs, Graphs, and SeuratCommands
    #   if (is.null(x = value) && inherits(x = x[[i]], what = 'Assay')) {
    #     objs.assay <- FilterObjects(
    #       object = x,
    #       classes.keep = c('DimReduc', 'SeuratCommand', 'Graph')
    #     )
    #     objs.assay <- Filter(
    #       f = function(o) {
    #         return(all(DefaultAssay(object = x[[o]]) == i) && !IsGlobal(object = x[[o]]))
    #       },
    #       x = objs.assay
    #     )
    #     for (o in objs.assay) {
    #       x[[o]] <- NULL
    #     }
    #   }
    #   # If adding a command, ensure it gets put at the end of the command list
    #   if (inherits(x = value, what = 'SeuratCommand')) {
    #     slot(object = x, name = slot.use)[[i]] <- NULL
    #     slot(object = x, name = slot.use) <- Filter(
    #       f = Negate(f = is.null),
    #       x = slot(object = x, name = slot.use)
    #     )
    #   }
    #   slot(object = x, name = slot.use)[[i]] <- value
    #   slot(object = x, name = slot.use) <- Filter(
    #     f = Negate(f = is.null),
    #     x = slot(object = x, name = slot.use)
    #   )
    }
    # CheckGC()
    gc(verbose = FALSE)
    return(x)
  }
)

cqo2[["newneame"]] <- c("a1", "b1")

# inherits is faster than method::is

# https://rdrr.io/cran/SeuratObject/src/R/zzz.R
.AddMetaData <- function(object, metadata, col.name = NULL) {
  if (is.null(x = col.name) && is.atomic(x = metadata)) {
    stop("'col.name' must be provided for atomic metadata types (eg. vectors)")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata) # `%||%`  in utilities.R https://rdrr.io/cran/Seurat/src/R/utilities.R
  if (is.null(x = col.name)) {
    stop("No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  return(object)
}

AddMetaData <- .AddMetaData # https://rdrr.io/cran/SeuratObject/src/R/seurat.R

AddMetaData(q0, c(1:2), col.name = "number") 


ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(
    x = as.character(x = field),
    split = ","
  )))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(
    strsplit(x = string, split = delim)[[1]][fields],
    collapse = delim
  ))
}

ExtractField("abc1")

# CreateSeuratObject.Assay <- function
# https://rdrr.io/cran/SeuratObject/src/R/seurat.R
# Assay <- setClass( https://rdrr.io/cran/SeuratObject/src/R/assay.R

CreateQboneObject <- function(
  data = data,
  project = 'QboneProject',
  assay = "Bone",
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  ...
){
  # if (is.null(assay)){
  #   stop('"assay =" is missing, with no default. Please provide a name of your assay ("Bone", "Image" etc.)')
  # }
  if (missing(x = data)) {
    stop("Must provide 'data'")
  }
  if (is.list(data)){
    if (!is.null(x = meta.data)) {
      if (nrow(x = meta.data) != length(data)) {
        stop("There is a mismatch between the number of Metadata and the number of input data.")
      }
    } else {
      warning("Metadata is not provided.")
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
  # Check assay key
  
  names(data) <- meta.data[,1]
  samples <- list(data)
  names(x = samples) <- "samples"
  assay.list <- list(samples)
  names(x = assay.list) <- assay
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    # X = colnames(x = meta.data),
    X = names(data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  if (any(is.na(x = idents))) {
    warning(
      "Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = length(data))
  }
  names(x = idents) <- names(data)
  
  object <- new(
    Class = 'Qbone',
    assays = assay.list,
    meta.data = meta.data,
    active.assay = assay,
    active.ident = idents,
    project.name = project,
    # version = packageVersion(pkg = 'Qbone'),
    ...
  )
  # if (!is.null(x = meta.data)) {
  #   object <- AddMetaData(object = object, metadata = meta.data)
  # }
}



class(c(data,data))
list1 = list(c("123", "345"), "234")
cqo1 = CreateQboneObject(data = list1,
                         meta.data = data.frame(name = c("a_1", "b_1"))
                         )

cqo2 = AddMetaData(cqo1, c("c", "d"), col.name = "names")

cqo2[["name"]] <- data.frame(name = c("a", "b"))

names(list1) <- c("name1", "name2")
list1 = list(q1@raw.data)
names(list1)
i = 2
list1[[i]]
list2 = list(list1 = list1, list2 = list1)
length(list1)
is.atomic(q1@raw.data[[1]])
is.recursive(q1@raw.data[[1]])

# Read10X 
# https://rdrr.io/cran/Seurat/src/R/preprocessing.R
ReadQbone <- function(
  data.dir,
  groupbyfolder = F, 
  data.column = 1
){
  full.data <- list()
  meta.file <- c()
  meta.group <- c()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    file.list = list.files(path = paste0(data.dir), pattern = ".csv$", recursive = TRUE)
    if (length(file.list) ==0){
      stop("Directory provided does not contain any .csv file")
    }
  }
  # read file & sort
  for (i in 1:length(file.list)){
    file1 =  drop(as.matrix(drop(as.matrix(
      data.table::fread(paste0(data.dir,"/",file.list[i]),
                        select=c(data.column)
      )
    ))))
    full.data[[i]] <- file1 # sort(file1)
    meta.file <- c(meta.file, gsub(".csv", "", basename(file.list[i])))
    if (groupbyfolder == T){
      meta.group <- c(meta.group, dirname(file.list[i]))
    }
  }
  if (groupbyfolder == T){
    meta.data = data.frame(id = meta.file,
                           group = meta.group)
  } else {
    meta.data = data.frame(id = meta.file)
  }
  object <- CreateQboneObject(data = full.data, meta.data = meta.data)
  return(object)
}

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"
q1 = ReadQbone(data.dir, groupbyfolder = T)

q1[['id']]

q1[['id']] <- c("03_VTK", "04_VTK")

Idents.Qbone <- function(object, ...) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'active.ident'))
}

Idents(cqo1)

Samples.Qbone <- function(x) {
  # names(x@assays[[x@active.assay]])
  # names(x@assays[[DefaultAssay(x)]])
  return(names(x@assays[[DefaultAssay(x)]][["samples"]]))
}
cqo2@active.assay
Samples(cqo2)

"Idents<-.Qbone" <- function(object, samples = NULL, drop = FALSE, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  samples <- samples %||% Samples(object)
  if (is.numeric(x = samples)) {
    samples <- Samples(object)[samples]
  }
  samples <- intersect(x = samples, y = Samples(object))
  samples <- match(x = samples, table = Samples(object))
  if (length(x = samples) == 0) {
    warning("Cannot find samples provided")
    return(object)
  }
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[[]])) {
    unlist(x = object[[value]], use.names = FALSE)[samples]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = samples))
  }
  new.levels <- if (is.factor(x = idents.new)) {
    levels(x = idents.new)
  } else {
    unique(x = idents.new)
  }
  old.levels <- levels(x = object)
  levels <- c(new.levels, old.levels)
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[samples] <- idents.new
  idents[is.na(x = idents)] <- 'NA'
  levels <- intersect(x = levels, y = unique(x = idents))
  names(x = idents) <- Samples(object)
  missing.samples <- which(x = is.na(x = names(x = idents)))
  if (length(x = missing.samples) > 0) {
    idents <- idents[-missing.samples]
  }
  idents <- factor(x = idents, levels = levels)
  slot(object = object, name = 'active.ident') <- idents
  if (drop) {
    object <- droplevels(x = object)
  }
  return(object)
}

Idents(cqo2)<- "names"

# IntegrateData <- function(
# https://rdrr.io/cran/Seurat/src/R/integration.R

# test speed ReadQbone ----
library("microbenchmark")

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"

nopeek <- function(
  data.dir,
  groupbyfolder = F, 
  data.column = 1,
  skip = 1,
  header = F,
  project = "QboneProject"
  # gene.column = 2,
  # cell.column = 1,
  # unique.features = TRUE,
  # strip.suffix = FALSE
) {
  full.data <- list()
  meta.file <- c()
  meta.group <- c()
  # if (groupbyfolder == T){
  #   meta.group <- c()
  # }
  # check dir/file existence
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    file.list = list.files(path = paste0(data.dir), pattern = ".csv$", recursive = TRUE)
    if (length(file.list) ==0){
      stop("Directory provided does not contain any .csv file")
    }
  }
  # read file & sort
  for (i in 1:length(file.list)){
    # peek.n10 = textConnection(peek_head(paste0(data.dir,"/",file.list[i]), n = 10, intern = TRUE)[-1])
    col.n = max(count.fields(paste0(data.dir,"/",file.list[i]), sep = ","))
    
    file1 =  drop(as.matrix(drop(as.matrix(
      read.csv(paste0(data.dir,"/",file.list[i]),
               skip = skip , header = header,
               colClasses = c(rep("NULL", c(1:col.n)[data.column]-1),
                              rep("numeric", 1),
                              rep("NULL", (col.n-c(1:col.n)[data.column]))))
    ))))
    full.data[[i]] <- sort(file1)
    meta.file <- c(meta.file, gsub(".csv", "", basename(file.list[i])))
    if (groupbyfolder == T){
      meta.group <- c(meta.group, dirname(file.list[i]))
    }
  }
  if (groupbyfolder == T){
    meta.data = data.frame(id = meta.file,
                           group = meta.group)
  } else {
    meta.data = data.frame(id = meta.file)
  }
  object <- new(
    Class = 'Qbone',
    raw.data = full.data,
    meta.data = meta.data,
    lasso.list = list(),
    project.name = project,
    # version = packageVersion(pkg = 'Qbone'),
    misc = list(),
    commands = list(),
    tools = list()
  )
  return(object)
}

ReadQbone(data.dir)
nopeek(data.dir)
microbenchmark(ReadQbone(data.dir),
               nopeek(data.dir))
# v1 read.csv
# Unit: seconds
#                 expr       min       lq      mean    median        uq       max neval
# ReadQbone(data.dir)  6.295032  6.41334  6.542009  6.547169  6.663608  6.859548   100
# nopeek(data.dir)    11.677127 11.89641 12.061280 12.019549 12.215060 12.640033   100

# v2 read_csv
# Unit: milliseconds                                                                                                               
# expr        min         lq       mean     median         uq        max neval
# ReadQbone(data.dir)   603.3834   671.3038   695.3638   687.6152   707.7549   954.6117   100
# nopeek(data.dir) 11749.8619 11855.3312 11999.2320 12017.4779 12083.4712 12329.0297   100

## readr::read_csv() ----

microbenchmark(read.csv(paste0(data.dir,"03_VTK_IO.csv")),
               readr::read_csv(paste0(data.dir,"03_VTK_IO.csv"))
)
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"
data.column = 1
skip = 1
header = F
peek.n10 = textConnection(peek_head(paste0(data.dir,"/03_VTK_IO.csv"), n = 10, intern = TRUE)[-1])
col.n = max(count.fields(peek.n10, sep = ","))

read1 = readr::read_csv(paste0(data.dir,"/03_VTK_IO.csv"), 
                        skip = skip , col_names = header,
                        col_types = paste0(
                          stringr::str_dup("-", c(1:col.n)[data.column]-1),
                          stringr::str_dup("n", 1),
                          stringr::str_dup("-", (col.n-c(1:col.n)[data.column]))
                        )
                        )

read2 = read.csv(paste0(data.dir,"/03_VTK_IO.csv"),
                 skip = skip , header = header,
                 colClasses = c(rep("NULL", c(1:col.n)[data.column]-1),
                                rep("numeric", 1),
                                rep("NULL", (col.n-c(1:col.n)[data.column]))))
all.equal(read1[,1], read2[,1])
identical(read1[,1], read2[,1])
file1 =  drop(as.matrix(drop(as.matrix(
  # readr::read_csv()
  readr::read_csv(paste0(data.dir,"/03_VTK_IO.csv"), 
                  skip = skip , col_names = header,
                  col_types = paste0(
                    stringr::str_dup("-", c(1:col.n)[data.column]-1),
                    stringr::str_dup("n", 1),
                    stringr::str_dup("-", (col.n-c(1:col.n)[data.column]))
                  )
)))))
class(file1)
file2 =  drop(as.matrix(drop(as.matrix(
  # readr::read_csv()
  read.csv(paste0(data.dir,"/03_VTK_IO.csv"),
           skip = skip , header = header,
           colClasses = c(rep("NULL", c(1:col.n)[data.column]-1),
                          rep("numeric", 1),
                          rep("NULL", (col.n-c(1:col.n)[data.column]))))
                  )
  )))
class(file2)
all.equal(file1, file2)
identical(file1, file2)
sum(file1 != file2)
file1.1 = file1[(file1 != file2)]
file2.1 = file2[(file1 != file2)]
sprintf("%.54f",file1.1[1])
sprintf("%.54f",file2.1[1])

## data.table::fread() ----

read3 = data.table::fread(paste0(data.dir,"/03_VTK_IO.csv"), 
                        # skip = skip , 
                        select=c(data.column)
)
file3 =  drop(as.matrix(drop(as.matrix(
  data.table::fread(paste0(data.dir,"/03_VTK_IO.csv"), 
                    # skip = skip , 
                    select=c(data.column)
  )
))))
all.equal(file1, file3)
identical(file1, file3)

microbenchmark(data.table::fread(paste0(data.dir,"/03_VTK_IO.csv"),
                                 select=c(data.column)),
               readr::read_csv(paste0(data.dir,"/03_VTK_IO.csv"),
                               skip = skip , col_names = header,
                               col_types = paste0(
                                 stringr::str_dup("-", c(1:col.n)[data.column]-1),
                                 stringr::str_dup("n", 1),
                                 stringr::str_dup("-", (col.n-c(1:col.n)[data.column]))
                               ))
               )
# Unit: milliseconds                                                                                                               
# expr
# data.table::fread()
# readr::read_csv()
#       min        lq      mean    median        uq       max neval
#  54.49732  56.41148  65.77239  62.13816  75.45741  94.98466   100
# 126.07687 151.37307 171.42465 161.03843 173.08058 467.76940   100

## arrow::read_feather()  ----



## vroom ----
# https://github.com/r-lib/vroom 


# Load data in ----
library("fpeek")
sapply(paste0(data.dir,"/",list1[1]), peek_count_lines)
sapply(paste0(data.dir,"/",list1[1]), peek_head)[[1]]

peek1 = peek_head(paste0(data.dir,"/",list1[1]), n = 10, intern = TRUE)[-1]
peek1 = textConnection(peek_head(paste0(data.dir,"/",list1[1]), n = 10, intern = TRUE)[-1])
max(count.fields(peek1, sep = ","))


count.fields(paste0(data.dir,"/",list1[1]), sep = ",", skip = 1185655)[1]
count.fields(paste0(data.dir,"/",list1[1]), sep = ",")[1]


## example ----
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv"
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group"
data.column = 1
skip = 1
header = F

i = 1

file.list = list.files(path = paste0(data.dir), pattern = ".csv$", recursive = TRUE)

peek.n10 = textConnection(peek_head(paste0(data.dir,"/",file.list[i]), n = 10, intern = TRUE)[-1])
col.n = max(count.fields(peek.n10, sep = ","))

file1 = read.csv(paste0(data.dir,"/",file.list[i]), skip = skip , header = header, 
                 colClasses = c(rep("NULL", c(1:col.n)[data.column]-1), 
                                rep("numeric", 1), 
                                rep("NULL", (col.n-c(1:col.n)[data.column]))))



  meta.data = data.frame(x1 = c(10:20), x2 = c(20:30))
  meta.data2 <- meta.data[c(6:9), , drop = FALSE]
  meta.data <- meta.data[intersect(x = meta.file, y = meta.file[-2]), , drop = FALSE]
  (x <- c(sort(sample(1:20, 9)), NA))
  (y <- c(sort(sample(3:23, 7)), NA))
  intersect(x, y, drop = FALSE)
  
  c(1,10,5,2,3)[order(c(1,10,5,2,3))]
  sort(c(1,10,5,2,3))
  
    # CreateSeuratObject.Assay 
    # https://rdrr.io/cran/SeuratObject/src/R/seurat.R
    #   object <- new(
    #                 Class = 'Seurat',
    # .AddMetaData 
    # https://rdrr.io/cran/SeuratObject/src/R/zzz.R
    
q1 = ReadQbone(data.dir)
# Thinning the data at the same time
