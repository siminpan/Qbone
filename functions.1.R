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

# test class ----
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

qbone <- setClass(
  Class = 'qbone',
  slots = c(
    raw.data = 'list',
    meta.data = 'data.frame',
    lasso.list = 'list',
    # active.assay = 'character',
    # active.ident = 'factor',
    # graphs = 'list',
    # neighbors = 'list',
    # reductions = 'list',
    # images = 'list',
    project.name = 'character',
    misc = 'list',
    version = 'package_version',
    commands = 'list',
    tools = 'list'
  )
)

# Load data in ----
library("fpeek")
sapply(paste0(data.dir,"/",list1[1]), peek_count_lines)
sapply(paste0(data.dir,"/",list1[1]), peek_head)[[1]]

peek1 = peek_head(paste0(data.dir,"/",list1[1]), n = 10, intern = TRUE)[-1]
peek1 = textConnection(peek_head(paste0(data.dir,"/",list1[1]), n = 10, intern = TRUE)[-1])
max(count.fields(peek1, sep = ","))


count.fields(paste0(data.dir,"/",list1[1]), sep = ",", skip = 1185655)[1]
count.fields(paste0(data.dir,"/",list1[1]), sep = ",")[1]

library("microbenchmark")

microbenchmark(count.fields(paste0(data.dir,"/",list1[1]), sep = ",", skip = 1185655)[1],
               count.fields(paste0(data.dir,"/",list1[1]), sep = ",")[1])
############
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


library("fpeek")
ReadQbone <- function(
  data.dir,
  groupbyfolder = F, 
  data.column = 1,
  skip = 1,
  header = F
  # gene.column = 2,
  # cell.column = 1,
  # unique.features = TRUE,
  # strip.suffix = FALSE
) {
  full.data <- list()
  meta.file <- c()
  if (groupbyfolder == T){
    meta.group <- c()
  }
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
  # read file
    for (i in 1:length(file.list)){
      peek.n10 = textConnection(peek_head(paste0(data.dir,"/",file.list[i]), n = 10, intern = TRUE)[-1])
      col.n = max(count.fields(peek.n10, sep = ","))
      
      file1 =  drop(as.matrix(drop(as.matrix(
        read.csv(paste0(data.dir,"/",file.list[i]),
                 skip = skip , header = header,
                 colClasses = c(rep("NULL", c(1:col.n)[data.column]-1),
                                rep("numeric", 1),
                                rep("NULL", (col.n-c(1:col.n)[data.column]))))
                        ))))
      full.data[[i]] <- file1
      meta.file <- c(meta.file, gsub(".csv", "", basename(file.list[i])))
      if (groupbyfolder == T){
        meta.group <- c(meta.group, dirname(file.list[i]))
      }
    }
    
    # if groupbyfolder ? dirname(file.list)
    # file name gsub(".csv", "", basename(file.list[i]))
    
    
    
    # CreateSeuratObject.Assay 
    # https://rdrr.io/cran/SeuratObject/src/R/seurat.R
    #   object <- new(
    #                 Class = 'Seurat',
    
    # barcode.loc <- file.path(run, 'barcodes.tsv')
    # gene.loc <- file.path(run, 'genes.tsv')
    # features.loc <- file.path(run, 'features.tsv.gz')
    # matrix.loc <- file.path(run, 'matrix.mtx')
    # # Flag to indicate if this data is from CellRanger >= 3.0
    # pre_ver_3 <- file.exists(gene.loc)
    # if (!pre_ver_3) {
    #   addgz <- function(s) {
    #     return(paste0(s, ".gz"))
    #   }
    #   barcode.loc <- addgz(s = barcode.loc)
    #   matrix.loc <- addgz(s = matrix.loc)
    # }
    # if (!file.exists(barcode.loc)) {
    #   stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    # }
    # if (!pre_ver_3 && !file.exists(features.loc) ) {
    #   stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    # }
    # if (!file.exists(matrix.loc)) {
    #   stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    # }
    data <- readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

# Thinning the data at the same time
