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
# https://rdrr.io/cran/SeuratObject/src/R/seurat.R
Qbone <- setClass(
  Class = 'Qbone',
  slots = c(
    raw.data = 'list',
    meta.data = 'data.frame',
    lasso.list = 'list',
    project.name = 'character',
    misc = 'list',
    # version = 'package_version',
    commands = 'list',
    tools = 'list'
  )
)

# Subset method ----
# https://rdrr.io/cran/SeuratObject/man/Seurat-methods.html

# Read Qbone ----
library("fpeek")
ReadQbone <- function(
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
    peek.n10 = textConnection(peek_head(paste0(data.dir,"/",file.list[i]), n = 10, intern = TRUE)[-1])
    col.n = max(count.fields(peek.n10, sep = ","))
    
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

# if groupbyfolder ? dirname(file.list)
# file name gsub(".csv", "", basename(file.list[i]))
# data raw.data thin.data
# Ident

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv"
q1 = ReadQbone(data.dir)

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group2"
q2 = ReadQbone(data.dir, groupbyfolder = F)
q3 = ReadQbone(data.dir, groupbyfolder = T)

## need utils and mthods ----
.AddMetaData <- function(object, metadata, col.name = NULL) {
  if (is.null(x = col.name) && is.atomic(x = metadata)) {
    stop("'col.name' must be provided for atomic metadata types (eg. vectors)")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name # %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    stop("No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  return(object)
}
rlang::`%||%`
AddMetaData <- .AddMetaData

AddMetaData(q2, c(1:15), col.name = "number") 

ReadQbone2 <- function(
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
    peek.n10 = textConnection(peek_head(paste0(data.dir,"/",file.list[i]), n = 10, intern = TRUE)[-1])
    col.n = max(count.fields(peek.n10, sep = ","))
    
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
    meta.data = data.frame(row.names = meta.file),
    lasso.list = list(),
    project.name = project,
    # version = packageVersion(pkg = 'Qbone'),
    misc = list(),
    commands = list(),
    tools = list()
  )
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv"
q1 = ReadQbone2(data.dir)

q1[['id']]

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
