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
    raw.data = 'list',
    thin.data = 'list',
    meta.data = 'data.frame',
    thin.meta = 'list',
    active.assay = 'character',
    active.ident = 'factor',
    lasso.list = 'list',
    project.name = 'character',
    misc = 'list',
    # version = 'package_version',
    commands = 'list',
    tools = 'list'
  ), 
  prototype = c(
    raw.data = list(),
    thin.data = list(),
    meta.data = data.frame(id = NULL),
    thin.meta = list(), # list(T, "0.1")
    active.assay = NA_character_,
    active.ident = factor(),
    lasso.list = list(),
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

# generic ----
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
    project.name = project,
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

data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group2"
q2 = ReadQbone0(data.dir, groupbyfolder = F)
q3 = ReadQbone0(data.dir, groupbyfolder = T)

# need utils and mthods ----

## Subset method ----
# https://rdrr.io/cran/SeuratObject/man/Seurat-methods.html
# "[[.Seurat"  in seurat.R
# https://rdrr.io/cran/SeuratObject/src/R/seurat.R

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

AddMetaData(q2, c(1:15), col.name = "number") 

# CreateSeuratObject.Assay <- function
# https://rdrr.io/cran/SeuratObject/src/R/seurat.R
# Assay <- setClass( https://rdrr.io/cran/SeuratObject/src/R/assay.R

CreateQboneObject <- function(
  data,
  project = 'QboneProject',
  assay = "Bone",
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  ...
) {
  # if (is.null(assay)){
  #   stop('"assay =" is missing, with no default. Please provide a name of your assay ("Bone", "Image" etc.)')
  # }
  if (missing(x = data)) {
    stop("Must provide 'data'")
  }
  
  object <- new(
    Class = 'Qbone',
    raw.data = full.data,
    meta.data = meta.data,
    project.name = project,
    # version = packageVersion(pkg = 'Qbone'),
  )
  
}

list1 = list(c("123", "345"), "234")
names(list1) <- c("name1", "name2")
list1 = list(q1@raw.data)
names(list1)
i = 2
list1[[i]]
list2 = list(list1 = list1, list2 = list1)
lengths(list2)
is.atomic(q1@raw.data[[1]])

# Read10X 
# https://rdrr.io/cran/Seurat/src/R/preprocessing.R
ReadQbone <- function(
  data.dir,
  groupbyfolder = F, 
  data.column = 1,
  # skip = 1,
  # header = F,
  # project = "QboneProject"
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
      # )
      #V3
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
