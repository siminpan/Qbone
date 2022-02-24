#' @include zzz.R
#' @include generics.R
#' @include utils.R
#' @include QboneData.R
#' @include QboneObject.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Class definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 readQbone ----
#' Read Qbone from files
#'
#' @param data.dir the folder directory where the raw csv file saved
#' @param groupbyfolder Use parent folder names as sample group information if set to T, default is F
#' @param data.column column number of sample data in the raw csv file
#'
#' @importFrom data.table fread
#' @export
#'
readQbone <- function(
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
  # read file
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
  object <- createQboneObject(
    data = full.data,
    meta.data = meta.data
    )
  return(object)
}

## 2.2 thinData ----
#' Thin the data
#' If the data set is too big and the calculation is too heavy for your computer you may use this function to thin the data to make the calculation faster.
#' Here we use thin the data proportionally on the sorted data. The new data set should reserve the same quantile feature as it was sequentially subsetted from ordered data set.
#'
#' @param object An Qboneobject.
#' @param new.assay.name new assay name assigned to the thined data
#' @param prop proportion to keep from the original data.
#'
#' @concept  Data thin
#' @export
#'
thinData <- function(
  object,
  new.assay.name = "Thin",
  prop = NULL
){
  if (is.null(prop)){
    stop("Proportion(prop) must be provided")
  } else if (prop <= 0 | prop >= 1){
    stop("Proportion(prop) must be between 0 and 1.")
  }
  orig.data = getQboneData(object, slot = 'data', assay = defaultAssay(object))
  new.data = vector("list", length(orig.data))
  names(new.data) <- names(orig.data)
  if (object@assays[[defaultAssay(object)]]@scale.data[["sort"]] == F){
    for (i in 1:length(orig.data)) {
      new.data[[i]] = orig.data[[i]][order(unlist(orig.data[[i]]))[seq(1, length(orig.data[[i]]), round(1/prop))]]
    }
  }else {
    for (i in 1:length(orig.data)) {
      new.data[[i]] = unlist(orig.data[[i]])[seq(1, length(orig.data[[i]]), round(1/prop))]
    }
  }
  new.qbonedata <- createQboneData(new.data,
                                   meta.assays = data.frame(id = names(orig.data)),
                                   sampleid.assays = 1,
                                   assay.name = new.assay.name,
                                   assay.orig = defaultAssay(object))
  new.qbonedata@scale.data <- append(object@assays[[defaultAssay(object)]]@scale.data,list(thin = c(T, prop)))
  object[[new.assay.name]] <- new.qbonedata
  defaultAssay(object) <- new.assay.name
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Qbone-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. R-defined generics ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. S4 methods ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Internal ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
