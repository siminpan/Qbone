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

## 2.2 trainingData ----
#'
#' @export
#'
trainingData <- function(

){

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
