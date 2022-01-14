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
getwd()


# matrix ----
A <- matrix(rnorm(1000000), 1000, 1000)
B <- matrix(rnorm(1000000), 1000, 1000)

sourceCpp("example.cpp")
sourceCpp("example2.cpp")

sourceCpp("cmatrix.cpp")

o1 = A%*%B
a1 = armaMatMult(A, B)
e1 = eigenMatMult(A, B)
e2 = eigenMapMatMult(A, B)

o3 = A%*%t(B)
e3 = eigenMapMatMulttrans(A,B)
e4 = eigenMatMulttrans(A,B)

# test result equality----
# https://stackoverflow.com/questions/23032387/how-to-compare-two-matrices-to-see-if-they-are-identical-in-r
identical(o3, e3)
all.equal(o3, e3)

max(abs(o1 - e2))

# test running time ----
library(microbenchmark)

microbenchmark(A%*%B, armaMatMult(A, B), eigenMatMult(A, B), eigenMapMatMult(A, B))
microbenchmark(A%*%B, eigenMatMult(A, B), eigenMapMatMult(A, B))
microbenchmark(t(A)%*%B, crossprod(A, B))
microbenchmark(A%*%t(B), eigenMatMulttrans(A,B), eigenMapMatMulttrans(A, B))
# Unit: milliseconds
# expr                             min        lq      mean    median        uq      max neval
# A %*% t(B)                 593.84899 601.65165 612.49599 606.54585 614.23624 772.8392   100
# eigenMatMulttrans(A, B)     30.91786  35.34698  51.72134  43.61956  57.19245 148.0706   100
# eigenMapMatMulttrans(A, B)  27.32339  28.11837  42.55237  35.25694  49.35140 154.2409   100

# on avarage 612.49599/42.55237 = 14.39393

x <- .Random.seed
