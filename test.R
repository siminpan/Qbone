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

.onUnload1 <- function (libpath) {
  library.dynam.unload("regexcite", libpath)
}

.onUnload2 <- function (libpath) {
  library.dynam.unload("Qbone", libpath)
}

.onUnload1("/home/span/Documents/rpackage/regexcite")
.onUnload2("/home/span/Documents/Qbone")

# matrix ----
A <- matrix(rnorm(1000000), 1000, 1000)
B <- matrix(rnorm(1000000), 1000, 1000)

sourceCpp("example.cpp")
sourceCpp("example2.cpp")

# // references  https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r

sourceCpp("~/Documents/rpackage/bonedeform/src/cmatrix.cpp")

o1 = A%*%B
a1 = armaMatMult(A, B)
e1 = eigenMatMult(A, B)
e2 = eigenMapMatMult(A, B)

o3 = A%*%t(B)
e3 = eigenMapMatMulttrans(A,B)
e4 = eigenMatMulttrans(A,B)

# test result equality----
# https://stackoverflow.com/questions/23032387/how-to-compare-two-matrices-to-see-if-they-are-identical-in-r
identical(o1, e1)
all.equal(o1, e1)

max(abs(o1 - e2))

# test running time ----
library("microbenchmark")

microbenchmark(t(A), eigenmt(A), eigenmapmt(A), times = 10L)
microbenchmark(t(A)%*%B,  eigenmapmm(t(A), B), eigenmapmtm(A, B), times = 10L)
microbenchmark(A%*%B, eigenmm(A, B), eigenmapmm(A, B), times = 10L)

microbenchmark(A%*%B, armaMatMult(A, B), eigenMatMult(A, B), eigenMapMatMult(A, B))

microbenchmark(eigenMapMatMult(A, B), eigenmapmm(A, B), times = 10L)

eigenMapMatMult(A, B)
eigenmapmm(A, B)

microbenchmark(A%*%B, eigenMapMatMult(A, B), .Call('_Qbone_eigenmapmm', PACKAGE = 'Qbone', A, B), times = 10L)


microbenchmark(A%*%B, eigenMatMult(A, B), eigenMapMatMult(A, B), eigenmm(A, B), eigenmapmm(A, B), times = 10L)

microbenchmark(t(A)%*%B, crossprod(A, B))
microbenchmark(A%*%t(B), eigenMatMulttrans(A,B), eigenMapMatMulttrans(A, B))
# Unit: milliseconds
# expr                             min        lq      mean    median        uq      max neval
# A %*% t(B)                 593.84899 601.65165 612.49599 606.54585 614.23624 772.8392   100
# eigenMatMulttrans(A, B)     30.91786  35.34698  51.72134  43.61956  57.19245 148.0706   100
# eigenMapMatMulttrans(A, B)  27.32339  28.11837  42.55237  35.25694  49.35140 154.2409   100

# on avarage 612.49599/42.55237 = 14.39393

bar1 = data.frame(expr = NULL, time = NULL, dim = NULL)
for (i in 10^seq(1,4)){
  set.seed(12345)
  A <- matrix(rnorm(i*i), i, i)
  set.seed(54321)
  B <- matrix(rnorm(i*i), i, i)
  foo1 = microbenchmark(A%*%B, eigenMapMatMult(A, B))
  foo1$dim = i
  bar1 = rbind(bar1, foo1)
}

save(bar1, file = "./test_running_time.1.RData")

## plot running time ----
library("ggplot2")

bar1$min = bar1$time/(10^9)

p1 = ggplot(bar1[which(bar1$dim < 10000),], aes(x=dim, y=min, color=expr)) +
  geom_point() +
  geom_smooth() +
  labs(title="Plot of Computational Time and Data Complexity",
       x ="Data Complexity \n (Matrix Dimensions)", y = "Process Time (seconds)") +
  scale_color_manual(name = "Algorithm", labels = c("pre-optimization", "optimized"),  values = c("red", "blue")) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), labels=c("0", "250 \n x 250", "500 \n x 500", "750 \n x 750", "1000 \n x 1000"))

p1

autoplot(bar1[which(bar1$dim < 10000),c(1,2)])

p2 = ggplot(bar1[,], aes(x=dim, y=min, color=expr)) +
  geom_point() +
  geom_smooth() +
  labs(title="Plot of Computational Time and Data Complexity",
       x ="Data Complexity \n (Matrix Dimensions)", y = "Process Time (seconds)") +
  scale_color_manual(name = "Algorithm", labels = c("pre-optimization", "optimized"),  values = c("red", "blue")) +
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), labels=c("0", "2500 \n x 2500", "5000 \n x 5000", "7500 \n x 7500", "10000 \n x 10000"))

p2

# nameing of package ----
library(available)

available("Q.regressR")

# Skelequant, Q-bone, iQ-Oss, SkeletoR, Q-regressR
