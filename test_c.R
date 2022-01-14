# http://adv-r.had.co.nz/Rcpp.html
# https://jbhender.github.io/Stats506/F17/Using_C_Cpp.html
# https://support.rstudio.com/hc/en-us/articles/200486088-Using-Rcpp-with-the-RStudio-IDE
# https://cran.rstudio.com/web/packages/Rcpp/vignettes/Rcpp-introduction.pdf
thisFile()
path0 = "/home/span/Documents/MADDNESS"
setwd(path0)

library("Rcpp")

cppFunction(
'int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}'
)
# add works like a regular R function
add

add(1, 2, 3)

# out juts 1 ----
one <- function() 1L

cppFunction(
'int oneC() {
  return 1;
}'
)

# Scalar input, scalar output ----
signR <- function(x) {
  if (x > 0) {
    1
  } else if (x == 0) {
    0
  } else {
    -1
  }
}

cppFunction(
'int signC(int x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}'
)

# MADDNESS ----
library("Rcpp")
# check AVX2 https://stackoverflow.com/questions/28939652/how-to-detect-sse-sse2-avx-avx2-avx-512-avx-128-fma-kcvi-availability-at-compile
path0 = "/home/span/Documents/MADDNESS/mini"
setwd(path0)
# /external
# /home/span/Documents/MADDNESS/mini/src/quantize/utils/eigen_utils.hpp
# /utils
# /home/span/Documents/MADDNESS/mini/src/quantize/bolt.cpp
sourceCpp(file = paste0(path0, "/src/quantize/mithral.cpp"))


path1 = "/home/span/Documents/MADDNESS/mini_2"
setwd(path1)
sourceCpp(file = paste0(path1, "/src/quantize/mithral.cpp"))
# https://stackoverflow.com/questions/13995266/using-3rd-party-header-files-with-rcpp
# https://knausb.github.io/2017/08/header-files-in-rcpp/
#!!!
# https://coderedirect.com/questions/436012/compilation-error-using-rcpp-with-typedef

m1 = matrix(rep(2,4),2)
m2 = matrix(rep(3,4),2)

mithral_encode(m1, m2)

# How to Set C Flag in R----
# https://stackoverflow.com/questions/21341106/overriding-system-defaults-for-c-compilation-flags-from-r
# https://stackoverflow.com/questions/32586455/how-to-change-and-set-rcpp-compile-arguments
#### https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Configure-and-cleanup

# mkdir ~/.R; touch ~/.R/Makevars
# Makevars
## CXX = g++
## CXXSTD = -std=c++11
## CXXFLAGS = -march=native

# https://stackoverflow.com/questions/58239904/building-an-r-package-which-uses-a-c-code-library

# Difference between CFLAGS and CXXFLAGS ----
# https://stackoverflow.com/questions/5541946/cflags-ccflags-cxxflags-what-exactly-do-these-variables-control

## How to define __AVX2__ without flag ----


# Package Devlopment ----
# https://r-pkgs.org/whole-game.html
library("devtools")
create_package("~/Documents/rpackage/regexcite")
