library(devtools)

document()

load_all()

load_all(recompile = T)

check()

devtools::check(manual=TRUE)

tinytex::install_tinytex()
tinytex::tlmgr_install("makeindex")
tinytex::tlmgr_update()
tinytex::reinstall_tinytex()
library("tinytex")
devtools::build_manual(pkg = ".")
# readQbone ----
q0 = new("Qbone")
q0 = new("Qbone", version = packageVersion(pkg = "Qbone"))
?dorem

list1 = list(c("123", "345"), "234")
meta.data = data.frame(name = c("a_1", "b_1"), names = c("c", "b"))
rownames(meta.data) = meta.data[,1]

qbonedata = createQboneData(list1, meta.data, sampleid = 2)

qbone1 = createQboneObject(qbonedata, meta.data = meta.data)

qbone2 = createQboneObject(list1, meta.data = meta.data, sampleid = 2)

qbone <- addMetaData(object = qbone, metadata = meta.data)
qbone[["name"]] <- meta.data

cqo2 = addMetaData(qbone, c("c", "d", "e"), col.name = "names")

data.frame(row.names = qbonedata@data)

idents(q2) <- c("c", "d")
idents(q2)

# lassolist ----
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group"
q1 = readQbone(data.dir, groupbyfolder = T)

qa1 = getQboneData(q1, slot = 'data', assay = defaultAssay(q1))


q2 = thinData(q1,prop=0.0001)
q3 = lassolist(q2)
q4 = quantlets(q3)
q5 = lassolist2(q2)

document()
q6 = lassolist2(q2)
all.equal(q6,q4)

## old ----
raw.dataset <- getQboneData(q2, slot = 'data', assay = defaultAssay(q2))
a1 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
a2 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
list1 = runlassolist(raw.dataset[[1]])

## assay.seed ----
set.seed(12345)
q7 = lassolist(q2,verbose = F, parallel = T)
q8 = lassolist(q2, assay.seed = q7@assays[["Lasso"]]@scale.data[["lassolist"]], verbose = F, parallel = T)
all.equal(q8@assays[["Lasso"]]@data,
          q7@assays[["Lasso"]]@data)

all.equal(q7@assays[["Lasso"]]@scale.data[["lassolist"]],
          q8@assays[["Lasso"]]@scale.data[["lassolist"]])

set.seed(12345)
q9 = lassolist(q2, parallel = T)
set.seed(12345)
q10 = lassolist2(q2, parallel = T)
all.equal(q9@assays[["Lasso"]]@data,
          q10@assays[["Lasso"]]@data)

all.equal(q9@assays[["Lasso"]]@scale.data[["lassolist"]],
          q10@assays[["Lasso"]]@scale.data[["lassolist"]])

.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]
q9 = lassolist(q2, parallel = T)
.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]
q10 = lassolist(q2, parallel = T)
all.equal(q9@assays[["Lasso"]]@data,
          q10@assays[["Lasso"]]@data)

all.equal(q9@assays[["Lasso"]]@scale.data[["lassolist"]],
          q10@assays[["Lasso"]]@scale.data[["lassolist"]])

all.equal(q7@assays[["Lasso"]]@scale.data[["lassolist"]],
          q10@assays[["Lasso"]]@scale.data[["lassolist"]])

# For this example, set the random seed
set.seed(423)
runif(3)
#> [1] 0.1089715 0.5973455 0.9726307

# Save the seed
set.seed(423)
oldseed <- .Random.seed
oldseed1 <- .Random.seed
.Random.seed <- oldseed1
oldseed2 <- .Random.seed

all.equal(oldseed,oldseed1)
all.equal(oldseed2,oldseed1)

.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]
runif(3)
.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]

# Get the same random numbers as before, after saving the seed
runif(3)



## performance ----
# library("biglasso")
library("doParallel")
library("foreach")
library("doMC")
library("glmnet")
registerDoMC(3)
do1 = cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = T)

do2 = cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3)


library("microbenchmark")
# lassolist2 is faster (0.1-0.9 second)
document()
q7 = lassolist(q2, parallel = T)
q7.2 = lassolist2(q2, verbose=T, parallel = T)
microbenchmark(lassolist(q2, verbose=F,  parallel = T), lassolist2(q2, verbose=F, parallel = T), times = 10L)
microbenchmark(lassolist2(q2, verbose=F), lassolist2(q2, verbose=T), times = 10L)


microbenchmark(eigenmapmt(CDFBETA), t(CDFBETA), times = 10L)
microbenchmark(generateBetaCDF(alpha = a1, beta = a2, index.p = grid.p), GENERATE_BETA_CDF(alpha = a1, beta = a2, index.p = grid.p), times = 10L)
microbenchmark(centeringFunction(BETASCDF0[ rowth, ], scale = TRUE), centering.function(BETASCDF0[ rowth, ], scale = TRUE), times = 10L)

microbenchmark(cv.biglasso(as.big.matrix(BETA_BASE_TOTAL_2), y, nfolds = 3),cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3), times = 10L)

microbenchmark(cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3),cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = T), times = 10L)

microbenchmark(lassolist2(q2, verbose=T, parallel = T),
               lassolist2(q2, verbose=F, parallel = T),
               times = 10L)

microbenchmark(try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y),
               try(eigenmapmm(smPsi_i_, eigenmapmm(ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps)), eigenmapmtm(smPsi_i_,y)))), # faster in large matrix
               times = 10L)
# Thin = 0.01
# Unit: milliseconds
# expr
# try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*%      t(smPsi_i_) %*% y)
# try(eigenmapmm(smPsi_i_, eigenmapmm(ginv(eigenmapmtm(smPsi_i_,      smPsi_i_), tol = sqrt(.Machine$double.eps)), eigenmapmtm(smPsi_i_,      y))))
# min          lq        mean      median         uq        max neval
# 16106.87076 16250.45060 16388.39732 16281.84498 16430.4827 17075.7220    10
# 25.95431    26.37465    34.45062    27.17549    35.3364    76.6718    10

all.equal(try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y),
          try(eigenmapmm(smPsi_i_, eigenmapmm(ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps)), eigenmapmtm(smPsi_i_,y))))
          )



microbenchmark(try(eigenmapmm(eigenmapmmt(eigenmapmm(smPsi_i_, ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps))), t(smPsi_i_)), y)),
               try(eigenmapmm(smPsi_i_, eigenmapmm(ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps)), eigenmapmtm(smPsi_i_,y)))), # faster
               times = 10L)

assignnew <- function(){
  gitsmp = ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps))
  try1 = try(eigenmapmtm(smPsi_i_,y))
  try2 = try(eigenmapmm(gitsmp, try1))
  try(eigenmapmm(smPsi_i_, try2))}

microbenchmark(assignnew(),
               try(eigenmapmm(smPsi_i_, eigenmapmm(ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps)), eigenmapmtm(smPsi_i_,y)))), # faster
               times = 10L)



microbenchmark(
  locc(
    leaveout.list = lasso_IncidenceVec_i_,
    remain.counts = lasso.counts.fit[[3]],
    remain.basis = lasso.counts.fit[[2]],
    Y.list = raw.dataset,
    maxim = length(p),
    alpha = alpha,
    beta = beta
  ),
  LCCC.FIT_1(
    leaveout.list = lasso_IncidenceVec_i_,
    remain.counts = lasso.counts.fit[[3]],
    remain.basis = lasso.counts.fit[[2]],
    Y.list = raw.dataset,
    maxim = length(p),
    alpha = alpha,
    beta = beta
  ),
               times = 10L)

microbenchmark(
  locc(
    leaveout.list = lasso_IncidenceVec_i_,
    remain.counts = lasso.counts.fit[[3]],
    remain.basis = lasso.counts.fit[[2]],
    Y.list = raw.dataset,
    maxim = length(p),
    alpha = alpha,
    beta = beta
  ),
  locc2(
    leaveout.list = lasso_IncidenceVec_i_,
    remain.counts = lasso.counts.fit[[3]],
    remain.basis = lasso.counts.fit[[2]],
    Y.list = raw.dataset,
    maxim = length(p),
    alpha = alpha,
    beta = beta
  ),
  times = 10L)

lasso.locc <- locc(
  leaveout.list = lasso_IncidenceVec_i_,
  remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]],
  Y.list = raw.dataset,
  maxim = length(p),
  alpha = alpha,
  beta = beta
)
### for each ----
library("parallel")
library("doParallel")

library("foreach")
library("doMC")
library("glmnet")
document()
load_all()
#### doParallel ----
cores = detectCores()
cl <- makeCluster(cores-1)
registerDoParallel(cl)

set.seed(12345)
q7 = lassolist(q2, verbose=T, parallel = T)
set.seed(12345)
q7.2 = lassolist2(q2, verbose=F, parallel = T)
all.equal(q7@assays[["Lasso"]]@data,
          q7.2@assays[["Lasso"]]@data)

stopCluster(cl)

#### doSNOW ----

document()

load_all()
# install()
# remove.packages("Qbone")
library("Qbone")
# detach("package:Qbone", unload=TRUE)

library("parallel")
library("doSNOW")
cores = detectCores()
cl <- makeCluster(cores-1)
registerDoSNOW(cl)
set.seed(12345)
q7 = lassolist(q2, verbose=T, parallel = T)
set.seed(12345)
q7.2 = lassolist2(q2, verbose=T, parallel = T)
set.seed(12345)
q7.3 = lassolist2(q2, verbose=F, parallel = T)


system.time(q7.0 <- lassolist(q2, verbose=F, parallel = T))
system.time(q7.1 <- lassolist(q2, verbose=T, parallel = T))
system.time(q7.2 <- lassolist2(q2, verbose=T, parallel = T))
system.time(q7.3 <- lassolist2(q2, verbose=F, parallel = T))

#### generateBetaCDF is slow ----

### Seeds for parallel----
# https://www.r-bloggers.com/2018/07/%F0%9F%8C%B1-setting-a-seed-in-r-when-using-parallel-simulation/
# https://stackoverflow.com/questions/58631433/how-to-set-seeds-when-using-parallel-package-in-r

set.seed(12345)
q9 = lassolist2(q2, verbose=T, parallel = T)
set.seed(12345)
q10 = lassolist2(q2, verbose=T, parallel = T)
all.equal(q9@assays[["Lasso"]]@data,
          q10@assays[["Lasso"]]@data)

.Random.seed <- q9@assays[["Lasso"]]@scale.data[["lassolist"]]
q9 = lassolist2(q2, verbose=T, parallel = T)
.Random.seed <- q9@assays[["Lasso"]]@scale.data[["lassolist"]]
q10 = lassolist2(q2, verbose=T, parallel = T)
all.equal(q9@assays[["Lasso"]]@data,
          q10@assays[["Lasso"]]@data)


stopCluster(cl)

runlassolist2 <- function(x, parallel = T, ...){
  y <- x
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  CDFBETA <- GENERATE_BETA_CDF(
    # alpha = a1, beta = a2,
    index.p = grid.p, ...)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
  # set.seed(12345)                                                     # || double check
  # lasso_fit <- biglasso(as.big.matrix(BETA_BASE_TOTAL_2), y, ncores = detectCores()-2) # , intercept = TRUE)
  # set.seed(12345)
  # .Random.seed <- assay.seed
  lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
  # set.seed(12345)
  # cvfit.lasso <- cv.biglasso(as.big.matrix(BETA_BASE_TOTAL_2), y, nfolds = 3) #, intercept = TRUE, nfolds = 3)
  # summary(cvfit.lasso)
  if (parallel){
    # message("parallel =", parallel)
    registerDoMC(3)
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = T)
  } else {
    # message("parallel =", parallel)
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = F)
  }
  # set.seed(12345)
  # cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 5, parallel = F)
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  # zeros2 <- as.vector(coef(lasso_fit2, s = cvfit.lasso2$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
  # selects2 <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros2 == FALSE ]
  # gc(verbose = FALSE)
  return(selects)
}

## verbose ----
testit <- function(x = sort(runif(20)), ...)
{
  pb <- txtProgressBar(...)
  for(i in c(0, x, 1)) {Sys.sleep(0.5); setTxtProgressBar(pb, i)}
  Sys.sleep(1)
  close(pb)
}
testit()
testit(runif(10))
testit(style = 3)

library(plyr)
laply(1:100, function(i) {Sys.sleep(0.05); i}, .progress=plyr::progress_text(style = 3) )

total <- 50
pb <- txtProgressBar(min = 0, max = total, style = 3)

lapply(1:total, function(i){
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
})

pb <- txtProgressBar(min = 0, max = length(orig.dataset), style = 3)
new.data <- list()
for (i in 1:length(orig.dataset)){
  bar1 <- runlassolist2(orig.dataset[[i]],
                        alpha = alpha,
                        beta = beta,
                        parallel = parallel)
  new.data <- append(new.data, list(bar1))
  setTxtProgressBar(pb, i)
}
close(pb)

# sort ----
library("parallel")
library("doParallel")
# library("foreach")
library("doMC")
library("glmnet")
document()
load_all()
#### doParallel ----
cores = detectCores()
cl <- makeCluster(cores-1)
registerDoParallel(cl)

q4 = lassolist(q2, verbose = F, parallel = T)
q5 = lassolist2(q2, verbose = T, parallel = T)
stopCluster(cl)

all.equal(q4@assays[["Lasso"]]@data,
          q5@assays[["Lasso"]]@data)

vector = q4@assays[["Lasso"]]@data[[1]]

sort(vector, method = "radix")
sort(vector, method = "quick")

library("microbenchmark")
microbenchmark(sort(vector, method = "radix"), sort(vector, method = "quick"))

# eigen mm ----
one <- rep(1, x.n)
one1 = matrix(rep(1, x.n), ncol = x.n)

x3 = one %*% object.x
x4 = eigenmapmm(one1, object.x)

# idents ----
# sloop::ftype(packageVersion(x = "0.0.0.9000"))
class(package_version(x = '2.99.0'))
idents(cqo2) <- c("a1", "b1")
samples(cqo2)

idents <- factor(x = unlist(x = lapply(
  # X = colnames(x = meta.data),
  X = names(data),
  FUN = ExtractField,
  field = field,
  delim = delim
)))

# PCA ----
p1024 <- signif(seq(0.001, 0.999, length = 1024), 4)
Qy = matrix(round(unlist(lapply(raw.dataset, quantile, probs=p1024, type=6)), 3), 1024)

# Beta ----
n <- length(raw.dataset)
leaveout.list <- vector("list", n) ### generate leave-one-out
for (i in 1:n) {
  lasso_fitIncd <- incidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
  leaveout.list[[i]] <- lasso_fitIncd
}

## agree.cci--  computes concordance correlation index in the version of matrix
##  input
## 	input.yhat: fitted values as the matrix
## 	input.y   : emprical true values as the matrix
##  output: concordance value for the matrtix

agree.cci <- function(input.yhat, input.y) {
  k <- dim(input.yhat)[2]
  values <- rep(NA, k)
  for (i in 1:k) {
    values[i] <- concordance(input.y[, i], input.yhat[, i])
  }
  return(values)
}

## concordance --  computes concordance correlation index
##  input
## 	yhat: fitted values
## 	y   : emprical true values
##  output: concordance value

concordance <- function(yhat, y) {
  2 * cov(yhat, y, use = "complete.obs") / (cov(yhat, yhat, use = "complete.obs") + cov(y, y, use = "complete.obs") + (mean(y, na.rm = TRUE) - mean(yhat, na.rm = TRUE))^2)
}


# other ----

Qy2 = matrix( round(unlist( lapply( raw.dataset,  quantiles_p )  ),3) , 1024 )

all.equal(Qy, Qy2)
object <- new(
  Class = 'Qbone',
  assays = assay.list,
  meta.data =  data.frame(row.names = names(data@data)),
  active.assay = assay,
  active.ident = idents,
  project.name = project,
  version = packageVersion(pkg = 'Qbone'))


mate = c("c", "d")
names(mate) = c("a_1", "b_1")
cqo2 = addMetaData(qbone, mate, col.name = "names")
assays(qbonedata)
assays <- FilterObjects(object = qbone, classes.keep = 'QboneData')
names(getQboneData(qbonedata))

foo1 <- function(verbose =T ){
  for (i in 1:100){
    if (verbose) {
      message("Calculating cluster ", i)
    }
    Sys.sleep(1)
  }

}
for (i in 1:100){
    message("Calculating cluster ", i)
  }

call1 = function(x,y){print(x)}

locc2 <- function(
  leaveout.list,
  remain.counts,
  remain.basis,
  Y.list,
  maxim = 1000,
  ...
) {
  n <- length(leaveout.list)
  active.set <- (remain.basis < maxim)
  if (tail(remain.counts[active.set], 1) == n - 1) {
    remain.counts[ length(remain.counts)] <- n - 2
  }
  feasible.long <- sum(active.set)
  max.long <- max(unlist(lapply(Y.list, length), use.names = FALSE))
  # checks <- matrix(NA, nrow = n, ncol = feasible.long)
  Values <- array(NA, c(max.long, n, feasible.long))

  pb <- txtProgressBar(min = 0, max = length(Y.list), style = 3)
  for (i in 1:n) {
    message("\n Computes the leave-one-out concordance correlation index for ", names(Y.list)[i])
    y <- Y.list[[i]]
    y.long <- length(y)
    grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
    CDFBETA <- generateBetaCDF( # a1, a2,
      grid.p, ...)
    NQ <- qnorm(grid.p, 0, 1)
    BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
    BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
    Psi <- cbind(rep(1, length(grid.p)), BETA_BASE_TOTAL_2)

    # pbj <- txtProgressBar(min = 0, max = length(feasible.long), style = 3)
    for (j in 1:feasible.long) { ###  j = 20

      colum_i_ <- leaveout.list[[i]][[1]]
      obs_i_ <- leaveout.list[[i]][[2]] ## length(IncidenceVec_i_)

      SET1_i_ <- sort(colum_i_[obs_i_ > remain.counts[active.set][j]])
      smPsi_i_ <- Psi[, (SET1_i_) + 1] ## dim(smPsi_i_)

      # Values[(1:y.long), i, j] <- try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y)
      # Values[(1:y.long), i, j] <- try(smPsi_i_ %*% (ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% (t(smPsi_i_) %*% y)))
      # Values[(1:y.long), i, j] <- try(smPsi_i_ %*% (ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps)) %*% eigenmapmtm(smPsi_i_,y)))
      ## dim(Qy)
      # gitsmp = ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps))
      # try1 = try(eigenmapmtm(smPsi_i_,y))
      # try2 = try(eigenmapmm(gitsmp, try1))
      # Values[(1:y.long), i, j] <- try(eigenmapmm(smPsi_i_, try2))
      Values[(1:y.long), i, j] <- try(eigenmapmm(smPsi_i_, eigenmapmm(ginv(eigenmapmtm(smPsi_i_, smPsi_i_), tol = sqrt(.Machine$double.eps)), eigenmapmtm(smPsi_i_,y))))
      # checks[i, j] <- length(SET1_i_)
      # setTxtProgressBar(pbj, i)
    }
    # close(pbj)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  outputs <- list(Values
                  # , checks
  )
}
