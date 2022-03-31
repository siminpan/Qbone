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

# lassoList ----
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group"
q1 = readQbone(data.dir, groupbyfolder = T)

qa1 = getQboneData(q1, slot = 'data', assay = defaultAssay(q1))

# save image ----
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.group"
q1 = readQbone(data.dir, groupbyfolder = T)
q2 = thinData(q1,prop=0.0001)
q3 = lassoList(q2)
q4 = preQuantlets(q3)
dxPlot(q4)
dxPlotRev(q4)
q5 = ecQuantlets(q4)
object = q5
q6 = qfrModel(q5, X1 = PX0)
save.image("~/Documents/Qbone.test.Rdata")

document()
q6 = lassolist2(q2)
all.equal(q6,q4)


covM <- function(object, int = T, added.cov = NULL){
  if (int == T){
    int <- rep(1, length(object@meta.data[["group"]]))
  } else {
    int <- rep(0, length(object@meta.data[["group"]]))
  }
  X <- cbind(int)
  table1 <- table(object@meta.data[["group"]])
  for (i in 1:length(names(table1))){
    assign(paste0("g", i), ifelse(object@meta.data[["group"]] == names(table1)[i], 1, 0))
    X <- cbind(X, get(paste0("g", i)))
  }
  colnames(X) = c("int", names(table1))
  X <- as.matrix(X)
  if (!is.null(added.cov)){
    X <- cbind(X, added.cov)
  }
  return(X)
}

X <- covM(object)


## old ----
raw.dataset <- getQboneData(q2, slot = 'data', assay = defaultAssay(q2))
a1 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
a2 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
list1 = runlassolist(raw.dataset[[1]])

## assay.seed ----
set.seed(12345)
q7 = lassoList(q2,verbose = F, parallel = T)
q8 = lassoList(q2, assay.seed = q7@assays[["Lasso"]]@scale.data[["lassoList"]], verbose = F, parallel = T)
all.equal(q8@assays[["Lasso"]]@data,
          q7@assays[["Lasso"]]@data)

all.equal(q7@assays[["Lasso"]]@scale.data[["lassolist"]],
          q8@assays[["Lasso"]]@scale.data[["lassolist"]])

set.seed(12345)
q9 = lassoList(q2, parallel = T)
set.seed(12345)
q10 = lassolist2(q2, parallel = T)
all.equal(q9@assays[["Lasso"]]@data,
          q10@assays[["Lasso"]]@data)

all.equal(q9@assays[["Lasso"]]@scale.data[["lassolist"]],
          q10@assays[["Lasso"]]@scale.data[["lassolist"]])

.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]
q9 = lassoList(q2, parallel = T)
.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]
q10 = lassoList(q2, parallel = T)
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
q7 = lassoList(q2, parallel = T)
q7.2 = lassolist2(q2, verbose=T, parallel = T)
microbenchmark(lassoList(q2, verbose=F,  parallel = T), lassolist2(q2, verbose=F, parallel = T), times = 10L)
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

microbenchmark(length((table(object@meta.data[["group"]]))), length(unique(object@meta.data[["group"]])))
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
q7 = lassoList(q2, verbose=T, parallel = T)
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
q7 = lassoList(q2, verbose=T, parallel = T)
set.seed(12345)
q7.2 = lassolist2(q2, verbose=T, parallel = T)
set.seed(12345)
q7.3 = lassolist2(q2, verbose=F, parallel = T)


system.time(q7.0 <- lassoList(q2, verbose=F, parallel = T))
system.time(q7.1 <- lassoList(q2, verbose=T, parallel = T))
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

q4 = lassoList(q2, verbose = F, parallel = T)
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

# F3A plotting ----
library("ggplot2")


quantlet.set <- LCCC.PLOT_1(
  plength = length(p), # p = p1024
  Values = lccc.lasso[[1]], # lasso.locc
  # checking defaultAssay(object):
  Y.list = raw.dataset, # raw.dataset <- getQboneData(object, slot = 'data', assay = object@assays[[defaultAssay(object)]]@assay.orig)
  output1 = lasso.counts.fit[[1]], # object@assays[[defaultAssay(object)]]@scale.data[["basis.columns"]]
  output2 = lasso.counts.fit[[2]], # object@assays[[defaultAssay(object)]]@scale.data[["remain.basis"]]
  output3 = lasso.counts.fit[[3]], # object@assays[[defaultAssay(object)]]@scale.data[["remain.counts"]]
  cutoff = 0.990 # Near-lossless values.
  )


texts1 <- expression(paste(bar(rho), " vs K"))
texts2 <- expression(paste(rho^0, " vs K"))

cex.lab <- 0.9 ### -> Size of labels for both x and y axis!
cex.axis <- 0.65 ### -> Size of coordinates for both x and y axis!
cex.main <- 0.7 ### -> Size of main topic ( USELESS since we will cut title )
cex.in <- 0.5 ### -> Size of all in boxplot (Ex: outlier )

yaxis.at <- c(-2, seq(0.85, 1.0, by = 0.05), 1.5)
yaxis.lab <- c("", "0.85", "0.90", "0.95", "1.00", "")
ylim <- c(0.84, 1.0)

numbasis <- sort(unique(c(lasso.x)), decreasing = TRUE)
plotxaxis <- seq(length(numbasis))

xaxis.at <- c(plotxaxis[1] - 15, plotxaxis, tail(plotxaxis, 1) + 15)
xaxis.lab <- c("", numbasis, "")
xaxis.lab1 <- c("", lasso.x1, "")


## par( mfrow=c(1,1), mar=c(4.5, 2.5, 3,2 ))

plot(0, type = "o", xlab = "", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim, xlim = c(0, length(numbasis) + 1))
points(plotxaxis, lasso.Chary1_i_, col = "red", pch = 19)
lines(plotxaxis, lasso.Chary1_i_, col = "red")

axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)

axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)
axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)
title(texts2)
## title(texts2 )

points(plotxaxis, lasso.Chary_i_, col = "blue", pch = 19)
lines(plotxaxis, lasso.Chary_i_, col = "blue")

axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)
axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)
axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)
title(texts1)


numbasis <- sort(unique(c(lasso.x)), decreasing = TRUE)
plotxaxis <- seq(length(numbasis))

xaxis.at <- c(plotxaxis[1] - 15, plotxaxis, tail(plotxaxis, 1) + 15)

plotdata1 = data.frame(plotxaxis = plotxaxis, lasso.Chary_i_ = lasso.Chary_i_, lasso.Chary1_i_ = lasso.Chary1_i_)

p1 <- ggplot(plotdata1, aes(x=plotxaxis, y=c(lasso.Chary_i_, lasso.Chary1_i_))) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0.85, 1))

p1

plotdata2 = data.frame(x = c(plotxaxis, plotxaxis), x2 = c(numbasis, numbasis), y = c(lasso.Chary_i_, lasso.Chary1_i_), group = c(rep("mean", length(lasso.Chary_i_)), rep("min", length(lasso.Chary1_i_))))


plotdata1 <- loccplotdata(object)
plotdata2 <- plotdata1[[3]]
plotdata3 <- plotdata2[which(plotdata2$y >= 0.85),]
lasso.Chary_i_ = plotdata2$y[plotdata2$group == "mean"]
lasso.Chary1_i_ = plotdata2$y[plotdata2$group == "min"]

p2 <- ggplot(plotdata2, aes(x=x, y=y, color= group)) +
  geom_point(data = plotdata3, aes(x=x, y=y, color= group)) +
  geom_line(data = plotdata3, aes(x=x, y=y, color= group)) +
  scale_y_continuous(limits = c(0.8, 1.0)) +
  scale_x_continuous(limits = c(-1, max(plotdata2$x)+0.5))+
  scale_color_manual(
    labels = c(expression(paste(bar(rho))), expression(paste(rho^0))),
    values = c("blue", "red")
  ) +
  labs(title = "Find a sparse yet near-lossless basis set",
       y = "Losslessness",
       color = NULL
  ) +
  geom_vline(xintercept = min(plotdata2$x[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001]),
             linetype="dotted",
             color = "black",
             size=0.5) +
  geom_hline(yintercept = plotdata2$cutoff[1],
             linetype="dotted",
             color = "black",
             size=0.5) +
  geom_label(aes(min(x[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001]),
                 0.855,
                 label=min(x[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001])
  )
  ) +
  geom_label(aes(unique(plotdata2$x)[2],
                 plotdata2$cutoff[1],
                 label="cutoff")
  )

p2
suppressWarnings(p2)
suppressWarnings(print(p2))

p2 +
  annotate(geom = "text", x = c(0,unique(plotdata2$x)), y = 0.84, label = c("C",unique(plotdata2$x)), size = 3) +
  annotate(geom = "text", x = c(0,unique(plotdata2$x)), y = 0.83, label = c(expression(K[C]),unique(plotdata2$x2)), size = 3, angle = 45) +
  coord_cartesian(ylim = c(0.85, 1.01), xlim = c(0, max(plotdata2$x)+0.5), expand = FALSE, clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


p2.5 = dxPlot(object)
p2.5

p3 <- ggplot(plotdata2, aes(x=x, y=y, color= group)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0.8, 1.0)) +
  scale_x_reverse(limits = c(max(plotdata2$x)+0.5, 0)) +
  scale_color_manual(
    labels = c(expression(paste(bar(rho))), expression(paste(rho^0))),
    values = c("blue", "red")
  ) +
  labs(title = "Find a sparse yet near-lossless basis set",
       y = "Losslessness",
       color = NULL
  ) +
  geom_vline(xintercept = min(plotdata2$x[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001]),
             linetype="dotted",
             color = "black",
             size=0.5) +
  geom_hline(yintercept = cutoff,
             linetype="dotted",
             color = "black",
             size=0.5) +
  geom_label(aes(min(x[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001]),
                 0.855,
                 label=max(x2[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001])
  )
  ) +
  geom_label(aes(unique(plotdata2$x)[2],
                 cutoff,
                 label="cutoff")
  ) +
  annotate(geom = "text", x = c(rev(unique(plotdata2$x)), 0), y = 0.84, label = c(expression(K[C]),rev(unique(plotdata2$x2))), size = 3, angle = 45) +
  coord_cartesian(ylim = c(0.85, 1.01), xlim = c(max(plotdata2$x)+0.5, 0), expand = FALSE, clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

p3

## off limit testing ----
z <- rnorm(1000, mean = 0, sd = 1)
dens <- density(z)
data <- tibble(x = dens$x, y = dens$y) %>%
  mutate(variable = case_when(
    (x >= -2 & x <= 0) ~ "On",
    (x >= 0.2 & x <= 1) ~ "Off",
    TRUE ~ NA_character_))

ggplot(data, aes(x, y)) + geom_line() +
  geom_area(data = filter(data, variable == 'On'), fill = 'grey') +
  geom_area(data = filter(data, variable == 'Off'), fill = 'light blue') +
  geom_text(aes(x = 7, y = 0.15, label = "test")) +
  theme(plot.margin = unit(c(1,7,1,1), "lines")) +
  scale_y_continuous(limits = c(0.1, 0.2),
                     expand = c(0, 0),
                     oob = scales::oob_squish) +
  scale_x_continuous(limits = c(-5, 5),
                     oob = scales::oob_keep) +
  coord_cartesian(clip = "off")

document()

dxPlot(object)
dxPlotRev(object)

p2 = dxPlot(object)
object@graphs <- list(p2)

# FIg 6 ----
mcmcinfer_object = mcmc_infer
p = p1024
edit=10
opt = 1
n.sup = 100
xdomain <- seq(xranges1[1], xranges1[2], length.out = n.sup)

plot( 0, type="n",    ylim=c(0,11), xlim=c(-0.2,0.3)  )
lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,1] , col="black", lty=1 , lwd=1) # NT
lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,3] , col="red", lty=1 , lwd=1)  # DRB
lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,5] , col="blue", lty=1 , lwd=1) # DKKMo
lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,7] , col="orange", lty=1 , lwd=1) # DKKMoDRB
title( "Predicted densities" , cex=1.5)

legend( 0.1, 11,
        c("NT", "Combination", "DKKMo", "DRB"),
        lty= c(1,1,1,1)  ,
        col =c("black","red","blue" ,"orange") ,
        cex = 1 , bty = "n", ncol=1)

# LCCC ----
lasso.locc <- locc(
  leaveout.list = lasso_IncidenceVec_i_,
  remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]],
  Y.list = raw.dataset,
  maxim = length(p),
  alpha = alpha,
  beta = beta
)

lasso.locc2 <- LCCC.FIT_1(
  leaveout.list = lasso_IncidenceVec_i_,
  remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]],
  Y.list = raw.dataset,
  maxim = length(p),
  alpha = alpha,
  beta = beta
)

all.equal(lasso.locc[[1]], lasso.locc2[[1]])

############# We choose # of basis ###############################################################################

lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x.idx <- seq(length(lasso.list1))[ (lasso.counts.fit[[2]] < 1000) ]
list.order <- lasso.counts.fit[[1]]
unlist(lapply(list.order, length))

# lasso.counts.fit[[1]]:
object@assays[["Quantiles"]]@scale.data[["basis.columns"]]

object@assays[["Quantiles"]]@scale.data[["basis.columns"]][[9]]
length(object@assays[["Quantiles"]]@scale.data[["basis.columns"]][[9]])

be <- 3
REDUCED_BASE9 <- BETA_BASE_TOTAL_2[, list.order[[be + 3 ] ]] # our choice in paper
REDUCED_BASE21 <- BETA_BASE_TOTAL_2[, list.order[[be + 13 ] ]] # Normal case in paper


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

# real bone project results ----
load("/home/span/Documents/MOSJ-3DCT/data/05.6_3Dpoints/quafunreg_raw.data.RData")
re1 = createQboneObject(raw.dataset)

names(re1@assays[["Bone"]]@data) <- c(paste0("DkkMo_", seq(10)),
                                      paste0("DkkMoDRB_", seq(9)),
                                      paste0("DRB_", seq(10)),
                                      paste0("NT_", seq(10))
                                      )
re1@assays[["Bone"]]@meta.assays <- data.frame(id = names(re1@assays[["Bone"]]@data), row.names = c(1:39))

re1 <- addMetaData(object = re1, metadata = re1@assays[["Bone"]]@meta.assays[["id"]]
                   , col = "id")
re1 <- addMetaData(object = re1, metadata = c(rep("DkkMo", c(10)),
                                              rep("DkkMoDRB", c(9)),
                                              rep("DRB", c(10)),
                                              rep("NT", c(10))
                                              ),
                   col = "group")

# no name if from list
library("devtools")
document()
re2 = thinData(re1,prop=0.01)
re3 = lassoList(re2)
re4 = preQuantlets(re3)
# object = re4
# dxPlot(re4)
# dxPlotRev(re4)
re5 = ecQuantlets(re4)
# object = re5
re6 = qfrModel(re5, X1 = PX0)
save.image(file = "/home/span/Documents/MOSJ-3DCT/data/05.6_3Dpoints/test.Qbone.RData")
# save(list=c("re1", "re2", "re3", "re4", "re5"), file = "/home/span/Documents/MOSJ-3DCT/data/05.6_3Dpoints/test.Qbone.RData")

load("/home/span/Documents/MOSJ-3DCT/data/05.6_3Dpoints/test.Qbone.RData")

load("/home/span/Documents/MOSJ-3DCT/data/05.6_3Dpoints/quafunreg_5.RData")


new.assay.name = "Q.F.Regression"
data.assay = defaultAssay(object)
X = NULL
delta2 = 0.95
H = NULL
pct.range = c(0.05, 0.95)
assay.seed2 = .Random.seed

# NT vs DkkMo, NT vs DkkMoDRB
# NT,DRB,Dkk,inter
mm_1 <- c(1,0.25,0,0.25)
mm_2 <- c(0,0.25,1,0.25)
mm_3 <- c(1,0.25,0.25,0)
mm_4 <- c(0,0.25,0.25,1)
mm_5 <- c(0.25,1,0.25,0)
mm_6 <- c(0.25,0,0.25,1)
mm_7 <- c(1,1,0.25,0.25)
mm_8 <- c(0.25,0.25,1,1)

PX0 <- rbind(mm_1,mm_2,mm_3,mm_4,mm_5,mm_6,mm_7,mm_8)
X1 <- PX0

# fit q4 in quafunreg.R ----
lasso.list = getQboneData(object, slot = 'data', assay = "Lasso.list")
raw.dataset = getQboneData(object, slot = 'data', assay = "Thin")
raw.dataset2 = getQboneData(object, slot = 'data', assay = "Thin")

all.equal(Values, locc)
all.equal(raw.dataset, raw.dataset2)
all.equal(lasso.locc, lasso.values)

# test hasArg(x1) in function ----

test1 <- function(...) {if(hasArg(x1)){print("ex")}else{F}}
test1(x1= NULL)
x
