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

## 2.1 lassolist ----
#' Run lasso list to build basis function
#'
#'
#'
#' @param object A Qboneobject
#' @param verbose Print a progress bar
#' @param new.assay.name New assay name assigned to the lassolist data
#' @param alpha Vector containing sequence of beta parameter for internal function \code{generateBetaCDF()}
#' @param beta  Vector containing sequence of beta parameter for internal function \code{generateBetaCDF()}
#' @param data.assay It is the name of the assay whose data will be used to compute the lasso list. Default is the data from the defaultAssay(object).
#' @param assay.seed assay information to add into the QboneData object scale.data. The default of \code{lassolist()} will save the random seed for the run. Use \code{.Random.seed <-  object@assays[["Lasso"]]@scale.data[["lassolist"]]} before run \code{lassolist()} for the same results.
#' @param parallel If TRUE, use parallel foreach to fit each fold in \code{cv.glmnet()}. Default use \code{registerDoMC()} to register parallel.
#' @param ... Arguments passed to other methods
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
lassolist <- function(
  object,
  verbose = TRUE,
  data.assay = defaultAssay(object),
  new.assay.name = "Lasso",
  alpha = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  beta = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  assay.seed = .Random.seed,
  parallel = T,
  ...
){
  if (data.assay == "Lasso"){
    stop('The default assay for this Qbone object is already "Lasso".')
  }
  .Random.seed <- assay.seed
  orig.dataset <- getQboneData(object, slot = 'data', assay = data.assay)
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(orig.dataset), style = 3)
    new.data <- list()
    for (i in 1:length(orig.dataset)){
      bar1 <- runlassolist(orig.dataset[[i]],
                           alpha = alpha,
                           beta = beta,
                           parallel = parallel, ...)
      new.data <- append(new.data, list(bar1))
      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    new.data <- lapply(
      X = orig.dataset,
      FUN = runlassolist,
      alpha = alpha,
      beta = beta,
      parallel = parallel
    )
  }
  new.qbonedata <- createQboneData(new.data,
                                   meta.assays = data.frame(id = names(orig.dataset)),
                                   sampleid.assays = 1,
                                   assay.name = new.assay.name,
                                   assay.orig = data.assay)
  new.qbonedata@scale.data <- append(object@assays[[data.assay]]@scale.data,list(lassolist = c(assay.seed)))
  object[[new.assay.name]] <- new.qbonedata
  defaultAssay(object) <- new.assay.name
  return(object)
}

## 2.2 commonBasis ----
#' Compute common basis
#'
#' @param object A Qboneobject
#' @param new.assay.name new assay name assigned to the commonBasis data
#'
#' @export
#'
commonBasis <- function(
  object,
  new.assay.name = "commonBasis",
  ...
){
  if(defaultAssay(object) != "Lasso"){
    warning('The default assay is not "Lasso" please double the defaultAssay() of this Qbone object. This step should be run on results of lassolist().')
  }
  message('Will compute the common basis based on "', defaultAssay(object), '" results from "', object@assays[["Lasso"]]@assay.orig, '" data.')
  raw.dataset <- getQboneData(object, slot = 'data', assay = object@assays[[defaultAssay(object)]]@assay.orig
)
  lasso.dataset <- getQboneData(object, slot = 'data', assay = defaultAssay(object))
  lasso.list1 <- lapply(lasso.dataset, catNorm) ## 1 term fix
  lasso.nonzero.obs1 <- lapply(lasso.list1, replist) ## D_i, generate 1 vector
  lasso.counts.fit <- countBasis(lasso.list1, lasso.nonzero.obs1)
  n <- length(raw.dataset)
  leaveout.list <- vector("list", n) ### generate leave-one-out
  for (i in 1:n) {
    lasso_fitIncd <- incidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
    leaveout.list[[i]] <- lasso_fitIncd
  }
  p1024 <- signif(seq(0.001, 0.999, length = 1024), 4)
  Qy = matrix(round(unlist(lapply(raw.dataset, quantile, probs=p1024, type=6)), 3), 1024)
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

## 6.1 generateBetaCDF ----
#' creates beta basis functions
#'
#'
#' @param alpha vector containing sequence of beta parameter
#' @param beta  vector containing sequence of beta parameter
#' @param index.p probability grids on (0,1)
#'
#' @return matrix containing # of beta parameters times the length of index.p
#'
#' @keywords internal
#'
#' @noRd
#'
generateBetaCDF <- function(
  alpha,
  beta,
  index.p,
  ...
  ) {
  # n1 <- as.character(alpha)
  # n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(beta) * length(alpha))
  for (i in 1:length(alpha)) { ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(beta)) {
      rowth <- (j - 1) * length(beta) + i
      BETASCDF0[rowth, ] <- pbeta(index.p, alpha[i], beta[j])
      BETASCDF0[rowth, ] <- round(centeringFunction(BETASCDF0[rowth, ], scale = TRUE), 7)
    }
  }
  ## name.mat = outer( paste0( "(", n1 , se="") ,paste0(",", n2, ")" , se="")  ,  FUN=paste ,sep="")
  ## matrix.rownames = as.vector(  name.mat )
  ## outputs= list( BETASCDF0 ,  name.mat , matrix.rownames )
  outputs <- BETASCDF0
  return(outputs)
}

## 6.2 centeringFunction ----
#' computes centering and scaling column by column
#'
#' @param raw.x any matrix with n times p
#' @param scale FALSE or TRUE (doing only centering if scale=FALSE )
#'
#' @return the matrix centerted or/and scaled
#'
#' @keywords internal
#'
#' @noRd
#'
centeringFunction <- function(
  raw.x,
  scale = FALSE
  ) {
  object.x <- as.matrix(raw.x)
  x.n <- dim(object.x)[1]
  one <- rep(1, x.n) # matrix(rep(1, x.n), ncol = x.n)
  meanx <- drop(eigenmapmm(matrix(rep(1, x.n), ncol = x.n), object.x)) / x.n
  mean1 <- as.matrix(meanx)
  n1 <- dim(mean1)[1]
  jj <- eigenmapmm(one,t(rep(1, dim(object.x)[2]))) # rep(1, x.n) %*% t(rep(1, dim(object.x)[2])) # one %*% two # eigenmapmm(one, two)
  mean.mat <- diag(meanx, nrow = n1, ncol = n1)
  cen.x <- object.x - eigenmapmm(jj, mean.mat)
  if (scale == FALSE) {
    return(cen.x)
  }
  if (scale == TRUE) {
    normx <- sqrt(drop(eigenmapmm(matrix(rep(1, x.n), ncol = x.n), (cen.x^2))))
    normx1 <- as.matrix(normx)
    s1 <- dim(normx1)[1]
    scale.mat <- diag(1 / normx, nrow = s1, ncol = s1)
    cen.scale.x <- eigenmapmm(cen.x, scale.mat)
    return(cen.scale.x)
  }
}

## 6.3 runlasso ----
#' Run lasso for list of data
#'
#' @param x any list of data
#' @param parallel If TRUE, use parallel foreach to fit each fold in \code{cv.glmnet(..., nfolds = 3)}. Default use \code{registerDoMC(3)} to register parallel.
#'
#' @return a list of lasso parameter
#'
#' @importFrom doMC registerDoMC
#' @keywords internal
#'
#' @noRd
#'
runlassolist <- function(
  x,
  parallel = T,
  ...){
  y <- x
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  CDFBETA <- generateBetaCDF(
    # alpha = a1, beta = a2,
    index.p = grid.p, ...)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
  # set.seed(12345)                                                     # || double check
  lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
  if (parallel){
    # message("parallel =", parallel)
    registerDoMC(3)
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = T)
  } else {
    # message("parallel =", parallel)
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = F)
  }
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[zeros == FALSE]
  # gc(verbose = FALSE)
  return(selects)
}

## 6.4 catNorm ----
#' Add 1 to all individual basis set
#'
#' Add Gaussian basis to all individual basis set
#'
#' @param vector vector:basis set
#'
#' @return vector:basis set including 1
#'
#' @keywords internal
#'
#' @noRd
#'
catNorm <- function(
  vector
  ) {
  unique(c(0, 1, sort(vector, method = "quick")))
}

## 6.5 replist ----
#' Generate one vectors with the length of each individual basis set
#'
#' @param vector basis set
#'
#' @return vector one vectors with the length of each individual basis set
#'
#' @keywords internal
#'
#' @noRd
#'
replist <- function(
  vector
  ) {
  rep(1, length(vector))
}

## 6.6 countBasis ----
#' Generate nested subsets of basis
#'
#' Number basis in each nested subsets, and  each count(# of votes) yielding each nested basis
#'
#' @param nonzero.list sparse set Di
#' @param nonzero.obs one vectors with the same length to Di
#'
#' @return list \code{[[1]]}:the nested subsets of basis, \code{[[2]]}:# basis in each nested subset \code{[[3]]}: each count(# of votes) yielding each nested basis.
#'
#' @keywords internal
#'
#' @noRd
#'
countBasis <- function(
  nonzero.list,
  nonzero.obs,
  ...
){
  colum <- unlist(nonzero.list) ## length(colum )
  obs <- unlist(nonzero.obs) ## length(obs)
  mytable_colum0 <- table(colum, obs)
  unique.colum <- as.numeric(rownames(margin.table(mytable_colum0, 1))) ## length( unique.colum )
  unique.count <- as.numeric(margin.table(mytable_colum0, 1)) ## length( unique.count )
  margin.count <- c(0, sort(unique(unique.count), method = "quick")) ## length(margin.count)
  long.set <- length(margin.count) - 1
  list.set <- vector("list", long.set)
  for (i in 1:long.set) {
    list.set[[i]] <- sort(unique.colum[unique.count > margin.count[i]], method = "quick")
  }
  outputs <- list(list.set, unlist(lapply(list.set, length)), margin.count[-length(margin.count)])
  return(outputs)
}

## 6.7 incidenceVec ----
#' Compute frequency for each selected basis
#'
#' @param setlistofobs basis set(CatNorm)
#' @param frqlistofobs frequency for each count(rep.list)
#'
#' @return list \code{[[1]]}:basis names across subjects, \code{[[2]]}:frequency for each basis across subjects.
#'
#' @keywords internal
#'
#' @noRd
incidenceVec <- function(
  setlistofobs,
  frqlistofobs
  ) {
  colum <- unlist(setlistofobs) ## length(colum )
  obs <- unlist(frqlistofobs) ## length( obs)
  mytable_colum0 <- table(colum, obs)
  unique.colum <- as.numeric(rownames(margin.table(mytable_colum0, 1)))
  unique.count <- as.numeric(margin.table(mytable_colum0, 1))
  outputs <- list(unique.colum, unique.count)
  return(outputs)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 7.1 lassolist2 ----
#' Run lasso list to build basis function
#'
#'
#' @param object A Qboneobject
#' @param verbose Print a progress bar
#' @param new.assay.name new assay name assigned to the lassolist data
#' @param alpha vector containing sequence of beta parameter for internal function \code{generateBetaCDF()}
#' @param beta  vector containing sequence of beta parameter for internal function \code{generateBetaCDF()}
#' @param assay.seed assay information to add into the QboneData object scale.data. The default of \code{lassolist()} will save the random seed for the run. \code{assay.seed = object@assays[["Lasso"]]@scale.data[["lassolist"]]} to run \code{lassolist()} for the same results.
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must register parallel before hand
#' @param ... Arguments passed to other methods
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach %dopar%
#'
#'
#' @export
#'
lassolist2 <- function(
  object,
  verbose = F,
  new.assay.name = "Lasso",
  alpha = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  beta = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  assay.seed = .Random.seed,
  parallel = T,
  ...
){
  .Random.seed <- assay.seed
  orig.dataset <- getQboneData(object, slot = 'data', assay = defaultAssay(object))
  # iterate.fxn <- ifelse(test = verbose, yes = pblapply, no = lapply)
  # new.data <- iterate.fxn(
  #   X = orig.dataset,
  #   FUN = runlassolist,
  #   alpha = a1,
  #   beta = a2
  # )
  if (parallel & verbose == F) {
    # .Random.seed <- assay.seed
    new.data <- foreach(i = 1:length(orig.dataset)) %dopar% { #, .packages=c("glmnet", "doMC")
      runlassolist2(orig.dataset[[i]],
                            alpha = alpha,
                            beta = beta,
                            parallel = parallel,
                            assay.seed = assay.seed, ...)
    }

  } else if ( parallel & verbose) {
    # .Random.seed <- assay.seed
    pb <- txtProgressBar(min = 0, max = length(orig.dataset), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    new.data <- foreach(i = 1:length(orig.dataset), .options.snow=opts) %dopar% { # , .packages="Qbone", .packages=c("glmnet", "doMC")
      runlassolist2(orig.dataset[[i]],
                    alpha = alpha,
                    beta = beta,
                    parallel = parallel,
                    assay.seed = assay.seed, ...)
    }
    close(pb)
  } else {stop()}
  # if (verbose) {
  #   # .Random.seed <- assay.seed
  #   pb <- txtProgressBar(min = 0, max = length(orig.dataset), style = 3)
  #   new.data <- list()
  #   for (i in 1:length(orig.dataset)){
  #     bar1 <- runlassolist2(orig.dataset[[i]],
  #                           alpha = alpha,
  #                           beta = beta,
  #                           parallel = parallel, ...)
  #     new.data <- append(new.data, list(bar1))
  #     setTxtProgressBar(pb, i)
  #   }
  #   close(pb)
  #   # stopCluster(cl)
  # } else {
  #   # .Random.seed <- assay.seed
  #   new.data <- lapply(
  #     X = orig.dataset,
  #     FUN = runlassolist2,
  #     alpha = alpha,
  #     beta = beta,
  #     parallel = parallel
  #   )
  # }
  new.qbonedata <- createQboneData(new.data,
                                   meta.assays = data.frame(id = names(orig.dataset)),
                                   sampleid.assays = 1,
                                   assay.name = new.assay.name,
                                   assay.orig = defaultAssay(object))
  new.qbonedata@scale.data <- append(object@assays[[defaultAssay(object)]]@scale.data,list(lassolist = c(assay.seed)))
  object[[new.assay.name]] <- new.qbonedata
  defaultAssay(object) <- new.assay.name
  return(object)
}

### 7.1.2 runlasso2 ----
#' Run lasso for list of data
#'
#' @param x any list of data
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must register parallel before hand
#'
#' @return a list of lasso parameter
#'
#' @importFrom doMC registerDoMC
#'
#' @keywords internal
#'
#' @noRd
#'
runlassolist2 <- function(x, parallel = parallel, assay.seed = assay.seed, ...){
  .Random.seed <- assay.seed
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
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[zeros == FALSE]
  # selects2 <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[zeros2 == FALSE]
  # gc(verbose = FALSE)
  return(selects)
}

### 7.1.3 generateBetaCDF ----
#' creates beta basis functions
#'
#'
#' @param alpha vector containing sequence of beta parameter
#' @param beta  vector containing sequence of beta parameter
#' @param index.p probability grids on (0,1)
#'
#' @return matrix containing # of beta parameters times the length of index.p
#'
#' @keywords internal
#'
#' @noRd
#'
GENERATE_BETA_CDF <- function(alpha, beta, index.p, ...) {
  n1 <- as.character(alpha)
  n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(n2) * length(n1))
  for (i in 1:length(n1)) { ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(n2)) {
      rowth <- (j - 1) * length(n2) + i
      BETASCDF0[rowth, ] <- pbeta(index.p, alpha[i], beta[j])
      BETASCDF0[rowth, ] <- round(centering.function(BETASCDF0[rowth, ], scale = TRUE), 7)
    }
  }
  ## name.mat = outer( paste0( "(", n1 , se="") ,paste0(",", n2, ")" , se="")  ,  FUN=paste ,sep="")
  ## matrix.rownames = as.vector(  name.mat )
  ## outputs= list( BETASCDF0 ,  name.mat , matrix.rownames )
  outputs <- BETASCDF0
  return(outputs)
}

### 7.1.4 centering.function ----
#' computes centering and scaling column by column
#'
#' @param raw.x any matrix with n times p
#' @param scale FALSE or TRUE (doing only centering if scale=FALSE )
#'
#' @return the matrix centerted or/and scaled
#'
#' @keywords internal
#'
#' @noRd
#'
centering.function <- function(raw.x, scale = FALSE) {
  object.x <- as.matrix(raw.x)
  x.n <- dim(object.x)[1]
  one <- rep(1, x.n)
  meanx <- drop(one %*% object.x) / x.n
  mean1 <- as.matrix(meanx)
  n1 <- dim(mean1)[1]
  jj <- rep(1, x.n) %*% t(rep(1, dim(object.x)[2]))
  mean.mat <- diag(meanx, nrow = n1, ncol = n1)
  cen.x <- object.x - jj %*% mean.mat
  if (scale == FALSE) {
    return(cen.x)
  }
  if (scale == TRUE) {
    normx <- sqrt(drop(one %*% (cen.x^2)))
    normx1 <- as.matrix(normx)
    s1 <- dim(normx1)[1]
    scale.mat <- diag(1 / normx, nrow = s1, ncol = s1)
    cen.scale.x <- cen.x %*% scale.mat
    return(cen.scale.x)
  }
}
