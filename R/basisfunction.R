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

# 2.1 lassolist ----
#' Run lasso list to build basis function
#'
#'
#'
#' @param object A Qboneobject
#' @param verbose Print a progress bar
#' @param new.assay.name new assay name assigned to the lassolist data
#' @param a1 vector containing sequence of beta parameter for internal function \code{generateBetaCDF()}
#' @param a2  vector containing sequence of beta parameter for internal function \code{generateBetaCDF()}
#' @param assay.seed assay information to add into the QboneData object scale.data. The default of \code{lassolist()} will save the random seed for the run. \code{assay.seed = object@assays[["Lasso"]]@scale.data[["lassolist"]]} to run \code{lassolist()} for the same results.
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must register parallel before hand
#' @param ... Arguments passed to other methods
#'
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @export
#'
lassolist <- function(
  object,
  verbose = TRUE,
  new.assay.name = "Lasso",
  a1 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  a2 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
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
  if (verbose) {
    # cl <- makeCluster(detectCores()-2)
    new.data <- lapply(
      X = orig.dataset,
      FUN = runlassolist,
      alpha = a1,
      beta = a2,
      parallel = parallel
    )
    # stopCluster(cl)
  } else {
    new.data <- lapply(
      X = orig.dataset,
      FUN = runlassolist,
      alpha = a1,
      beta = a2,
      parallel = parallel
    )
  }
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
generateBetaCDF <- function(alpha, beta, index.p, ...) {
  # n1 <- as.character(alpha)
  # n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(beta) * length(alpha))
  for (i in 1:length(alpha)) { ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(beta)) {
      rowth <- (j - 1) * length(beta) + i
      BETASCDF0[ rowth, ] <- pbeta(index.p, alpha[i], beta[j])
      BETASCDF0[ rowth, ] <- round(centeringFunction(BETASCDF0[ rowth, ], scale = TRUE), 7)
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
centeringFunction <- function(raw.x, scale = FALSE) {
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
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must register parallel before hand
#'
#' @return a list of lasso parameter
#'
#' @importFrom doMC registerDoMC
#' @keywords internal
#'
#' @noRd
#'
runlassolist <- function(x, parallel = T, ...){
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
    registerDoMC(3)
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = T)
  } else {
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = F)
  }
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
  # gc(verbose = FALSE)
  return(selects)
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
  # iterate.fxn <- ifelse(test = verbose, yes = parLapply, no = lapply)
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(orig.dataset), style = 3)
    new.data <- list()
    for (i in 1:length(orig.dataset)){
      bar1 <- runlassolist2(orig.dataset[[i]],
                            alpha = alpha,
                            beta = beta,
                            parallel = parallel, ...)
      new.data <- append(new.data, list(bar1))
      setTxtProgressBar(pb, i)
    }
    close(pb)
    # stopCluster(cl)
  } else {
    new.data <- lapply(
      X = orig.dataset,
      FUN = runlassolist2,
      alpha = alpha,
      beta = beta,
      parallel = parallel
    )
  }
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

## 7.2 runlasso2 ----
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
  lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
  # set.seed(12345)
  # cvfit.lasso <- cv.biglasso(as.big.matrix(BETA_BASE_TOTAL_2), y, nfolds = 3) #, intercept = TRUE, nfolds = 3)
  # summary(cvfit.lasso)
  if (parallel){
    registerDoMC(3)
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = T)
  } else {
    cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3, parallel = F)
  }
  # set.seed(12345)
  # cvfit.lasso2 <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 5, parallel = T)
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  # zeros2 <- as.vector(coef(lasso_fit2, s = cvfit.lasso2$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
  # selects2 <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros2 == FALSE ]
  # gc(verbose = FALSE)
  return(selects)
}

## 7.3 generateBetaCDF ----
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
      BETASCDF0[ rowth, ] <- pbeta(index.p, alpha[i], beta[j])
      BETASCDF0[ rowth, ] <- round(centering.function(BETASCDF0[ rowth, ], scale = TRUE), 7)
    }
  }
  ## name.mat = outer( paste0( "(", n1 , se="") ,paste0(",", n2, ")" , se="")  ,  FUN=paste ,sep="")
  ## matrix.rownames = as.vector(  name.mat )
  ## outputs= list( BETASCDF0 ,  name.mat , matrix.rownames )
  outputs <- BETASCDF0
  return(outputs)
}

## 7.4 centering.function ----
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
  normx <- sqrt(drop(one %*% (cen.x^2)))
  normx1 <- as.matrix(normx)
  s1 <- dim(normx1)[1]
  scale.mat <- diag(1 / normx, nrow = s1, ncol = s1)
  cen.scale.x <- cen.x %*% scale.mat
  if (scale == FALSE) {
    return(cen.x)
  }
  if (scale == TRUE) {
    return(cen.scale.x)
  }
}
