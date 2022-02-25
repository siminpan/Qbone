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
#' @param assay.seed assay information to add into the QboneData object scale.data. The default of \code{lassolist()} will save the random seed for the run. \code{.Random.seed <- boject@assays[["Lasso"]]@scale.data[["lassolist"]] } to reset the seed for the same results.
#' @param ... Arguments passed to other methods
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom pbapply pblapply
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
  ...
){
  orig.dataset <- getQboneData(object, slot = 'data', assay = defaultAssay(object))
  iterate.fxn <- ifelse(test = verbose, yes = pblapply, no = lapply)
  new.data <- iterate.fxn(
    X = orig.dataset,
    FUN = runlassolist,
    alpha = a1,
    beta = a2
  )
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
#'
#' @return a list of lasso parameter
#'
#' @keywords internal
#'
#' @noRd
#'
runlassolist <- function(x, ...){
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
  cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3)
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
  # gc(verbose = FALSE)
  return(selects)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
