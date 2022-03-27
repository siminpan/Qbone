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

# Function logic based on R Object ----

## 2.1 lassoList ----
#' Use penalized regression (lasso) to find a sparse subset of dictionary
#' elements
#'
#' First construct overcomplete dictionary (Beta CDF). Then uses penalized
#' regression (lasso) to find a sparse subset of dictionary elements.
#'
#' @param object A Qboneobject
#' @param verbose Print a progress bar
#' @param new.assay.name New assay name assigned to the lassolist data
#' @param data.assay It is the name of the assay whose data will be used to
#' compute the lasso list. Default is the data from the defaultAssay(object).
#' @param alpha Vector containing sequence of beta parameter for internal
#' function \code{generateBetaCDF()}
#' @param beta  Vector containing sequence of beta parameter for internal
#' function \code{generateBetaCDF()}
#' @param assay.seed assay information to add into the QboneData object
#' scale.data. The default of \code{lassolist()} will save the random seed for
#'  the run.
#'  Use \code{.Random.seed <-  object@assays[["Lasso.list"]]@scale.data[["lassolist"]] }
#'  before run \code{lassolist()} for the same results.
#' @param parallel If TRUE, use parallel foreach to fit each fold in
#' \code{cv.glmnet()}. Default use \code{registerDoMC()} to register parallel.
#' There is another function \code{lassolist_parallel()} for overall parallel
#' computing for all sample.
#' @param ... Arguments passed to other methods
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
lassoList <- function(
  object,
  verbose = TRUE,
  new.assay.name = "Lasso.list",
  data.assay = defaultAssay(object),
  alpha = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  beta = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  assay.seed = .Random.seed,
  parallel = T,
  ...
){
  # Check data.assay
  if (data.assay == "Lasso.list"){
    stop('The default assay for this Qbone object is already "Lasso.list".')
  }
  # Assign random seed
  .Random.seed <- assay.seed
  # Get data
  orig.dataset <- getQboneData(object, slot = 'data', assay = data.assay)

  # build overcomplete dictionary
  # (The overcomplete dictionary is too big not practical to store)
  # message("\n Building overcomplete dictionary (Beta CDF)")
  # betaDictionary <- lapply(
  #   X = orig.dataset,
  #   FUN = buildDictionary,
  #   alpha = alpha,
  #   beta = beta
  # )
  # dictionary.qbonedata <- createQboneData(betaDictionary,
  #                                  meta.assays = data.frame(id = names(orig.dataset)),
  #                                  sampleid.assays = 1,
  #                                  assay.name = "oDictionary",
  #                                  assay.orig = data.assay)
  # object[["oDictionary"]] <- dictionary.qbonedata
  # betaDictionary <- getQboneData(object, slot = 'data', assay = "oDictionary")

  # Penalized regression (lasso)
  message("\n This step may take a while as it involves k-fold cross-validation.")
  if (verbose){
    message("\n Using penalized regression (lasso) to find a sparse subset of dictionary elements.")
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
  # Build new Qbone object
  new.qbonedata <- createQboneData(new.data,
                                   meta.assays = data.frame(id = names(orig.dataset)),
                                   sampleid.assays = 1,
                                   assay.name = new.assay.name,
                                   assay.orig = data.assay)
  new.qbonedata@scale.data <- append(object@assays[[data.assay]]@scale.data,
                                     list(lassolist = c(assay.seed)))
  object[[new.assay.name]] <- new.qbonedata
  defaultAssay(object) <- new.assay.name
  return(object)
}

## 2.2 preQuantlets ----
#' Get the pre-quantlets basis functions
#'
#' Take union set and find near-lossless subset as confirmed by
#' cross-validated concordance correlation coefficient (CVCCC)
#'
#' @param object A Qboneobject
#' @param new.assay.name New assay name assigned to the quantlets data
#' @param data.assay It is the name of the assay whose data will be used
#'  to compute the lasso list. Default is the data from the
#'  defaultAssay(object).
#' @param p Vector of length P  in (0,1)	Probability grids.
#' @param alpha Vector containing sequence of beta parameter for internal
#'  function \code{generateBetaCDF()}
#' @param beta  Vector containing sequence of beta parameter for internal
#'  function \code{generateBetaCDF()}
#' @param ... Arguments passed to other methods
#'
#' @export
#'
preQuantlets <- function(
  object,
  new.assay.name = "Pre.Quantiles",
  data.assay = defaultAssay(object),
  p = signif(seq(0.001, 0.999, length = 1024), 4),
  alpha = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  beta = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1)),
  ...
){
  # Check data.assay
  if(data.assay != "Lasso.list"){
    warning('The default assay is not "Lasso.list" please double the defaultAssay() of this Qbone object. This step should be run on results of lassoList().')
  }
  message('Will compute the quantlets basis functions based on "', data.assay, '" results from "', object@assays[["Lasso.list"]]@assay.orig, '" data.', "\n This step may take a while.")
  # Get data
  raw.dataset <- getQboneData(object, slot = 'data', assay = object@assays[[data.assay]]@assay.orig)
  lasso.dataset <- getQboneData(object, slot = 'data', assay = data.assay)

  # Generate nested subsets of basis
  lasso.list1 <- lapply(lasso.dataset, catNorm) ## 1st term fix
  lasso.nonzero.obs1 <- lapply(lasso.list1, replist) ## D_i, generate 1 vector
  lasso.counts.fit <- countBasis(lasso.list1, lasso.nonzero.obs1)

  # Compute frequency for each selected basis
  n <- length(raw.dataset)
  lasso_IncidenceVec_i_ <- vector("list", n)
  for (i in 1:n){
    lasso_fitIncd <- incidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
    lasso_IncidenceVec_i_[[i]] <- lasso_fitIncd
  }

  # Leave-one-out concordance correlation
  lasso.locc <- locc(
    leaveout.list = lasso_IncidenceVec_i_,
    remain.counts = lasso.counts.fit[[3]],
    remain.basis = lasso.counts.fit[[2]],
    Y.list = raw.dataset,
    maxim = length(p),
    alpha = alpha,
    beta = beta
  )

  # Pre-quantlets basis functions, unify grid set as p1024 to avoid problems like gramshumit, memory error, irregar, empiricalQ.
  betaCDF <- generateBetaCDF(alpha = alpha,
                             beta = beta,
                             index.p = p)
  NQ <- qnorm(p)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(betaCDF))
  reduced_BASE <- BETA_BASE_TOTAL_2

  q.raw.dataset <- lapply(raw.dataset,
                          function(x, p){quantile(x, p, type = 7)},
                          p = p
                          ) # Empirical quantiles for each subject.

  # Build new Qbone object
  new.qbonedata <- createQboneData(q.raw.dataset,
                                   meta.assays = data.frame(id = names(raw.dataset)),
                                   sampleid.assays = 1,
                                   assay.name = new.assay.name,
                                   assay.orig = data.assay)
  new.qbonedata@scale.data <- append(object@assays[[data.assay]]@scale.data,
                                     list(betaCDF = reduced_BASE, # Pre-quantlets basis functions, which implies the selected beta cdf functions. (P by K matrix. P: desired number of probability grids, and K: number of basis).
                                          locc = lasso.locc, # LOCC value for each choice of the basis function (Z: the length of the possible choices for basis).
                                          basis.columns = lasso.counts.fit[[1]], # Corresponding basis columns for each possible choice.
                                          remain.basis = lasso.counts.fit[[2]], # Number of the basis for each possible choice.
                                          remain.counts = lasso.counts.fit[[3]] # Frequency for each possible choice.
                                          )
                                     )
  object[[new.assay.name]] <- new.qbonedata
  defaultAssay(object) <- new.assay.name
  return(object)
}

## 2.3 ecQuantlets ----
#' Get the Compute Empirical Coefficients for quantlets
#'
#' Compute empirical coefficients for quantlets
#'
#' @param object A Qboneobject
#' @param new.assay.name New assay name assigned to the quantlets data
#' @param data.assay It is the name of the assay whose data will be used
#'  to compute the lasso list. Default is the data from the
#'  defaultAssay(object).
#' @param k number of basis coefficients to keep. This is based on the number
#'  of dictionary elements (C) to keep. Default will pick the one chosen from
#'   \code{dxPlot()}.
#' @param sparsity Sparsity regularization parameter.
#' @param ... Arguments passed to other methods
#'
#' @importFrom pracma gramSchmidt
#'
#' @export
#'

ecQuantlets <- function(
  object,
  new.assay.name = "Empirical.Coefficients",
  data.assay = defaultAssay(object),
  k = NULL,
  sparsity = 0.001,
  ...
){
  # Check data.assay
  if(data.assay != "Pre.Quantiles"){
    warning('The default assay is not "Pre.Quantiles" please double the defaultAssay() of this Qbone object. This step should be run on results of preQuantlets().')
  }
  # Get data and compute the reduce basis
  orig.dataset <- raw.dataset <- getQboneData(object, slot = 'data',
                                              assay = object@assays[[
                                                "Lasso.list"
                                                #object@assays[[defaultAssay(object)]]@assay.orig
                                              ]]@assay.orig
  )
  if (is.null(k)){
    # plotdata1 <- loccplotdata(object, ...)
    plotdata2 <- loccplotdata(object, ...)[[3]]
    lasso.mean_i_ = plotdata2$y[plotdata2$group == "mean"]
    lasso.min_i_ = plotdata2$y[plotdata2$group == "min"]
    basis.columns.no = min(plotdata2$x[c(lasso.mean_i_ - lasso.min_i_) > sparsity])
    k = max(plotdata2$x2[c(lasso.mean_i_ - lasso.min_i_) > sparsity])
    basis.columns.select = object@assays[["Pre.Quantiles"]]@scale.data[["basis.columns"]][[basis.columns.no]]
    reduceBasis = object@assays[["Pre.Quantiles"]]@scale.data[["betaCDF"]][, basis.columns.select]
  } else {
    if (k %in% unlist(lapply(object@assays[["Pre.Quantiles"]]@scale.data[["basis.columns"]], length), use.names = F)){
      basis.columns.no = which(k ==unlist(lapply(object@assays[["Pre.Quantiles"]]@scale.data[["basis.columns"]], length), use.names = F))
      basis.columns.select = object@assays[["Pre.Quantiles"]]@scale.data[["basis.columns"]][[basis.columns.no]]
      reduceBasis = object@assays[["Pre.Quantiles"]]@scale.data[["betaCDF"]][, basis.columns.select]
    } else {
      stop("Based on overcomplete dictionary, k must be one of the following numbers ",
           paste0(unlist(lapply(object@assays[["Pre.Quantiles"]]@scale.data[["basis.columns"]], length), use.names = F), sep = ", "))
    }
  }
  # Update Pre.Quantiles with Reduce.basis in Qbone object
  object@assays[[data.assay]]@scale.data <- append(object@assays[[data.assay]]@scale.data,
                                     list(C = basis.columns.no,
                                          k = k,
                                          basis.columns.select = basis.columns.select,
                                          reduceBasis = reduceBasis
                                          )
                                     )
  names(object@assays) <- ifelse(names(object@assays) == data.assay, "Reduce.basis", names(object@assays))
  object@assays[["Reduce.basis"]]@assay.name <- "Reduce.basis"
  object@assays[["Reduce.basis"]]@assay.orig <- data.assay
  defaultAssay(object) <- "Reduce.basis"
  # Orthogonalization
  basis.columns.select.no = object@assays[["Reduce.basis"]]@scale.data[["basis.columns"]][[basis.columns.no + 10]]
  reduceBasis.norms <- object@assays[["Reduce.basis"]]@scale.data[["betaCDF"]][, basis.columns.select.no]
  gramS <- gramSchmidt(reduceBasis, tol = .Machine$double.eps^0.5)
  norms <- gramSchmidt(as.matrix(reduceBasis.norms), tol = .Machine$double.eps^0.5)$Q
  # Denoising
  quantlets <- waveletDe(gramS$Q[, -1], filter.number = 2, family = "DaubExPhase")[, , 2]
  quantlets.ns <- cbind(norms, centering.function(quantlets, scale = TRUE))
  # Compute Empirical Coefficients
  # eQ = getQboneData(object, slot = 'data', assay = defaultAssay(object))
  # eQy = matrix(round(unlist(eQ, use.names = F),3), 1024)
  eQy = matrix(
    round(
      unlist(getQboneData(object, slot = 'data', assay = defaultAssay(object)),
             use.names = F
             ),
      3
      ),
    dim(object@assays[["Reduce.basis"]]@scale.data[["reduceBasis"]])[1]
    )
  emp.raw.dataset <- vector("list", length(orig.dataset))
  for (i in 1:length(orig.dataset)) {
    emp.raw.dataset[[i]] <- eQy[, i]
  }
  empCoefs <- empCoefs(quantlets.ns, emp.raw.dataset)
  # Build new Qbone object
  empCoefs.list <- split(empCoefs, seq(nrow(empCoefs)))
  names(empCoefs.list) <- object@meta.data[["id"]]
  # get it back
  # empCoefs2 <- matrix(unlist(empCoefs.list), nrow=length(empCoefs.list), byrow=TRUE)
  new.qbonedata <- createQboneData(empCoefs.list,
                                   meta.assays = data.frame(id = names(orig.dataset)),
                                   sampleid.assays = 1,
                                   assay.name = new.assay.name,
                                   assay.orig = data.assay)
  new.qbonedata@scale.data <- append(object@assays[[defaultAssay(object)]]@scale.data,
                                     list(quantlets = quantlets.ns))
  object[[new.assay.name]] <- new.qbonedata
  defaultAssay(object) <- new.assay.name
  return(object)
}

# Function logic based on Dic ----
## 2.4 ocDic ----

## 2.5 unionDic ----

## 2.6 cmnDic ----

## 2.7 orthogDic ----

## 2.8 denoiseDic ----

## 2.9 quantletsDic ----

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
#' @param ... Arguments passed to other methods
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
  ){
  # n1 <- as.character(alpha)
  # n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(beta) * length(alpha))
  for (i in 1:length(alpha)){ ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(beta)){
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
  ){
  object.x <- as.matrix(raw.x)
  x.n <- dim(object.x)[1]
  one <- rep(1, x.n) # matrix(rep(1, x.n), ncol = x.n)
  meanx <- drop(eigenmapmm(matrix(rep(1, x.n), ncol = x.n), object.x)) / x.n
  mean1 <- as.matrix(meanx)
  n1 <- dim(mean1)[1]
  jj <- eigenmapmm(one,t(rep(1, dim(object.x)[2]))) # rep(1, x.n) %*% t(rep(1, dim(object.x)[2])) # one %*% two # eigenmapmm(one, two)
  mean.mat <- diag(meanx, nrow = n1, ncol = n1)
  cen.x <- object.x - eigenmapmm(jj, mean.mat)
  if (scale == FALSE){
    return(cen.x)
  }
  if (scale == TRUE){
    normx <- sqrt(drop(eigenmapmm(matrix(rep(1, x.n), ncol = x.n), (cen.x^2))))
    normx1 <- as.matrix(normx)
    s1 <- dim(normx1)[1]
    scale.mat <- diag(1 / normx, nrow = s1, ncol = s1)
    cen.scale.x <- eigenmapmm(cen.x, scale.mat)
    return(cen.scale.x)
  }
}

## 6.3 buildDictionary ----
#' Build dictionary
#'
#' (The overcomplete dictionary is too big not practical to store)
#' Construct an dictionary that contains bases spanning the space
#' of Gaussian quantile functions plus a large number of Beta cumulative
#'density functions.
#'
#' @param x any list of data
#' @param ... Arguments passed to other methods
#'
#' @return matrix containing dictionary (Beta CDF)
#'
#' @keywords internal
#'
#' @noRd
#'
buildDictionary <- function(
  x,
  ...
){
  y <- x
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  betaCDF <- generateBetaCDF(
    # alpha = a1, beta = a2,
    index.p = grid.p, ...)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(betaCDF))
  return(BETA_BASE_TOTAL_2)
}

## 6.4 runlasso ----
#' Penalized regression (lasso)
#'
#' First construct overcomplete dictionary (Beta CDF). Then uses penalized regression (lasso) to find a sparse subset of dictionary elements.
#'
#' @param x any list of data
#' @param parallel If TRUE, use parallel foreach to fit each fold in \code{cv.glmnet(..., nfolds = 3)}. Default use \code{registerDoMC(3)} to register parallel.
#' @param ... Arguments passed to other methods
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
  betaCDF <- generateBetaCDF(
    # alpha = a1, beta = a2,
    index.p = grid.p, ...)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(betaCDF))
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

## 6.5 catNorm ----
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
  ){
  unique(c(0, 1, sort(vector, method = "quick")))
}

## 6.6 replist ----
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
  ){
  rep(1, length(vector))
}

## 6.7 countBasis ----
#' Generate nested subsets of basis
#'
#' Number basis in each nested subsets, and  each count(number of votes) yielding each nested basis
#'
#' @param nonzero.list sparse set Di
#' @param nonzero.obs one vectors with the same length to Di
#' @param ... Arguments passed to other methods
#'
#' @return list \code{[[1]]}: the nested subsets of basis,
#'              \code{[[2]]}: number of basis in each nested subset
#'              \code{[[3]]}: each count(number of votes) yielding each nested basis.
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
  for (i in 1:long.set){
    list.set[[i]] <- sort(unique.colum[unique.count > margin.count[i]], method = "quick")
  }
  outputs <- list(list.set, unlist(lapply(list.set, length)), margin.count[-length(margin.count)])
  return(outputs)
}

## 6.8 incidenceVec ----
#' Compute frequency for each selected basis
#'
#' @param setlistofobs basis set(CatNorm)
#' @param frqlistofobs frequency for each count(rep.list)
#'
#' @return list \code{[[1]]}:basis names across subjects,
#'              \code{[[2]]}:frequency for each basis across subjects.
#'
#' @keywords internal
#'
#' @noRd
incidenceVec <- function(
  setlistofobs,
  frqlistofobs
  ){
  colum <- unlist(setlistofobs) ## length(colum )
  obs <- unlist(frqlistofobs) ## length( obs)
  mytable_colum0 <- table(colum, obs)
  unique.colum <- as.numeric(rownames(margin.table(mytable_colum0, 1)))
  unique.count <- as.numeric(margin.table(mytable_colum0, 1))
  outputs <- list(unique.colum, unique.count)
  return(outputs)
}

## 6.9 locc ----
#' Computes the leave-one-out concordance correlation index
#'
#' @param leaveout.list leave-one-out list (it may be output of incidenceVec function)
#' @param remain.counts vector containing c value in paper  (it may be output[[2]] of countBasis)
#' @param remain.basis vector containing k(c) value in paper  (it may be output[[3]] of countBasis)
#' @param Y.list raw.dataset
#' @param maxim scalar value (length of probability grids)
#' @param ... Arguments passed to other methods
#'
#' @return Predicted values ( length(probability grids) times n times length(c)  array output ) based on leaveout.list
#'
#' @importFrom MASS ginv
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @keywords internal
#'
#' @noRd
#'
locc <- function(
  leaveout.list,
  remain.counts,
  remain.basis,
  Y.list,
  maxim = 1000,
  ...
  ){
  n <- length(leaveout.list)
  active.set <- (remain.basis < maxim)
  if (tail(remain.counts[active.set], 1) == n - 1){
    remain.counts[ length(remain.counts)] <- n - 2
  }
  feasible.long <- sum(active.set)
  max.long <- max(unlist(lapply(Y.list, length), use.names = FALSE))
  # checks <- matrix(NA, nrow = n, ncol = feasible.long)
  Values <- array(NA, c(max.long, n, feasible.long))

  pb <- txtProgressBar(min = 0, max = length(Y.list), style = 3)
  for (i in 1:n){
    # message("\n Computes the leave-one-out concordance correlation index for ", names(Y.list)[i])
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
    for (j in 1:feasible.long){ ###  j = 20

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
  return(outputs)
}

## 6.10 waveletDenose ----
#' Conducts Wavelet Denoising
#'
#' @param object n by p matrix
#' @param ... Arguments passed to other methods
#'
#' @return array[,,1] includes denoised values based on the usual wavelet
#' method. array[,,2] includes denoised values based on nondecimal method
#'
#' @importFrom wavethresh wd accessD nlevelsWT threshold wst wr AvBasis
#' @importFrom stats mad
#'
#' @keywords internal
#'
#' @noRd
#'
waveletDe <- function(
  object,
  filter.number = 2,
  family = "DaubExPhase",
  ...
){
  object <- as.matrix(object)
  n <- dim(object)[1]
  p <- dim(object)[2]
  outputs <- array(NA, c(n, p, 2))
  for (i in 1:p) {
    ywd <- wd(object[, i], filter.number = filter.number, family = family)
    fineCoefs <- accessD(ywd, lev = nlevelsWT(ywd) - 1)
    sigma <- mad(fineCoefs)
    utDJ <- sigma * sqrt(2 * log(1024))
    ywdT <- threshold(ywd, policy = "manual", value = utDJ)
    outputs[, i, 1] <- wr(ywdT)

    ywst <- wst(object[, i], filter.number = filter.number, family = family)
    fineWSTCoefs <- accessD(ywst, lev = nlevelsWT(ywst) - 1)
    sigmaWST <- mad(fineWSTCoefs)
    utWSTDJ <- sigmaWST * sqrt(2 * log(1024))
    ywstT <- threshold(ywst, policy = "manual", value = utWSTDJ)
    outputs[, i, 2] <- AvBasis(ywstT)
  }
  return(outputs)
}

## 6.11 empCoefs----
#' Compute Empirical Coefficients
#'
#' @param quantlet quanvelts basis
#' @param raw.dataset input dataset
#' @param ... Arguments passed to other methods
#'
#' @return mats
#'
#' @keywords internal
#'
#' @noRd
#'

empCoefs <- function(
  quantlet,
  raw.dataset,
  ...
  ){
  n <- length(raw.dataset)
  mats <- matrix(NA, nrow = n, ncol = (dim(quantlet)[2] + 1))
  for (i in 1:n) {
    y <- raw.dataset[[i]]
    y.long <- length(y)
    set.seed(123 + i)
    grids <- sort(sample(dim(quantlet)[1], y.long))
    Psi <- cbind(rep(1, length(grids)), quantlet[grids, ])
    mats[i, ] <- solve(t(Psi) %*% Psi) %*% t(Psi) %*% y
  }
  return(mats)
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
#' @param assay.seed assay information to add into the QboneData object scale.data. The default of \code{lassolist()} will save the random seed for the run. \code{assay.seed = object@assays[["Lasso.list"]]@scale.data[["lassolist"]]} to run \code{lassolist()} for the same results.
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
  new.assay.name = "Lasso.list",
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
  if (parallel & verbose == F){
    # .Random.seed <- assay.seed
    new.data <- foreach(i = 1:length(orig.dataset)) %dopar% { #, .packages=c("glmnet", "doMC")
      runlassolist2(orig.dataset[[i]],
                            alpha = alpha,
                            beta = beta,
                            parallel = parallel,
                            assay.seed = assay.seed, ...)
    }

  } else if ( parallel & verbose){
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
  # if (verbose){
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
GENERATE_BETA_CDF <- function(alpha, beta, index.p, ...){
  n1 <- as.character(alpha)
  n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(n2) * length(n1))
  for (i in 1:length(n1)){ ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(n2)){
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
centering.function <- function(raw.x, scale = FALSE){
  object.x <- as.matrix(raw.x)
  x.n <- dim(object.x)[1]
  one <- rep(1, x.n)
  meanx <- drop(one %*% object.x) / x.n
  mean1 <- as.matrix(meanx)
  n1 <- dim(mean1)[1]
  jj <- rep(1, x.n) %*% t(rep(1, dim(object.x)[2]))
  mean.mat <- diag(meanx, nrow = n1, ncol = n1)
  cen.x <- object.x - jj %*% mean.mat
  if (scale == FALSE){
    return(cen.x)
  }
  if (scale == TRUE){
    normx <- sqrt(drop(one %*% (cen.x^2)))
    normx1 <- as.matrix(normx)
    s1 <- dim(normx1)[1]
    scale.mat <- diag(1 / normx, nrow = s1, ncol = s1)
    cen.scale.x <- cen.x %*% scale.mat
    return(cen.scale.x)
  }
}

### 7.2.1 LCCC.FIT_1 ----
#' computes the leave-one-out concordance correlation index
#'
#' @importFrom MASS ginv
#' @keywords internal
#'
#' @noRd
#'
LCCC.FIT_1 <- function(leaveout.list, remain.counts, remain.basis, Y.list, maxim = 1000, ...){
  n <- length(leaveout.list)
  active.set <- (remain.basis < maxim)
  if (tail(remain.counts[active.set], 1) == n - 1){
    remain.counts[ length(remain.counts)] <- n - 2
  }
  feasible.long <- sum(active.set)
  max.long <- max(unlist(lapply(Y.list, length)))
  checks <- matrix(NA, nrow = n, ncol = feasible.long)
  Values <- array(NA, c(max.long, n, feasible.long))


  for (i in 1:n){ ### i = 35
    y <- Y.list[[i]]
    y.long <- length(y)
    grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
    CDFBETA <- GENERATE_BETA_CDF( # a1, a2,
      grid.p, ...)
    NQ <- qnorm(grid.p, 0, 1)
    BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
    BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
    Psi <- cbind(rep(1, length(grid.p)), BETA_BASE_TOTAL_2)

    for (j in 1:feasible.long){ ###  j = 20

      colum_i_ <- leaveout.list[[i]][[1]]
      obs_i_ <- leaveout.list[[i]][[2]] ## length(IncidenceVec_i_)

      SET1_i_ <- sort(colum_i_[  obs_i_ > remain.counts[active.set ][j]   ])
      smPsi_i_ <- Psi[, (SET1_i_) + 1 ] ## dim(smPsi_i_)

      Values[(1:y.long), i, j] <- try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y) ## dim(Qy)
      # checks[i, j] <- length(SET1_i_)
    }
  }
  outputs <- list(Values, checks)
}
