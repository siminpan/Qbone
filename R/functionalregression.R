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

## 2.1 qfrModel ----
#' Fit the quantile functional regression model
#'
#' fFit the quantlet space model in using Markov chain Monte Carlo (MCMC)
#'
#' @param object A Qboneobject
#' @param new.assay.name New assay name assigned to the quantlets data
#' @param data.assay It is the name of the assay whose data will be used
#'  to compute the lasso list. Default is the data from the
#'  defaultAssay(object).
#' @param X Covariates (N by A matrix, N is the number of observations)
#' @param X1 New covariates for inferencing MCMC fit result.
#' @param delta2 Cutoff percentage of the energy for the functional coefficients
#' @param H Number of clustering groups. Cluster basis indices based on their eigen-values.
#' @param pct.range Percentage range from original data set. Use as the range of the domain for the density estimate. Default is c(0.05, 0.95). You can override it by assign \code{xranges = c(min, max)}
#' @param ... Arguments passed to other methods
#'
#' @importFrom pracma gramSchmidt
#'
#' @export
#'
qfrModel <- function(
  object,
  new.assay.name = "Q.F.Regression",
  data.assay = defaultAssay(object),
  X = NULL,
  X1 = NULL,
  delta2 = 0.95,
  H = NULL,
  pct.range = c(0.05, 0.95),
  assay.seed2 = .Random.seed,
  ...
){
  # Check data.assay
  if (data.assay != "Empirical.Coefficients"){
    warning('The default assay is not "Empirical.Coefficients" please double check the defaultAssay() of this Qbone object. This step should be run on results of ecQuantlets().')
  }
  if (is.null(X)){
    message('Covariates X was not provided, will atomically generate covariate matrix based on Qbone object metadata group information. Please double check the covariates matrix X')
    X <- covM(object)
  }
  if (is.null(H)){
    message('Number of clustering groups H was not provided, will atomically generate covariate matrix based on Qbone object metadata group information. Please double check H.')
    H <- length(unique(object@meta.data[["group"]]))
  }
  # Get data
  empCoefs.list <- getQboneData(object, slot = 'data', assay = data.assay)
  empCoefs <- matrix(unlist(empCoefs.list, use.names = F), nrow=length(empCoefs.list), byrow=TRUE)
  reduceBasis <- object@assays[[data.assay]]@scale.data[["reduceBasis"]]
  orthogBasis <- object@assays[[data.assay]]@scale.data[["orthogBasis"]]
  quantlets <- object@assays[[data.assay]]@scale.data[["quantlets"]]
  # preMCMC
  emp_fit <- preMCMC(empCoefs = empCoefs,
                     reduceBasis = reduceBasis,
                     orthogBasis = orthogBasis,
                     quantlets = quantlets,
                     X = X,
                     delta2 = delta2,
                     H = H,
                     ...)
  # MCMC
  # mcmc_fit <- MCMC(X, emp_fit$sd_l2, emp_fit$cluster, emp_fit$TB00, ...)
  mcmc_fit <- MCMC(X = X,
                   empCoefs = emp_fit$sd_l2,
                   r = emp_fit$cluster,
                   TB00 = emp_fit$TB00,
                   ...
                   )
  # MCMC inference
  # orig.dataset <- getQboneData(object, slot = 'data', assay = object@assays[[data.assay]]@assay.orig)
  if ("thin" %in% names(object@assays[[data.assay]]@scale.data)){
    orig.dataset <- getQboneData(object, slot = 'data', assay = object@assays[["Thin"]]@assay.name)
  } else {
    orig.dataset <- getQboneData(object, slot = 'data', assay = object@assays[[1]]@assay.name)
  }
  if (hasArg(xranges)){
    xranges1 = xranges
  } else {
    # xranges1 = c(min(sapply(orig.dataset, min)), max(sapply(orig.dataset, max)))
    xranges1 = quantile(sapply(orig.dataset, quantile, probs = seq(0, 1, 1/10)), probs = pct.range)
  }
  message("Inferencing MCMC model, this may take a while.")
  mcmc_infer <- inferenceMCMC(mcmcEmpCoef = mcmc_fit[[1]],
                              BackTransfor = emp_fit$sdPhi,
                              X = X,
                              X1 = X1,
                              p = object@assays[[data.assay]]@scale.data[["p"]],
                              xranges = xranges1,
                              ...)
  new.qbonedata <- object@assays[[data.assay]]
  new.qbonedata@assay.orig <- new.qbonedata@assay.name
  new.qbonedata@assay.name <- new.assay.name
  new.qbonedata@scale.data <- append(object@assays[[data.assay]]@scale.data,
                                     list(mcmc_infer = mcmc_infer)
  )
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

## 6.1 preMCMC ----
#' Getting data ready for Markov chain Monte Carlo (MCMC).
#'
#' Computes scaled empirical coefficients, clustering, and initial
#' coefficients input
#'
#' @param empCoefs empirical coefficients
#' @param reduceBasis  reduced basis function
#' @param orthogBasis orthogonal basis function
#' @param quantlets quantlets basis function
#' @param X Covariates
#' @param delta2 Cutoff percentage of the energy for the functional coefficients
#' @param H Number of clustering groups.
#' @param ... Arguments passed to other methods
#'
#' @return matrix containing # of beta parameters times the length of index.p
#'
#' @keywords internal
#'
#' @noRd
#'
preMCMC <- function(
  empCoefs,
  reduceBasis,
  orthogBasis,
  quantlets,
  X,
  delta2,
  H,
  ...){
  if (dim(reduceBasis)[1] == dim(quantlets)[1]){
    reduceBasis_s <- reduceBasis
  } else {
    biggrids.long <- dim(quantlets)[1]

    # Assign random seed
    .Random.seed <- assay.seed2

    biggrids <- sort(sample(dim(reduceBasis)[1], biggrids.long), method = "quick")
    reduceBasis_s <- reduceBasis[biggrids, ]
  }


  # upper_R_T <- t(reduceBasis_s) %*% quantlets
  upper_R_T <- eigenmapmtm(reduceBasis_s, quantlets) # dim(lambda_v)  ## each column lambda_v
  singluar <- round(diag(upper_R_T), 8) ## upper_R_T%*%t(upper_R_T)
  est_eigen <- c(200, singluar)

  # r <- basisCluster(reduceBasis, orthogBasis, H) # reduceBasis = reduceBasis, orthogBasis = orthogBasis, H = H
  r <- basisCluster(reduceBasis = reduceBasis,
                    orthogBasis = orthogBasis,
                    H = H,
                    ...)

  d_l2 <- empCoefs
  K <- dim(d_l2)[2]
  Tp <- dim(quantlets)[1]

  Phi <- cbind(rep(1, Tp), quantlets)

  # scale_f2 <- diag(t(d_l2) %*% d_l2)^{
  #   -0.5
  # }
  scale_f <- diag(eigenmapmtm(d_l2, d_l2))^{-0.5}

  # sd_l22 <- d_l2 %*% diag(scale_f, K, K)
  sd_l2 <- eigenmapmm(d_l2, diag(scale_f, K, K)) ## diag( t(sd_l2)%*%sd_l2 )

  # sdPhi2 <- Phi %*% diag(scale_f^{
  #   -1
  # }, K, K)
  sdPhi <- eigenmapmm(Phi, diag(scale_f^{-1}, K, K))

  Px <- dim(X)[2]

  B00 <- solve(t(X) %*% X) %*% t(X) %*% sd_l2
  beta_LSE <- B00
  ## 1 -> smallest!  64 -> largest!

  ### ORD_beta_LSE = matrix(NA, ncol= Px, nrow= K)   ## j=1

  ### for(j in 1:Px){ORD_beta_LSE[, j] = order(abs(beta_LSE[j,]))}

  JP_XI <- matrix(NA, nrow = Px, ncol = K)
  JP_XE <- matrix(NA, nrow = Px, ncol = K) ##  d_0_2[i,][JP_I[i,]]
  JP_XB <- matrix(NA, nrow = Px, ncol = K)
  B00_2 <- B00^2
  for (i in 1:Px){ ##  i =1
    JP_XI[i, ] <- c(1, 2, (order(B00_2[i, -c(1, 2)], decreasing = TRUE) + 2))
    ## JP_XI[i,] = c(1, (order(B00_2[i, -c(1)], decreasing = TRUE) + 1))
    JP_XE[i, ] <- cumsum(B00_2[i, JP_XI[i, ]]) / sum(B00_2[i, JP_XI[i, ]])
    JP_XB[i, JP_XI[i, ]] <- JP_XE[i, ]
  }

  if (Px == length(delta2)){
    cuts <- delta2
  } else {
    cuts <- rep(delta2, Px)
  }

  Set.off <- vector("list", Px)
  Set.on <- vector("list", Px)
  for (i in 1:Px){ ## i=3
    Set.off[[i]] <- seq(K)[JP_XB[i, ] > cuts[i]]
    Set.on[[i]] <- seq(K)[JP_XB[i, ] <= cuts[i]]
  }

  TB00 <- matrix(0, nrow = Px, ncol = K)
  lambda_TB00 <- matrix(0, nrow = K, ncol = Px)
  for (i in 1:Px){ ## i=3
    if (is.null(Set.off[[i]])){
      lambda_TB00[, i] <- 0
    } else {
      lambda_TB00[  Set.off[[i]], i] <- 100000
      lambda_TB00[ -Set.off[[i]], i] <- 0
    }
    TB00[i, ] <- ifelse(lambda_TB00[, i] == 0, B00[i, ], 0)
  }

  outputs <- list(sd_l2, sdPhi, scale_f, d_l2, Phi, TB00, r)
  names(outputs) <- list("sd_l2", "sdPhi", "scale_f", "d_l2", "Phi", "TB00", "cluster")
  return(outputs)
}

## 6.2 covM ----
#' Generate covariate matrix from Qbone object atomically.
#'
#' atomically generate covariate matrix based on Qbone object metadata group information.
#'
#' @param object Qbone object
#' @param added.cov probability grids on (0,1)
#' @param ... Arguments passed to other methods
#'
#' @return matrix containing # of beta parameters times the length of index.p
#'
#' @keywords internal
#'
#' @noRd
#'
covM <- function(
  object,
  added.cov = NULL,
  ...
  ){
  # if (int == T){
    int <- rep(1, length(object@meta.data[["group"]]))
  # } else {
  #   int <- rep(0, length(object@meta.data[["group"]]))
  # }
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
  X <- X[,-1]
  return(X)
}

## 6.3 basisCluster ----
#' Cluster basis indices based on their eigen-values
#'
#' @param reduceBasis  reduced basis function
#' @param orthogBasis orthogonal basis function
#' @param H Number of clustering groups.
#' @param ... Arguments passed to other methods
#'
#' @return matrix containing # of beta parameters times the length of index.p
#'
#' @importFrom stats hclust cutree
#' @keywords internal
#'
#' @noRd
#'
basisCluster <- function(
  reduceBasis,
  orthogBasis,
  H = NULL,
  ...
  ){
  lambda_v <- eigenmapmtm(reduceBasis, orthogBasis)
  est_eigen <- round(diag(eigenmapmtm(lambda_v, lambda_v)), 8)
  if (length(est_eigen) >= 6){
    H1 <- H - 1
    hc <- hclust(dist(est_eigen)^2, "cen")
    r <- cutree(hc, k = H1)
  } else {
    r <- seq(length(est_eigen))
  }
  return(c(1, r + 1))
}

## 6.4 MCMC ----
#' Markov chain Monte Carlo (MCMC) computations
#'
#' @param X Covariates
#' @param empCoefs Empirical coefficients
#' @param r Clustering indicators
#' @param TB00  initial values of functional coefficients
#' @param noi number of MCMC iteration
#' @param burn Number of burn-in
#' @param ep1 prior for nu. Default is inverse gamma dist
#' @param keep.zero T or F to keep 0
#' @param ... Arguments passed to other methods
#'
#' @return MCMC samples
#'
#' @importFrom stats hclust cutree
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#'
#' @noRd
#'
MCMC <- function(
  X,
  empCoefs,
  r,
  TB00,
  noi = 2000,
  burn = 200,
  ep1 = 0.0064,
  keep.zero = FALSE,
  ...
){
  # Set up
  N <- dim(X)[1]
  Px <- dim(X)[2]
  K <- dim(empCoefs)[2]

  H <- length(unique(r))
  B00 <- solve(t(X) %*% X) %*% eigenmapmtm(X, empCoefs)
  A00 <- TB00
  A00 <- ifelse(A00 == 0, 0, 1)
  C00 <- matrix(rep(r, Px), nrow = Px, ncol = K, byrow = TRUE)
  vec_C00 <- as.vector(C00)
  PA00 <- matrix(NA, ncol = K, nrow = Px)
  vec_PA00 <- as.vector(PA00)

  SSE_K <- diag(eigenmapmtm((empCoefs - X %*% B00), (empCoefs - X %*% B00)))
  sigma2_K <- SSE_K/N

  nu0 <- ep1
  g <- rep(100000, K)

  one <- 1

  # Fully Bayesian

  set.seed(noi + 2022) # || double check

  MCMC_BETA <- matrix(0, ncol = Px * K, nrow = noi)
  MCMC_OMEGA <- matrix(0, ncol = K, nrow = noi)

  MCMC_ALPHA <- matrix(0, ncol = Px * K, nrow = noi)


  MCMC_PHI <- matrix(NA, ncol = length(unique(r)), nrow = burn)
  MCMC_GAM <- matrix(NA, ncol = length(unique(r)), nrow = burn)

  MCMC_PA <- matrix(0, ncol = Px * K, nrow = burn)
  MCMC_TAU <- matrix(0, ncol = Px * K, nrow = burn)
  PA00 <- matrix(NA, ncol = K, nrow = Px)
  TAU00 <- matrix(NA, ncol = K, nrow = Px)
  message("Start MCMC iteration: ")
  pb <- txtProgressBar(min = 0, max = noi, style = 3)
  for (it in 1:noi){

    if (it <= burn){

      if (it == 1){
        sigma2_K0 <- sigma2_K
        temp.B00 <- B00
        MCMC_BETA[ 1, ] <- as.vector(B00)
        temp.A00 <- A00
        temp.C00 <- C00
        MCMC_PA[ 1, ] <- as.vector(PA00)
        temp.PA00 <- PA00
        temp.TAU00 <- TAU00
      } else {
        sigma2_K0 <- MCMC_OMEGA[(it - 1), ]
        temp.B00 <- matrix(MCMC_BETA[(it - 1), ], ncol = K, nrow = Px)
        temp.A00 <- matrix(MCMC_ALPHA[(it - 1), ], ncol = K, nrow = Px)
      }
      # Empirical Bayes part

      for (j in 1:Px){

        if ((Px - 1) != 1){
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1))
          X_j_ <- X[, -j]
        } else {
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1), byrow = TRUE)
          X_j <- matrix(X[, -j], ncol = (Px - 1), nrow = N, byrow = TRUE)
        }

        B_j <- temp.B00[  j, ]
        X_j <- X[, j]

        V_j_all_h <- (sum(X_j^2)/sigma2_K0)^{-1}
        eta_j_all_h <- B_j/sqrt(V_j_all_h)

        # cluster by engen  ##indicator by trashold by energy
        cluster_j <- C00[j, ]
        nonzeor_j <- A00[j, ]

        Ob_j_all_h <- rep(NA, K)

        for (h in 1:length(unique(r))){
          hh <- unique(r)[h]

          crruent_cluster_j <- (cluster_j == hh)

          nonzero_in_crruent_cluster_j <- (crruent_cluster_j) & (nonzeor_j == 1)
          zero_in_crruent_cluster_j <- (crruent_cluster_j) & (nonzeor_j == 0)

          total.obs_in_crruent_cluster_j <- sum(crruent_cluster_j)
          effective.obs_in_crruent_cluster_j <- sum(nonzero_in_crruent_cluster_j)
          prob_nonzero_in_crruent_cluster_j <- effective.obs_in_crruent_cluster_j / total.obs_in_crruent_cluster_j

          freedom_in_crruent_cluster_j <- effective.obs_in_crruent_cluster_j - 1

          MCMC_PHI[it, h] <- prob_nonzero_in_crruent_cluster_j

          if (prob_nonzero_in_crruent_cluster_j == 1){
            nonzeor_j[crruent_cluster_j] <- 1
            temp.TAU00[j, crruent_cluster_j] <- Inf
            MCMC_GAM[it, h] <- Inf
            temp.PA00[j, crruent_cluster_j] <- 1
          } else if (prob_nonzero_in_crruent_cluster_j == 0){
            if (keep.zero == TRUE){
              nonzeor_j[ crruent_cluster_j] <- 0
              temp.TAU00[j, crruent_cluster_j ] <- 0
              MCMC_GAM[it, h] <- 0
              temp.PA00[j, crruent_cluster_j ] <- 0
            } else {
              nonzeor_j[ crruent_cluster_j] <- 1
              temp.TAU00[j, crruent_cluster_j ] <- Inf
              MCMC_GAM[it, h] <- Inf
              temp.PA00[j, crruent_cluster_j ] <- 1
            }
          } else {
            if (freedom_in_crruent_cluster_j >= 1){
              G_j_h <- max(0, sum((eta_j_all_h^2)[nonzero_in_crruent_cluster_j]) / (sum(nonzero_in_crruent_cluster_j) - 1))
              MCMC_GAM[it, h] <- G_j_h

              temp.TAU00[j, crruent_cluster_j ] <- V_j_all_h[crruent_cluster_j] * G_j_h
              prob_j_h <- prob_nonzero_in_crruent_cluster_j

              if (sum(nonzero_in_crruent_cluster_j) >= 1){
                prob0 <- prob_j_h / (1 - prob_j_h) / sqrt(1 + G_j_h) * exp(0.5 * eta_j_all_h[nonzero_in_crruent_cluster_j]^2 * (G_j_h / (1 + G_j_h)))
                Ob_j_all_h[nonzero_in_crruent_cluster_j] <- ifelse(prob0 != Inf, prob0, 10^128)
                temp.PA00[j, nonzero_in_crruent_cluster_j] <- Ob_j_all_h[nonzero_in_crruent_cluster_j] / (Ob_j_all_h[nonzero_in_crruent_cluster_j] + 1)
              }
              if (sum(zero_in_crruent_cluster_j) >= 1){
                Ob_j_all_h[zero_in_crruent_cluster_j] <- prob_j_h / (1 - prob_j_h) / sqrt(1 + G_j_h) * exp(0.5 * eta_j_all_h[zero_in_crruent_cluster_j]^2 * (G_j_h / (1 + G_j_h)))
                temp.PA00[j, zero_in_crruent_cluster_j] <- Ob_j_all_h[zero_in_crruent_cluster_j] / (Ob_j_all_h[zero_in_crruent_cluster_j] + 1)
              }
            }
            if (freedom_in_crruent_cluster_j == 0){
              G_j_h <- Inf
              MCMC_GAM[it, h] <- G_j_h
              temp.TAU00[j, crruent_cluster_j] <- V_j_all_h[crruent_cluster_j] * G_j_h
              prob_j_h <- prob_nonzero_in_crruent_cluster_j

              if (sum(nonzero_in_crruent_cluster_j) >= 1){
                Ob_j_all_h[  nonzero_in_crruent_cluster_j] <- Inf
                temp.PA00[j, nonzero_in_crruent_cluster_j] <- 1
              }
              if (sum(zero_in_crruent_cluster_j) >= 1){
                Ob_j_all_h[  zero_in_crruent_cluster_j] <- Inf
                temp.PA00[j, zero_in_crruent_cluster_j] <- 1
              }
            }
          } ## close probability condition
        }
      }

      MCMC_TAU[it, ] <- as.vector(temp.TAU00)
      MCMC_PA [it, ] <- as.vector(temp.PA00)
    }

    # Main Estimation

    for (j in 1:Px){

      for (k in 1:K){

        if ((Px - 1) != 1){
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1))
          X_j_ <- X[, -j]
        } else {
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1), byrow = TRUE)
          X_j <- matrix(X[, -j], ncol = (Px - 1), nrow = N, byrow = TRUE)
          X_j_ <- 1 # || double check
        }

        B_j <- temp.B00[  j, ]
        X_j <- X[, j]

        hat_xSx <- sum(X_j^2 / sigma2_K0[k])
        hat_xSy1 <- sum(X_j * empCoefs[, k] / sigma2_K0[k])

        hat_xSy2 <- sum(X_j * eigenmapmm(X_j_, B_j_[, k]) / sigma2_K0[k]) #### t(X_j)%*%X_j_%*% B_j_[, k]/sigma2_K0[k]
        hat_xSy <- hat_xSy1 - hat_xSy2

        if (temp.TAU00[j, k] <= 0.000000001){
          Vb <- 0
          Eb <- eigenmapmm(Vb, (hat_xSy))
          E <- matrix(rnorm(one, 0, 1), nrow = 1, ncol = one) # || double check seed
          beta_j <- t(c(Eb))
        } else {
          Vb <- g[j] * solve(hat_xSx + 1 / temp.TAU00[j, k]) / (g[j] + 1)
          Eb <- eigenmapmm(Vb, (hat_xSy))
          E <- matrix(rnorm(one, 0, 1), nrow = 1, ncol = one) # || double check seed
          beta_j <- t(t(E %*% chol(Vb)) + c(Eb))
        } ##

        u <- runif(1, 0, 1)  # || double check seed
        if (u <= temp.PA00[j, k]){
          beta_j <- beta_j
        } else {
          beta_j <- 0
        }
        ### posit =   (k-1)*Px + j
        ### MCMC_BETA[ it,  posit ] = beta_j
        temp.B00[j, k ] <- beta_j
      } ##   close ## j =1
    } ##  close ## k =1
    MCMC_BETA[it, ] <- as.vector(temp.B00)
    temp.A00 <- temp.B00
    temp.A00 <- ifelse(temp.A00 == 0, 0, 1)
    MCMC_ALPHA[it, ] <- as.vector(temp.A00)

    est <- matrix(MCMC_BETA[it, ], nrow = Px, ncol = K)
    SSE <- diag(eigenmapmtm((empCoefs - X %*% est), (empCoefs - X %*% est)))
    MCMC_OMEGA[it, ] <- 1 / rgamma(K, (nu0 + N) / 2, (nu0 + SSE) / 2) # || double check seed
    setTxtProgressBar(pb, it)
  }
  close(pb)
  #

  POST_BETA_FULL <- MCMC_BETA[-c(1:burn), ]
  POST_OMEGA_FULL <- MCMC_OMEGA[-c(1:burn), ]
  POST_ALPHA <- MCMC_ALPHA[-c(1:burn), ]
  POST_TAU <- MCMC_TAU[-c(1:burn), ]

  outputs <- list(POST_BETA_FULL)
  names(outputs) <- list("postbeta")
  return(outputs)
}

## 6.5 inferenceMCMC ----
#' Inference in data space based on the output produced from MCMC function
#'
#' @param mcmcEmpCoef MCMC samples for coefficients. Output from \code{MCMC()}
#' @param BackTransfor P by K matrix. Basis function to transform back to the
#' original data space (P: desired number of probability grids, and K: number
#' of basis)
#' @param X Covariates.
#' @param signifit Scalar in (0,1)	The level of size alpha (default=0.975
#' allow for the inference at the level of .05)
#' @param X1 New covariates
#' @param p Probability grids
#' @param n.sup Positive integer. Number of grids points of density estimates
#' @param xranges Vector consisting of c(min, max) as the range of the domain
#' for the density estimate
#' @param ... Arguments passed to other methods
#'
#' @return list of MCMC inference result:
#'              est       Estimations of intercept and covariates in the
#'              projected space.
#'              estCIu		Upper point confidence intervals for the intercept
#'              and covariates in the data space
#'              estCIl		Lower point confidence intervals for the intercept
#'              and covariates in the data space
#'              DataEst   Estimations of intercept and covariates in the data
#'              space.
#'              jointCI   jointCI[,,1] includes upper joint confidence
#'              intervals for the intercept and covariates in the data space
#'              whereas jointCI[,,2] includes lower joint confidence intervals
#'              for the intercept and covariates in the data space.
#'              local_p	  SimBas for covariates (logical values where true
#'              means significant and false means insignificant for the given
#'              level of size alpha).
#'              global_p	Posterior probability scores for covariates
#'              mu_G		  Point estimation of mu for new covariates X1 is given
#'              in the first row where the second and third rows give the upper
#'              and lower confidence intervals.
#'              sigma_G	  Point estimation of variance for new covariates X1 is
#'              given in the first row where the second and third rows give the
#'              upper and lower confidence intervals.
#'              mu3_G	    Point estimation of skewness for new covariates X1 is
#'              given in the first row where the second and third rows give the
#'              upper and lower confidence intervals.
#'              mu4_G     Point estimation of kurtosis for new covariates X1 is
#'              given in the first row where the second and third rows give the
#'              upper and lower confidence intervals.
#'              mu_diff    Results for the mean different testing for two
#'              consecutive subjects
#'              (8th column means  posterior probability scores).
#'              sigma_diff  Results for the variance different testing for two
#'              consecutive subjects
#'              (8th column means  posterior probability scores).
#'              mu3_diff    Results for the  skewness  different testing for two
#'              consecutive subjects
#'              (8th column means posterior probability scores).
#'              mu4_diff    Results for the  kurtosis  different testing for two
#'              consecutive subjects
#'              (8th column means posterior probability scores).
#'              den_G       Density estimates for the new covariates.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#'
#' @noRd
#'
inferenceMCMC <- function(
  mcmcEmpCoef,
  BackTransfor = NULL,
  X,
  signifit = 0.975,
  X1 = NULL,
  p,
  n.sup = 100,
  xranges = c(-10, 60),
  ...
  ){
  BN <- dim(mcmcEmpCoef)[1]
  Px <- dim(X)[2]
  Tp <- dim(BackTransfor)[1]
  K <- dim(BackTransfor)[2]

  ITF_BETA_FULL <- array(NA, c(BN, Tp, dim(X)[2]))

  for (i in 1:(BN)){
    ITQ_BAYES <- eigenmapmmt(BackTransfor, matrix(mcmcEmpCoef[i, ], ncol = K, nrow = Px))
    for (j in 1:Px){
      ITF_BETA_FULL[i, , j] <- ITQ_BAYES[, j]
    }
  }

  MEAN_BETA_FULL <- apply(mcmcEmpCoef, 2, mean)
  BETA_FULL <- matrix(MEAN_BETA_FULL, nrow = Px)

  QUAN_BETA_FULL <- apply(ITF_BETA_FULL, c(2, 3), quantile, probs = c(1 - signifit, signifit))

  QUAN_ucl_BETA_FULL <- QUAN_BETA_FULL[2, , ]
  QUAN_lcl_BETA_FULL <- QUAN_BETA_FULL[1, , ]

  DataSp_est <- eigenmapmmt(BackTransfor, BETA_FULL)

  MCMC_F_BASIS_MEAN <- matrix(0, nrow = Tp, ncol = Px)
  MCMC_F_BASIS_VAR <- matrix(0, nrow = Tp, ncol = Px)
  MCMC_F_BASIS_CI <- array(0, c(Tp, Px, 2))
  MCMC_F_BASIS_Z <- matrix(0, nrow = BN, ncol = Px)

  for (j in 1:Px){
    M_j <- apply(ITF_BETA_FULL[, , j], 2, mean)
    S_j <- (apply(ITF_BETA_FULL[, , j], 2, var))^{0.5}
    Z_j <- eigenmapmm((ITF_BETA_FULL[, , j] - matrix(M_j, nrow = BN, ncol = Tp, byrow = TRUE)), diag(S_j^{-1}, ncol = Tp, nrow = Tp))
    abs_Z_J <- apply(abs(Z_j), 1, max)
    qalpha <- quantile(abs_Z_J, signifit - 0.025)
    I_j_995 <- M_j + qalpha * S_j
    I_j_005 <- M_j - qalpha * S_j
    MCMC_F_BASIS_MEAN[, j] <- M_j
    MCMC_F_BASIS_VAR[, j] <- S_j
    MCMC_F_BASIS_Z[, j] <- abs_Z_J

    MCMC_F_BASIS_CI[, j, 1 ] <- I_j_995
    MCMC_F_BASIS_CI[, j, 2 ] <- I_j_005
  }

  #
  PMAP_FULL0 <- matrix(0, nrow = Tp, ncol = Px - 1)
  for (j in 2:Px){
    M_j <- apply(ITF_BETA_FULL[, , j], 2, mean)
    S_j <- (apply(ITF_BETA_FULL[, , j], 2, var))^{0.5}
    for (k in 1:BN){
      PMAP_FULL0[, (j - 1)] <- PMAP_FULL0[, (j - 1)] + (abs(M_j / S_j) <= MCMC_F_BASIS_Z[k, j]) + 0
    }
  }
  PMAP_FULL <- PMAP_FULL0 / BN

  LPMAP_FULL <- apply(PMAP_FULL, 2, function(x) (x <= (1 - signifit) * 2))

  GPMAP_FULL <- apply(PMAP_FULL, 2, min)
  #
  if (!is.null(X1)){
    PX0 <- X1

    STAT_MU <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    for (i in 1:BN){
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px){
        ITQ_BAYES0[, j] <- ITF_BETA_FULL[i, , j]
      }
      ITQ_BAYES <- eigenmapmmt(ITQ_BAYES0, PX0)
      STAT_MU[i, ] <- eigenmapmtm(ITQ_BAYES, rep(1 / Tp, Tp))
    }

    MU_BAYES <- rbind(apply(STAT_MU, 2, mean), apply(STAT_MU, 2, quantile, probs = c(1 - signifit, signifit)))

    STAT_VAR <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    for (i in 1:BN){
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px){
        ITQ_BAYES0[, j] <- ITF_BETA_FULL[i, , j]
      }

      ITQ_BAYES <- eigenmapmmt(ITQ_BAYES0, PX0)
      CT_ITQ_BAYES <- ITQ_BAYES - matrix(rep(MU_BAYES[1, ], each = Tp, nrow = Tp), ncol = dim(PX0)[1])
      STAT_VAR[i, ] <- (eigenmapmtm(CT_ITQ_BAYES^2, rep(1 / Tp, Tp)))^{0.5}
    }

    VAR_BAYES <- rbind(apply(STAT_VAR, 2, mean), apply(STAT_VAR, 2, quantile, probs = c(1 - signifit, signifit)))

    STAT_MU3 <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    STAT_MU4 <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    for (i in 1:BN){
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px){
        ITQ_BAYES0[, j] <- ITF_BETA_FULL[i, , j]
      }
      ITQ_BAYES <- eigenmapmmt(ITQ_BAYES0, PX0)

      CT_ITQ_BAYES <- ITQ_BAYES - matrix(rep(MU_BAYES[1, ], each = Tp), nrow = Tp, ncol = dim(PX0)[1])
      ST_ITQ_BAYES <- eigenmapmm((ITQ_BAYES - matrix(rep(MU_BAYES[1, ], each = Tp), nrow = Tp, ncol = dim(PX0)[1])), diag(VAR_BAYES[1, ]^{-1}, ncol = dim(PX0)[1], nrow = dim(PX0)[1]))
      STAT_MU3[i, ] <- eigenmapmtm((ST_ITQ_BAYES^3), rep(1 / Tp, Tp))
      STAT_MU4[i, ] <- eigenmapmtm((CT_ITQ_BAYES^4), rep(1 / Tp, Tp)) / VAR_BAYES[1, ]^4
    }

    MU3_BAYES <- rbind(apply(STAT_MU3, 2, mean), apply(STAT_MU3, 2, quantile, probs = c(1 - signifit, signifit)))

    MU4_BAYES <- rbind(apply(STAT_MU4, 2, mean), apply(STAT_MU4, 2, quantile, probs = c(1 - signifit, signifit)))

    #

    spx <- dim(PX0)[1]
    STAT_MU_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)){
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1
      TEST_STATS <- STAT_MU[, j ] - STAT_MU[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{0.5}
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_MU_DIFF[k.n, 1] <- i
      STAT_MU_DIFF[k.n, 2] <- j
      STAT_MU_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_MU_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0), 0, 1)
      STAT_MU_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_MU_DIFF[k.n, 9] <- ifelse(STAT_MU_DIFF[k.n, 3] == 0, 0, 1)
      STAT_MU_DIFF[k.n, 10] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }

    STAT_VAR_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)){
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1
      TEST_STATS <- STAT_VAR[, j] - STAT_VAR[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{0.5}
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_VAR_DIFF[k.n, 1] <- i
      STAT_VAR_DIFF[k.n, 2] <- j
      STAT_VAR_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_VAR_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0), 0, 1)
      STAT_VAR_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_VAR_DIFF[k.n, 9] <- ifelse(STAT_VAR_DIFF[k.n, 3] == 0, 0, 1)
      STAT_VAR_DIFF[k.n, 10 ] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }

    STAT_MU3_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)){
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1
      TEST_STATS <- STAT_MU3[, j] - STAT_MU3[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{0.5}
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_MU3_DIFF[k.n, 1] <- i
      STAT_MU3_DIFF[k.n, 2] <- j
      STAT_MU3_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_MU3_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0), 0, 1)
      STAT_MU3_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_MU3_DIFF[k.n, 9] <- ifelse(STAT_MU3_DIFF[k.n, 3] == 0, 0, 1)
      STAT_MU3_DIFF[k.n, 10 ] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }

    STAT_MU4_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)){
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1

      TEST_STATS <- STAT_MU4[, j] - STAT_MU4[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{0.5}
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_MU4_DIFF[k.n, 1] <- i
      STAT_MU4_DIFF[k.n, 2] <- j
      STAT_MU4_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_MU4_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0), 0, 1)
      STAT_MU4_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_MU4_DIFF[k.n, 9] <- ifelse(STAT_MU4_DIFF[k.n, 3] == 0, 0, 1)
      STAT_MU4_DIFF[k.n, 10] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }
    #

    xdomain <- seq(xranges[1], xranges[2], length.out = n.sup)
    I_DENSITY_FULL <- array(NA, c(BN, n.sup, dim(X1)[1]))
    DENSITY_FULL <- array(NA, c(BN, (n.sup - 1), dim(X1)[1]))

    for (i in 1:BN){
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px){
        ITQ_BAYES0[, j] <- ITF_BETA_FULL[i, , j]
      }
      ITQ_BAYES <- eigenmapmmt(ITQ_BAYES0, X1)
      for (h in 1:dim(X1)[1]){
        for (k in 1:n.sup){
          if (sum(ITQ_BAYES[, h] <= xdomain[k]) == 0){
            I_DENSITY_FULL[i, k, h] <- 0
          } else {
            I_DENSITY_FULL[i, k, h] <- max(p[ITQ_BAYES[, h] <= xdomain[k]])
          }
          DENSITY_FULL[i, , h] <- diff(I_DENSITY_FULL[i, , h], 1) / diff(xdomain, 1)
        }
      }
    }
    MEAN_DENS_FULL <- apply(DENSITY_FULL, c(2, 3), mean)
    QUNT_DENS_FULL <- apply(DENSITY_FULL, c(2, 3), quantile, probs = c(1 - signifit, signifit))

    #
  } else {
    MU_BAYES <- NULL
    VAR_BAYES <- NULL
    MU3_BAYES <- NULL
    MU4_BAYES <- NULL
    STAT_MU_DIFF <- NULL
    STAT_VAR_DIFF <- NULL
    STAT_MU3_DIFF <- NULL
    STAT_MU4_DIFF <- NULL
    QUNT_DENS_FULL <- NULL
    MEAN_DENS_FULL <- NULL
  }
  xdomain = seq(xranges[1], xranges[2], length.out = n.sup)

  outputs <- list(
    BETA_FULL, QUAN_ucl_BETA_FULL, QUAN_lcl_BETA_FULL, DataSp_est,
    MCMC_F_BASIS_CI, LPMAP_FULL, GPMAP_FULL, MU_BAYES, VAR_BAYES, MU3_BAYES, MU4_BAYES,
    STAT_MU_DIFF, STAT_VAR_DIFF, STAT_MU3_DIFF, STAT_MU4_DIFF, QUNT_DENS_FULL, MEAN_DENS_FULL,
    X, X1, xdomain
  )

  names(outputs) <- list(
    "est", "estCIu", "estCIl", "DataEst",
    "jointCI", "local_p", "global_p", "mu_G", "sigma_G", "mu3_G", "mu4_G",
    "mu_diff", "sigma_diff", "mu3_diff", "mu4_diff", "denCI_G", "den_G",
    "covariates.X", "covariates.X1", "xdomain"
  )

  return(outputs)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
