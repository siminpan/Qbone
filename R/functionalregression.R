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
#' @param delta2 Cutoff percentage of the energy for the functional coefficients
#' @param H Number of clustering groups. Cluster basis indices based on their eigen-values.
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
  delta2 = 0.95,
  H = NULL,
  assay.seed2 = .Random.seed,
  ...
){
  # Check data.assay
  if(data.assay != "Empirical.Coefficients"){
    warning('The default assay is not "Empirical.Coefficients" please double the defaultAssay() of this Qbone object. This step should be run on results of ecQuantlets().')
  }
  empCoefs.list <- getQboneData(object, slot = 'data', assay = data.assay)
  empCoefs <- matrix(unlist(empCoefs.list, use.names = F), nrow=length(empCoefs.list), byrow=TRUE)
  reduceBasis <- object@assays[[data.assay]]@scale.data[["reduceBasis"]]
  quantlets <- object@assays[[data.assay]]@scale.data[["quantlets"]]
  emp_fit <- preMCMC(empCoefs, reduceBasis, quantlets, # X = X, delta2 = delta2, H = H,
                     ...)

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
#' @param empCoefs vector containing sequence of beta parameter
#' @param reduceBasis  vector containing sequence of beta parameter
#' @param quantlets probability grids on (0,1)
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

  if (length(est_eigen) >= 7){
    H1 <- H
    hc <- hclust(dist(est_eigen)^2, "cen")
    r <- cutree(hc, k = H1)
  } else {
    H1 <- H
    r <- seq(length(est_eigen))
  }

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
  BETA_LSE <- B00
  ## 1 -> smallest!  64 -> largest!

  ### ORD_BETA_LSE = matrix(NA, ncol= Px, nrow= K)   ## j=1

  ### for(j in 1:Px){ORD_BETA_LSE[, j] = order(abs(BETA_LSE[j,]))}

  JP_XI <- matrix(NA, nrow = Px, ncol = K)
  JP_XE <- matrix(NA, nrow = Px, ncol = K) ##  d_0_2[i,][JP_I[i,]]
  JP_XB <- matrix(NA, nrow = Px, ncol = K)
  B00_2 <- B00^2
  for (i in 1:Px){ ##  i =1
    JP_XI[i, ] <- c(1, 2, (order(B00_2[i, -c(1, 2) ], decreasing = TRUE) + 2))
    ## JP_XI[i,] =   c(1,   (order(    B00_2[i, -c(1) ] , decreasing = TRUE ) + 1 )  )
    JP_XE[i, ] <- cumsum(B00_2[i, JP_XI[i, ]]) / sum(B00_2[i, JP_XI[i, ]])
    JP_XB[i, JP_XI[i, ]  ] <- JP_XE[i, ]
  }

  if (Px != length(delta2)){
    cuts <- rep(delta2, Px)
  }
  if (Px == length(delta2)){
    cuts <- delta2
  }
  # cuts = c( 0.667, 0.667,  0.675, 0.652, 0.667, 0.652, 0.65  )
  Set.off <- vector("list", Px)
  Set.on <- vector("list", Px)
  for (i in 1:Px){ ## i=3
    Set.off[[i]] <- seq(K)[ JP_XB[i, ] > cuts[i]   ]
    Set.on[[i]] <- seq(K)[ JP_XB[i, ] <= cuts[i]   ]
  }



  TB00 <- matrix(0, nrow = Px, ncol = K)
  Lambda_TB00 <- matrix(0, nrow = K, ncol = Px)
  for (i in 1:Px){ ## i=3
    if (is.null(Set.off[[i]]) == TRUE){
      Lambda_TB00[, i] <- 0
    }
    if (is.null(Set.off[[i]]) != TRUE){
      Lambda_TB00[  Set.off[[i]], i] <- 100000
      Lambda_TB00[ -Set.off[[i]], i] <- 0
    }
    TB00[i, ] <- ifelse((Lambda_TB00[, i] == 0), B00[i, ], 0)
  }

  outputs <- list(sd_l2, sdPhi, scale_f, d_l2, Phi, TB00, r)
  names(outputs) <- list("sd_l2", "sdPhi", "scale_f", "d_l2", "Phi", "TB00", "cluster")
  return(outputs)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
