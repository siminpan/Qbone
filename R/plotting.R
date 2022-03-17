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

## 2.1 nllplot ----
#' Plot near-losslessness parameters.
#'
#' Plots the near-losslessness parameters minimum concordance (\eqn{\rho^{o}})
#' and average (\eqn{\bar{\rho}}) against the number of basis coefficients
#' function of \eqn{K_{c}} in the reduced set.
#'
#' @param object A Qboneobject
#' @param ... Arguments passed to other methods
#'
#' @importFrom ggplot2 ggplot geom_point geom_line
#'
#' @export
#'
nllplot <- function(
  object,
  ...
){
  message("Getting data for plotting.")
  plotdata1 <- loccplotdata(object, ...)
  plotdata2 <- plotdata1[[3]]
  p2 <- ggplot(plotdata2, aes(x=x, y=y, color= group)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(0.8, 1.0)) +
    scale_x_continuous(limits = c(-1, max(plotdata2$x)+0.5)) +
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
                   label=min(x[c(lasso.Chary_i_ - lasso.Chary1_i_) > 0.001])
    )
    ) +
    geom_label(aes(unique(plotdata2$x)[2],
                   cutoff,
                   label="cutoff")
    ) +
    annotate(geom = "text", x = c(0,unique(plotdata2$x)), y = 0.84, label = c("C",unique(plotdata2$x)), size = 3) +
    annotate(geom = "text", x = c(0,unique(plotdata2$x)), y = 0.83, label = c(expression(K[C]),unique(plotdata2$x2)), size = 3, angle = 45) +
    coord_cartesian(ylim = c(0.85, 1.01), xlim = c(0, max(plotdata2$x)+0.5), expand = FALSE, clip = "off") +
    theme_bw() +
    theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  return(suppressWarnings(print(p2)))
}

## 2.2 nllplotrev ----
#' Plot number of basis coefficients.
#'
#' Reversed nllplot useful when plotting with other basis like principal
#' components.
#' minimum concordance (\eqn{\rho^{o}}) and average (\eqn{\bar{\rho}})
#' for quantlets basis, varying with the number of basis coefficients.
#'
#' @param object A Qboneobject
#'
#' @importFrom ggplot2 ggplot geom_point
#'
#' @export
#'
nllplotrev <- function(){}

## 2.3 qbasisplot ----
#' Plot first n quantlet basis functions
#'
#' Plot first n quantlet basis functions
#'
#' @param object A Qboneobject
#' @param n Number of first n basis functions to plot, default = 16.
#'
#' @importFrom ggplot2 ggplot geom_point
#'
#' @export
#'
qbasisplot <- function(){}

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

## 6.1 concordCIDX ----
#' computes concordance correlation index
#'
#' @param yhat fitted values
#' @param y empirical true values
#'
#' @return concordance value
#'
#' @keywords internal
#'
#' @noRd
#'
concordCIDX <- function(yhat, y){
  2 * cov(yhat, y, use = "complete.obs") / (cov(yhat, yhat, use = "complete.obs") + cov(y, y, use = "complete.obs") + (mean(y, na.rm = TRUE) - mean(yhat, na.rm = TRUE))^2)
}

## 6.2 concordCIMX ----
#' computes concordance correlation index as in matrix
#'
#' @param yhat fitted values
#' @param y empirical true values
#'
#' @return concordance value for the matrix
#'
#' @keywords internal
#'
#' @noRd
#'
concordCIMX <- function(input.yhat, input.y){
  k <- dim(input.yhat)[2]
  values <- rep(NA, k)
  # pb <- txtProgressBar(min = 0, max = length(k), style = 3)
  for (i in 1:k) {
    values[i] <- concordCIDX(input.y[, i], input.yhat[, i])
  }
  # close(pb)
  return(values)
}

## 6.3 loccplotdata ----
#' Get data for plotting locc result
#'
#' @param object Qbone object
#' @param p Vector of length P  in (0,1)	Probability grids.
#' @param cutoff Near-lossless values#'
#' @param ... Arguments passed to other methods
#'
#' @return concordance value for the matrix
#'
#' @keywords internal
#'
#' @noRd
#'
loccplotdata <- function(
  object,
  p = signif(seq(0.001, 0.999, length = 1024), 4),
  cutoff = 0.990,
  ...
  ){
  plength = length(p)
  locc <- object@assays[[defaultAssay(object)]]@scale.data[["locc"]][[1]]
  # checking defaultAssay(object):
  raw.dataset <- getQboneData(object, slot = 'data',
                              assay = object@assays[[
                                object@assays[[defaultAssay(object)]]@assay.orig
                                ]]@assay.orig
                              )
  basis.columns <- object@assays[[defaultAssay(object)]]@scale.data[["basis.columns"]]
  remain.basis <- object@assays[[defaultAssay(object)]]@scale.data[["remain.basis"]]
  remain.counts <- object@assays[[defaultAssay(object)]]@scale.data[["remain.counts"]]
  n <- length(raw.dataset)
  y.mat <- matrix(NA, nrow = max(unlist(lapply(raw.dataset, length))), ncol = length(raw.dataset))
  for (i in 1:n) {
    y.mat[(1:length(raw.dataset[[i]])), i] <- raw.dataset[[i]]
  }


  lasso.long <- dim(locc)[3]
  lasso.locc <- matrix(NA, nrow = n, ncol = lasso.long)

  for (j in 1:lasso.long) {
    # message("\n Computing concordance correlation index for sample ", object@assays[[defaultAssay(object)]]@meta.assays[["id"]][[2]], " for plotting.")
    lasso.locc[, j] <- concordCIMX(y.mat, locc[, , j])
  }

  lasso.Chary_i_ <- apply(lasso.locc, 2, mean)
  lasso.Chary1_i_ <- apply(lasso.locc, 2, min)

  lasso.x <- remain.basis[ (remain.basis < plength) ][1:lasso.long]
  lasso.x1 <- remain.counts[ (remain.basis < plength) ][1:lasso.long]

  id2 <- (lasso.x) <= (length(raw.dataset) - 1) ## length(lasso.x)
  id3 <- lasso.Chary1_i_ >= cutoff
  id1 <- id3 & id2
  if (sum(id1) == 0) {
    id <- id2
    this <- min(seq(length(id))[id])
  }
  if (sum(id1) != 0) {
    id <- id1
    this <- max(seq(length(id))[id])
  }

  quantlet.set <- basis.columns[[  this   ]] [-1]

  numbasis <- sort(unique(c(lasso.x)), decreasing = TRUE)
  plotxaxis <- seq(length(numbasis))
  plotdata = data.frame(x = c(plotxaxis, plotxaxis),
                        x2 = c(numbasis, numbasis),
                        y = c(lasso.Chary_i_, lasso.Chary1_i_),
                        group = c(rep("mean", length(lasso.Chary_i_)), rep("min", length(lasso.Chary1_i_)))
                        )

  outputs <- list(quantlet.set, lasso.Chary1_i_[ this ], plotdata)
  names(outputs) <- list("quantlet", "lccc", "plotdata")
  return(outputs)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
