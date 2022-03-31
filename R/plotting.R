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

## 2.1 dxPlot ----
#' Plot near-losslessness parameters.
#'
#' Plots the near-losslessness parameters minimum concordance (\eqn{\rho^{o}})
#' and average (\eqn{\bar{\rho}}) against the number of basis coefficients
#' function of \eqn{K_{c}} in the reduced set.
#'
#' @param object A Qboneobject
#' @param sparsity Sparsity regularization parameter.
#' @param ... Arguments passed to other methods
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_y_continuous scale_x_continuous scale_color_manual labs geom_vline geom_hline geom_label annotate coord_cartesian theme_bw theme element_blank unit
#'
#' @export
#'
dxPlot <- function(
  object,
  sparsity = 0.001,
  ...
){
  # message("Getting data for plotting.")
  plotdata1 <- loccplotdata(object, ...)
  plotdata2 <- plotdata1[[3]]
  plotdata3 <- plotdata2[which(plotdata2$y >= 0.85),]
  if (nrow(plotdata1[[3]]) != nrow(plotdata3)){
    message(nrow(plotdata1[[3]]) - nrow(plotdata3), " points removed from ", expression(rho^0), " as it's lower than 0.85")
  }
  lasso.mean_i_ = plotdata2$y[plotdata2$group == "mean"]
  lasso.min_i_ = plotdata2$y[plotdata2$group == "min"]
  p2 <- ggplot(plotdata2, aes(x=x, y=y, color= group)) +
    geom_point(data = plotdata3, aes(x=x, y=y, color= group)) +
    geom_line(data = plotdata3, aes(x=x, y=y, color= group)) +
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
    geom_vline(xintercept = min(plotdata2$x[c(lasso.mean_i_ - lasso.min_i_) > sparsity]),
               linetype="dotted",
               color = "black",
               size=0.5) +
    geom_hline(yintercept = plotdata2$cutoff,
               linetype="dotted",
               color = "black",
               size=0.5) +
    geom_label(aes(min(x[c(lasso.mean_i_ - lasso.min_i_) > sparsity]),
                   0.865,
                   label=paste0("C = ", min(x[c(lasso.mean_i_ - lasso.min_i_) > sparsity]), "\n",
                                expression(K[C]), " = ", max(x2[c(lasso.mean_i_ - lasso.min_i_) > sparsity]))
    )
    ) +
    geom_label(aes(unique(plotdata2$x)[2],
                   plotdata2$cutoff,
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

## 2.2 dxPlotRev ----
#' Plot number of basis coefficients.
#'
#' Reversed dxPlot useful when plotting with other basis like principal
#' components.
#' minimum concordance (\eqn{\rho^{o}}) and average (\eqn{\bar{\rho}})
#' for quantlets basis, varying with the number of basis coefficients.
#'
#' @param object A Qboneobject#'
#' @param sparsity Sparsity regularization parameter.
#' @param ... Arguments passed to other methods
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_y_continuous scale_x_reverse scale_color_manual labs geom_vline geom_hline geom_label annotate coord_cartesian theme_bw theme element_blank unit
#'
#' @export
#'
dxPlotRev <- function(
  object,
  sparsity = 0.001,
  ...
){
  # message("Getting data for plotting.")
  plotdata1 <- loccplotdata(object, ...)
  plotdata2 <- plotdata1[[3]]
  plotdata3 <- plotdata2[which(plotdata2$y >= 0.85),]
  if (nrow(plotdata1[[3]]) != nrow(plotdata3)){
    message(nrow(plotdata1[[3]]) - nrow(plotdata3), " points removed from ", expression(rho^0), " as it's lower than 0.85")
  }
  lasso.mean_i_ = plotdata2$y[plotdata2$group == "mean"]
  lasso.min_i_ = plotdata2$y[plotdata2$group == "min"]
  p3 <- ggplot(plotdata2, aes(x=x, y=y, color= group)) +
    geom_point(data = plotdata3, aes(x=x, y=y, color= group)) +
    geom_line(data = plotdata3, aes(x=x, y=y, color= group)) +
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
    geom_vline(xintercept = min(plotdata2$x[c(lasso.mean_i_ - lasso.min_i_) > sparsity]),
               linetype="dotted",
               color = "black",
               size=0.5) +
    geom_hline(yintercept = plotdata2$cutoff,
               linetype="dotted",
               color = "black",
               size=0.5) +
    geom_label(aes(min(x[c(lasso.mean_i_ - lasso.min_i_) > sparsity]),
                   0.858,
                   label=paste0(expression(K[C]), " = ", max(x2[c(lasso.mean_i_ - lasso.min_i_) > sparsity]))
    )
    ) +
    geom_label(aes(unique(plotdata2$x)[2],
                   plotdata2$cutoff,
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
  return(suppressWarnings(print(p3)))
}

## 2.3 qbasisPlot ----
#' Plot first n quantlet basis functions
#'
#' Plot first n quantlet basis functions
#'
#' @param object A Qboneobject
#' @param n Number of first n basis functions to plot, default = 16.
#' @param ... Arguments passed to other methods
#'
#' @importFrom ggplot2 ggplot geom_line geom_smooth theme_bw ggtitle theme element_text element_blank
#' @importFrom gridExtra grid.arrange
#'
#' @export
#'
qbasisPlot <- function(
  object,
  n = 16,
  ...
){
  # use base R
  # par(mfrow = c(4, 4), mar = c(4, 2, 3, 2))
  # plot(object@assays[[defaultAssay(object)]]@scale.data[["p"]],
  #      rep(1, length(object@assays[[defaultAssay(object)]]@scale.data[["p"]])),
  #      type = "l", lty = 1, lwd = 0.2, main = bquote("Quantlet" ~ psi[.(1)]))
  # for (v in 1:(n-1)) {
  #   if (v >= 7) {
  #     ylims <- c(-0.05, 0.05)
  #   }
  #   if (v < 7) {
  #     ylims <- c(-0.2, 0.2)
  #   }
  #   plot(object@assays[[defaultAssay(object)]]@scale.data[["p"]],
  #        object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]][, v],
  #        type = "l", lty = 1, lwd = 0.2, main = bquote("Quantlet" ~ psi[.(v + 1)]),
  #        ylim = ylims, xlab = "")
  # }
  # use ggplot2
  p0 <- ggplot(data.frame(x = object@assays[[defaultAssay(object)]]@scale.data[["p"]],
                          y = rep(1, length(object@assays[[defaultAssay(object)]]@scale.data[["p"]]))),
               aes(x=x, y = y)) +
    geom_line() +
    ggtitle(bquote("Quantlet" ~ psi[.(1)])) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
  list <- list(p0)
  for (v in 1:(n-1)) {
    if (v >= 7) {
      ylims <- c(-0.05, 0.055)
    } else {
      ylims <- c(-0.2, 0.2)
    }
    p1 <- ggplot(data.frame(x = object@assays[[defaultAssay(object)]]@scale.data[["p"]],
                            y = object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]][, v]),
                 aes(x=x, y = y)) +
      geom_line() +
      # geom_smooth(formula = y ~ s(x, bs = "cs"), se=F, method = 'gam', color = "black", size=0.5) +
      scale_y_continuous(limits = ylims) +
      ggtitle(bquote("Quantlet" ~ psi[.(v + 1)])) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5))
    list <- append(list, list(p1))
  }
  return(
    suppressWarnings(
      grid.arrange(grobs = list,
                   nrow = 4)
    )
  )
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
#' @param p Probability grids.
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
  p = object@assays[[defaultAssay(object)]]@scale.data[["p"]],
  cutoff = 0.990,
  ...
  ){
  plength = length(p)
  locc <- object@assays[[defaultAssay(object)]]@scale.data[["locc"]][[1]]
  # checking defaultAssay(object):
  raw.dataset <- getQboneData(object, slot = 'data',
                              assay = object@assays[[
                                "Lasso.list"
                                #object@assays[[defaultAssay(object)]]@assay.orig
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

  lasso.mean_i_ <- apply(lasso.locc, 2, mean)
  lasso.min_i_ <- apply(lasso.locc, 2, min)

  lasso.x <- remain.basis[ (remain.basis < plength) ][1:lasso.long]
  lasso.x1 <- remain.counts[ (remain.basis < plength) ][1:lasso.long]

  id2 <- (lasso.x) <= (length(raw.dataset) - 1) ## length(lasso.x)
  id3 <- lasso.min_i_ >= cutoff
  id1 <- id3 & id2
  if (sum(id1) == 0) {
    id <- id2
    this <- min(seq(length(id))[id])
  } else {
    id <- id1
    this <- max(seq(length(id))[id])
  }

  quantlet.set <- basis.columns[[this]][-1]

  numbasis <- sort(unique(c(lasso.x)), decreasing = TRUE)
  plotxaxis <- seq(length(numbasis))
  plotdata = data.frame(x = c(plotxaxis, plotxaxis),
                        x2 = c(numbasis, numbasis),
                        y = c(lasso.mean_i_, lasso.min_i_),
                        group = c(rep("mean", length(lasso.mean_i_)), rep("min", length(lasso.min_i_))),
                        cutoff = cutoff
                        )

  outputs <- list(quantlet.set, lasso.min_i_[ this ], plotdata)
  names(outputs) <- list("quantlet", "lccc", "plotdata")
  return(outputs)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
