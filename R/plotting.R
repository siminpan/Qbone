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

## 2.2 dxPlotRev ----   || double check issue with Qbone.test.prop0.0001.Rdata
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
    # if (v >= 7) {
    #   ylims <- c(-0.05, 0.055)
    # } else {
    #   ylims <- c(-0.2, 0.2)
    # }
    ylims <- c(min(object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]][, v]) - 0.05,
               max(object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]][, v]) + 0.05)
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

## 2.4 qbasisPlot3D ----
#' 3D Plot of first n quantlet basis functions
#'
#' Plot first n quantlet basis functions in 3D Plots
#'
#' @param object A Qboneobject
#' @param n Number of first n basis functions to plot, default = 16.
#' @param ... Arguments passed to other methods
#'
#' @importFrom plotly plot_ly layout
#' @importFrom magrittr %>%
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
qbasisPlot3D <- function(
  object,
  n = 16,
  ...
){
  df = data.frame(x = rep(object@assays[[defaultAssay(object)]]@scale.data[["p"]],dim(object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]])[2]+1),
                  y = c(rep(1, length(object@assays[[defaultAssay(object)]]@scale.data[["p"]])),
                        object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]]),
                  z = rep(1:(dim(object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]])[2]+1), each = dim(object@assays[[defaultAssay(object)]]@scale.data[["quantlets"]])[1])
  )
  axx <- list(
    title = "Percentile"
  )

  axy <- list(
    title = "Quantlets"
  )

  axz <- list(
    title = " "
  )
  p <- plot_ly(df[which(df$z <= n),], x = ~x, y = ~z, z = ~y, split = ~z,
               type = "scatter3d", mode = "lines", color = ~z,
               colors = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$z)))) %>%
    layout(title = "Quantlet &#936;",
           scene = list(xaxis=axx,yaxis=axy,zaxis=axz),
           legend = list(title=list(text='<b> Quantlets </b>')))
  return(p)
}

## 2.5 pdPlot ----
#' Predicted Density Plot
#'
#' Plot predicted densities
#'
#' @param object A Qboneobject
#' @param plot.col Columns to plot from \code{qfrModel()} results \code{object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["den_G"]] }. Corresponding to \code{X1} (new covariates) agrement in \code{qfrModel()}.
#' @param group.names Group name for plot legend for argument \code{plot.col}.
#' @param mean.diff T or F to add mean difference testing for two consecutive subjects posterior probability scores
#' @param var.diff T or F to add variance difference testing for two consecutive subjects posterior probability scores
#' @param skewed.diff T or F to add skewness difference testing for two consecutive subjects posterior probability scores
#' @param kurtosis.diff T or F to add kurtosis difference testing for two consecutive subjects posterior probability scores
#' @param ... Arguments passed to other methods
#'
#' @importFrom ggplot2 ggplot geom_line geom_smooth theme_bw ggtitle theme element_text element_blank annotate
#' @importFrom gridExtra grid.arrange
#'
#' @export
#'
pdPlot <- function(
  object,
  plot.col,
  group.names,
  mean.diff = F,
  var.diff = F,
  skewed.diff = F,
  kurtosis.diff = F,
  ...
){
  # x aixs
  data.x = data.frame(x = object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["xdomain"]][-1])
  # Density estimates
  data = as.data.frame(object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["den_G"]])[, c(plot.col)]
  if (hasArg(group.names)){
    colnames(data) = group.names
  }
  data.m = reshape2::melt(data)
  data.m = cbind(data.x, data.m)
  # Plot
  p <- ggplot(data.m,
              aes(x=x, y = value, color = variable)) +
    geom_line() +
    ggtitle("Predicted Densities") +
    # theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
    # geom_vline(xintercept = 0,
    #            linetype="dotted",
    #            color = "black",
    #            size=0.5)
  # Mean difference
  if (mean.diff){
    data.mu = as.data.frame(round(object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["mu_diff"]] ,3))
    for (i in 1:length(plot.col)){
      if (plot.col[i] %in% data.mu$V1){
        p <- p +
          annotate(x = -Inf, y = Inf, hjust = 0, vjust = 2 + (i - 1) * 1.5,
                   # x= quantile(data.m$x, probs = c(0.001)),
                   # y= quantile(data.m$value, probs = c(0.99)) + (i - 1) * quantile(data.m$value, probs = c(0.15)),
                   # hjust=0, vjust=0,
                   size = 3,
                   geom = "text",
                   label = paste0("Mean Shift ", colnames(data)[i], " (p=", data.mu[which(plot.col[i] == data.mu$V1),8], ")"))
      }
    }
  }
  # Variance difference
  if (var.diff){
    data.var = as.data.frame(round(object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["sigma_diff"]] ,3))
    for (i in 1:length(plot.col)){
      if (plot.col[i] %in% data.var$V1){
        p <- p +
          annotate(x = Inf, y = Inf, hjust = 1, vjust = 2 + (i - 1) * 1.5,
                   # x= quantile(data.m$x, probs = c(0.75)),
                   # y= quantile(data.m$value, probs = c(0.99)) + (i - 1) * quantile(data.m$value, probs = c(0.15)),
                   # hjust=0, vjust=0,
                   size = 3,
                   geom = "text",
                   label = paste0("Var Shift ", colnames(data)[i], " (p=", data.var[which(plot.col[i] == data.var$V1),8], ")"))
      }
    }
  }
  # Skewness difference
  if (skewed.diff){
    data.skewed = as.data.frame(round(object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["mu3_diff"]] ,3))
    for (i in 1:length(plot.col)){
      if (plot.col[i] %in% data.skewed$V1){
        p <- p +
          annotate(x = -Inf, y = -Inf, hjust = 0, vjust = -2 - (i - 1) * 1.5,
                   # x= quantile(data.m$x, probs = c(0.75)),
                   # y= quantile(data.m$value, probs = c(0.99)) + (i - 1) * quantile(data.m$value, probs = c(0.15)),
                   # hjust=0, vjust=0,
                   size = 3,
                   geom = "text",
                   label = paste0("Skewness Shift ", colnames(data)[i], " (p=", data.skewed[which(plot.col[i] == data.skewed$V1),8], ")"))
      }
    }
  }
  # Kurtosis difference
  if (kurtosis.diff){
    data.kurtosis = as.data.frame(round(object@assays[["Q.F.Regression"]]@scale.data[["mcmc_infer"]][["mu4_diff"]] ,3))
    for (i in 1:length(plot.col)){
      if (plot.col[i] %in% data.kurtosis$V1){
        p <- p +
          annotate(x = Inf, y = -Inf, hjust = 1, vjust = -2 - (i - 1) * 1.5,
                   # x= quantile(data.m$x, probs = c(0.75)),
                   # y= quantile(data.m$value, probs = c(0.99)) + (i - 1) * quantile(data.m$value, probs = c(0.15)),
                   # hjust=0, vjust=0,
                   size = 3,
                   geom = "text",
                   label = paste0("Kurtosis Shift ", colnames(data)[i], " (p=", data.kurtosis[which(plot.col[i] == data.kurtosis$V1),8], ")"))
      }
    }
  }
  return(suppressWarnings(print(p)))
}

## 2.6 quantileFPlot3D ----
#' 3D Plot of first n quantlet with observed and predicted quantile function
#'
#'
#'
#' @param object A Qboneobject
#' @param n Number of first n basis functions to plot, default = 16.
#' @param group Group to be plotted. All group will be plotted if NULL. Default is NULL.
#' @param data.assay It is the name of the assay whose data will be plotted
#' @param plot Plot "Observed", "Quantlets", "Predicted" or "All" which is "Observed -> Quantlets -> Predicted". Default is "All"
#' @param ... Arguments passed to other methods
#'
#' @importFrom plotly plot_ly layout subplot
#' @importFrom magrittr %>%
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
quantileFPlot3D <- function(
  object,
  n = 16,
  group = NULL,
  data.assay = defaultAssay(object),
  plot = "All",
  ...
){
  if (!(plot %in% c("Observed", "Quantlets", "Predicted", "All"))){
    stop('plot argument should be either "Observed", "Quantlets", "Predicted" or "All"')
  }
  # Get Data
  empCoefs.list <- getQboneData(object, slot = 'data', assay = data.assay)
  empCoefs <- matrix(unlist(empCoefs.list, use.names = F), nrow=length(empCoefs.list), byrow=TRUE)
  reduceBasis <- object@assays[[data.assay]]@scale.data[["reduceBasis"]]
  orthogBasis <- object@assays[[data.assay]]@scale.data[["orthogBasis"]]
  quantlets <- object@assays[[data.assay]]@scale.data[["quantlets"]]
  p <- object@assays[[data.assay]]@scale.data[["p"]]
  ## quantlets data
  df = data.frame(x = rep(p,dim(quantlets)[2]+1),
                  y = c(rep(1, length(p)),
                        quantlets),
                  z = rep(1:(dim(quantlets)[2]+1), each = dim(quantlets)[1])
  )
  ## Predicted
  df2 = eigenmapmmt(empCoefs, cbind(rep(1, length(p)),
                                    quantlets)
  )
  df2 = cbind(data.frame(id = object@meta.data[["id"]],
                         group = object@meta.data[["group"]]),
              df2)

  df3.m = reshape2::melt(df2)
  df3.m$variable = rep(unique(df$x), each = nrow(df2))
  df3.m$value = reRange(df3.m$value, range = c(min(df$y), max(df$y)))
  df3.m$group = paste("Predicted", df3.m$group)
  ## Observed
  observed = lapply(object@assays[["Thin"]]@data,
                    quantile,
                    prob = p)
  df4 = data.frame(variable = rep(p, length(observed)),
                   value = unlist(observed, use.names = F),
                   group = rep(object@meta.data[["group"]], each = length(p)),
                   id = rep(object@meta.data[["id"]], each = length(p))
  )

  df4$group = paste("Observed", df4$group)
  df4$value = reRange(df4$value, range = c(min(df$y), max(df$y)))
  if (is.null(group)){
    ## plot n basis
    mycolors0 = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$z)))
    fig0 <- plot_ly(df[which(df$z <= n),], x = ~x, y = ~z, z = ~y, split = ~z,
                    type = "scatter3d", mode = "lines", color= ~z,
                    colors = mycolors0)
    ## plot Predicted
    f <- list(
      family = "Courier New, monospace",
      size = 18,
      color = "black")
    a <- list(
      text = "Predicted",
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 1,
      showarrow = FALSE
    )
    mycolors1 = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df3.m$id)))
    fig1 <- plot_ly(df3.m,
                    x = ~variable, y = ~value, split = ~id,
                    type = "scatter", mode = "lines", color= ~id,
                    colors = mycolors1) %>%
      layout(annotations = a)
    ## plot Observed
    b <- list(
      text = "Observed",
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 1,
      showarrow = FALSE
    )
    mycolors2 = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df4$id)))
    fig2 <- plot_ly(df4,
                    x = ~variable, y = ~value, split = ~id,
                    type = "scatter", mode = "lines", color= ~id,
                    colors = mycolors2) %>%
      layout(annotations = b)
  } else {
    ## plot n basis
    mycolors0 = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$z)))
    fig0 <- plot_ly(df[which(df$z <= n),], x = ~x, y = ~z, z = ~y, split = ~z,
                    type = "scatter3d", mode = "lines", color= ~z,
                    colors = mycolors0)
    ## plot Predicted
    mycolors1 = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df3.m$id)))
    fig1 <- plot_ly(df3.m[grep(paste0(" ",unique(object@meta.data[["group"]])[group],"$"), df3.m$group), ],
                    x = ~variable, y = ~group, z = ~value, split = ~id,
                    type = "scatter3d", mode = "lines", color= ~id,
                    colors = mycolors1)
    ## plot Observed
    mycolors2 = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df4$id)))
    fig2 <- plot_ly(df4[grep(paste0(" ",unique(object@meta.data[["group"]])[group],"$"), df4$group), ],
                    x = ~variable, y = ~group, z = ~value, split = ~id,
                    type = "scatter3d", mode = "lines", color= ~id,
                    colors = mycolors2)
  }
  ## plot together
  axy <- list(
    title = "Observed -> Quantlets -> Predicted"
  )

  fig3 <- subplot(fig2, fig0, fig1) %>%
    layout(title = paste0('From Observed to Predicted \n Plot with First ', n, ' Basis'),
           scene = list(yaxis=axy)
    )
  # Output
  if (plot == "Observed"){
    return(suppressWarnings(print(fig2)))
  } else if (plot == "Quantlets"){
    return(suppressWarnings(print(fig0)))
  } else if (plot == "Predicted"){
    return(suppressWarnings(print(fig1)))
  } else if (plot == "All"){
    return(suppressWarnings(print(fig3)))
  }

}

## 2.7 histogram3D ----
#' Histogram Plot for each sample in 3D layout
#'
#'
#'
#' @param object A Qboneobject
#' @param title Title of the plot.
#' @param binbreaks bin width breaks for the Histogram. Defulat is 1000
#' @param data.assay It is the name of the assay whose data will be plotted
#' @param plotting  Plotting type, default is "scatter3d".
#' @param ... Arguments passed to other methods
#'
#' @importFrom plotly plot_ly layout subplot
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
histogram3D <- function(
  object,
  title = "Histogram",
  binbreaks = 1000,
  data.assay = defaultAssay(object),
  plotting = "scatter3d",
  ...
){
  # Get data
  df = data.frame(x = unlist(object@assays[[data.assay]]@data, use.names = F),
                  # y = unlist(mapply(seq, 1, unlist(lapply(object@assays[[data.assay]]@data, length), use.names = F)), use.names = F),
                  z = rep(names(object@assays[[data.assay]]@data), unlist(lapply(object@assays[[data.assay]]@data, length), use.names = F))
  )
  df$cut = cut(df$x, breaks = binbreaks)
  # Get histogram infor
  df1 <- df %>%
    group_by(z,cut) %>%
    summarise(Freq = n())
  df1$group = object@meta.data[["group"]][c(match(df1$z, object@meta.data[["id"]]))]
  # Plot
  mycolors = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df1$group)))
  if (plotting == "scatter3d"){
    fig <- plot_ly(df1, x = ~z, y = ~cut, z = ~Freq, split = ~group, type = plotting, mode = 'lines', colors = mycolors)
  } else if (plotting == "histogram"){
    fig <- plot_ly(df, x = ~x, split = ~z, type = plotting, colors = mycolors) %>%
      layout(barmode = "overlay")
  }
  fig1 <- fig  %>%
    layout(title = title)
  # Output
  return(suppressWarnings(print(fig1)))
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

## 6.4 reRange ----
#' rescale the data in diserd range
#'
#' @param x input data
#' @param range new data range
#'
#' @return rescaled value
#'
#' @keywords internal
#'
#' @noRd
#'
reRange <- function(x, range) {
  (x - min(x))/(max(x)-min(x)) * (range[2] - range[1]) + range[1]
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. TESTING ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
