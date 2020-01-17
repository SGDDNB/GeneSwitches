#' @title Plot gene expression
#'
#' @description This function plot gene expression on two-dimensional space
#'
#' @param sce SingleCellExperiment
#' @param geneofi one gene of interest
#' @param reduction dimensional reduction method
#' @param downsample if do random downsampling of zeros
#' @param ds_cutoff only do downsampling if zero percentage is over this cutoff
#' @param zero_ratio downsampling zeros to this proportion
#' @param ptsize point size
#' @param fitting if plot logistic regression fitting
#' @return
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gridExtra grid.arrange
#' @export
#'
plot_gene_exp <- function(sce, gene, reduction, fitting = FALSE, bin_width = 0.1,
                          downsample = FALSE, ds_cutoff = 0.7, zero_ratio = 0.7, ptsize = 0.7){
  glmdata = data.frame(expvalue = assays(sce)$expdata[gene, ],
                      timedata = sce$Pseudotime,
                      State = as.numeric(assays(sce)$binary[gene, ]),
                      Dim1 = reducedDim(sce, reduction)[,1],
                      Dim2 = reducedDim(sce, reduction)[,2])

  if (downsample == TRUE & round(sum(glmdata$expvalue == 0)/nrow(glmdata),3) > ds_cutoff) {
    glmdata <- downsample_zeros(glmdata, ratio_ds = zero_ratio)
    geneofi <- paste0(geneofi, " (downsampled--",zero_ratio*100,"%)")
  }

  # Define plot theme
  plotTheme = theme(plot.title = element_text(size = 24, face = "bold"),
                    text = element_text(size = 18),
                    panel.background = element_rect(fill = "white", colour = NA),
                    plot.background  = element_rect(fill = "white", colour = NA),
                    axis.line = element_line(colour = "black", size = 1),
                    axis.text = element_text(size = 15, colour = "black"),
                    axis.ticks = element_line(colour = "black", size = 1),
                    legend.key = element_rect(fill = NA),
                    legend.key.width = unit(14, "pt"),
                    legend.position = "right")

  # Plot 1
  p1 <- ggplot(glmdata, aes(Dim1, Dim2, color = expvalue)) +
    geom_point(size = 0.7, alpha = 1.0) +
    scale_shape_manual(guide = FALSE, values = 16) +
    scale_color_gradientn("",colours = c("grey85", brewer.pal(9, "OrRd"))) +
    plotTheme + ggtitle(gene)

  # Plot3
  p3 <- ggplot(glmdata, aes(timedata, State)) + geom_point() +
    xlab("Pseudotime") + ylab("Probability") +
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
    geom_hline(yintercept = 0.5, linetype="dashed", color = "#E69F00", size=1) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    scale_x_continuous(breaks = seq(0, max(glmdata$timedata), 10)) +
    plotTheme + ggtitle(" ")

  if (fitting == TRUE & "root" %in% names(rowData(sce))) {
    # Plot 2
    ggGaussian = data.frame(x = seq(min(glmdata$expvalue) - 0.1,
                                    max(glmdata$expvalue) + 0.1, length.out = 500))
    ggGaussian$y1 = (nrow(glmdata)/10) * rowData(sce)[gene,]$lambda1 *
      dnorm(ggGaussian$x, mean = rowData(sce)[gene,]$mu1,
            sd = rowData(sce)[gene,]$sigma1)
    ggGaussian$y2 = (nrow(glmdata)/10) * rowData(sce)[gene,]$lambda2 *
      dnorm(ggGaussian$x, mean = rowData(sce)[gene,]$mu2,
            sd = rowData(sce)[gene,]$sigma2)
    p2 <- ggplot() + xlab("Log Expression") +
      geom_histogram(data = glmdata, aes(x=expvalue, fill = as.factor(State)),
                     binwidth = bin_width, alpha=0.7, position = "identity") +
      geom_path(data = ggGaussian, aes(x,y1), color = "#999999", size=.7) +
      geom_path(data = ggGaussian, aes(x,y2), color = "chocolate2", size=.7) +
      geom_vline(xintercept = rowData(sce)[gene,]$root, size=1, linetype="dashed") +
      scale_fill_manual("State", values = c("grey", "#E69F00"),
                        labels = c("Off", "On")) +
      plotTheme + theme(legend.position = "top")
    # Combine
    return(grid.arrange(p1, p2, p3, nrow = 1))
  } else if (fitting == TRUE) {
    return(grid.arrange(p1, p3, nrow = 1))
  } else {
    plot(p1)
  }
}
