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
plot_gene_exp <- function(sce, geneofi, reduction,
                          downsample = FALSE, ds_cutoff = 0.7,
                          zero_ratio = 0.7, ptsize = 0.7, fitting = FALSE){
  i <- which(rownames(rowData(sce)) == geneofi)
  glmdata <- cbind(State = as.numeric(assays(sce)$binary[i, ]),
                   expvalue = as.numeric(assays(sce)$expdata[i, ]),
                   timedata = sce$Pseudotime,
                   Dim1 = reducedDim(sce, reduction)[,1],
                   Dim2 = reducedDim(sce, reduction)[,2])
  glmdata <- as.data.frame(glmdata)

  if (downsample == TRUE & round(sum(glmdata$State == 0)/nrow(glmdata),3) > ds_cutoff) {
    glmdata <- downsample_zeros(glmdata, ratio_ds = zero_ratio)
    geneofi <- paste0(geneofi, " (downsampled--",zero_ratio*100,"%)")
  }

  p1 <- ggplot(glmdata, aes(Dim1, Dim2, color = expvalue)) +
    ggtitle(geneofi) +
    geom_point(size = ptsize, alpha = 1.0) +
    scale_shape_manual(guide = FALSE, values = 16) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          text = element_text(size = 12.5, family = "Helvetica", face = "bold"),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(colour = "black", size = 1),
          legend.key = element_rect(fill = NA),
          legend.key.width = unit(10, "pt"),
          legend.position = "right") +
    scale_color_gradientn("",colours = c("grey85", brewer.pal(9, "OrRd")))
  ##State
  # expdata$State <- as.numeric(as.character(expdata$State))
  if (fitting == TRUE) {
    p2 <- ggplot(glmdata, aes(x=timedata, y=State)) + ggtitle(geneofi) +
      geom_point() +geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE)
    return(grid.arrange(p1, p2, nrow = 1))
  } else {
    plot(p1)
  }
}
