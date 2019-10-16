#' @title plot monocle2 trajectory colored by State
#'
#' @description This function plots monocle2 trajectory with "State" colors
#'
#' @param monocle2_obj monocle2 output object
#' @import Biobase
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @return
#'
#' @export
#'
plot_monocle_State <- function(monocle2_obj){
  mcells <- pData(monocle2_obj)
  mcells$dim1 <- monocle2_obj@reducedDimS[1,]
  mcells$dim2 <- monocle2_obj@reducedDimS[2,]

  p1 <- ggplot(mcells, aes(dim1, dim2, color = State)) +
    ylab("Component 2") + xlab("Component 1") +
    geom_point(size = 1, alpha = 1.0) +
    scale_shape_manual(guide = FALSE, values = 16) +
    theme(text = element_text(size = 12, family = "Helvetica"),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = NA),
          legend.key.width = unit(12, "pt"),
          legend.text = element_text(size = 10,colour = "black"),
          legend.title = element_text(size = 11, colour = "black"),
          legend.position = "top")

  p2 <- ggplot(mcells, aes(dim1, dim2, color = Pseudotime)) +
    ylab("Component 2") + xlab("Component 1") +
    geom_point(size = 1.7, alpha = 1.0) +
    scale_shape_manual(guide = FALSE, values = 16) +
    theme(text = element_text(size = 12, family = "Helvetica"),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = NA),
          legend.key.width = unit(20, "pt"),
          legend.text = element_text(size = 10,colour = "black"),
          legend.title = element_text(size = 11, colour = "black"),
          legend.position = "top")
  return(grid.arrange(p1, p2, nrow = 1))
}

#' @title Convert monocle2 output into GeneSwitches object
#'
#' @description This function converts monocle2 output into GeneSwitches object
#'
#' @param monocle2_obj monocle2 output object
#' @param states a vector of states (path) that are interested in
#' @param logexpdata log-normal gene expression
#' @import Biobase
#' @return
#'
#' @export
#'
convert_monocle2 <- function(monocle2_obj, states, expdata){
  allcells <- pData(monocle2_obj)
  # extract cells and log-normal expression in certain path
  cells <- allcells[allcells$State %in% states,]
  expd <- expdata[,rownames(cells)]
  expd <- expd[apply(expd > 0,1,sum) >= 3,]
  # create GeneSwitches object
  sce <- SingleCellExperiment(assays = List(expdata = expd))
  # pass pseudotime info
  colData(sce)$Pseudotime <- cells$Pseudotime
  # pass reduced dims info
  rd <- t(cardiac_monocle2@reducedDimS)[rownames(cells),]
  colnames(rd) <- c("Component 1", "Component 2")
  reducedDims(sce) <- SimpleList(monocleRD=rd)

  return(sce)
}


#' @title Convert slingshot output into GeneSwitches object
#'
#' @description This function converts slingshot output into GeneSwitches object
#'
#' @param sce_slingshot slingshot SingleCellExperiment output object
#' @param pseudotime_idx name of desired pseudotime path to apply GeneSwitches
#' @return
#'
#' @export
#'
convert_slingshot <- function(sce_slingshot, pseudotime_idx, assayname = "expdata"){
  allcells <- as.data.frame(colData(sce_slingshot))
  # extract cells and log-normal expression in certain path
  cells <- allcells[!is.na(allcells[,pseudotime_idx]),]
  expd <- assays(sce_slingshot)[[assayname]][,rownames(cells)]
  expd <- expd[apply(expd > 0,1,sum) >= 3,]
  # create GeneSwitches object
  sce <- SingleCellExperiment(assays = List(expdata = expd))
  # pass pseudotime info
  colData(sce)$Pseudotime <- cells[,pseudotime_idx]
  # pass reduced dims info
  for (i in 1:length(reducedDims(sce_slingshot))) {
    reducedDims(sce)[[i]] <- reducedDims(sce_slingshot)[[i]][rownames(cells),]
  }
  names(reducedDims(sce)) <- names(reducedDims(sce_slingshot))

  return(sce)
}
