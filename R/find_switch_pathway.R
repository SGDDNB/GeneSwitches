#' @title Find significantly changed pathways and switching timepoint
#'
#' @description Apply hypergeometric test to determine significantly changed
#'     pathways and determine the switching timepoint for each pathway
#'
#' @param N expressed genes in pathways
#' @param pathways a list of pathways with genes
#' @param toplotgl_ptw swiching genes to plot
#' @param sig_FDR FDR cut off for significant pathways
#' @param direction switching direction, up or down
#' @return A dataframe of significantly changed pathways with switching timepoints
#'
phyper_pathway <- function(N, pathways, toplotgl_ptw, sig_FDR = 0.05, direction = c("up", "down")){
  pv <- c()
  swt_pw <- c()
  for (i in 1:length(pathways)) {
    gp <- pathways[[i]]
    genestoplot <- intersect(rownames(toplotgl_ptw), gp)
    m <- length(gp) ## Number of "marked" elements, i.e. genes associated to this biological process
    n <- N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
    x <- length(genestoplot) ## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process
    k <- nrow(toplotgl_ptw) ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
    pv <- c(pv, phyper(q=x-1, m=m, n=n, k=k, lower.tail = FALSE))

    toplotgl <- toplotgl_ptw[genestoplot,]
    swt_pw <- c(swt_pw, median(toplotgl$switch_at_time))
  }
  switch_pw <- data.frame(switch_at_time = swt_pw, pvalue = pv)
  rownames(switch_pw) <- names(pathways)
  switch_pw$FDR <- p.adjust(switch_pw$pvalue, method = "BH")

  switch_pw_sig <- switch_pw[switch_pw$FDR < sig_FDR,]
  if (nrow(switch_pw_sig) > 0) {
    switch_pw_sig$direction <- direction
  }
  return(switch_pw_sig)
}

#' @title Merge redundant pathways
#'
#' @description Merge significant pathways that are with same genes over certain ratio
#'
#' @param switch_pw a data frame with significantly changed pathways
#' @param pathways a list of pathways with genes
#' @param ratio cutoff ratio for merging redundant pathways
#' @return A dataframe of merged pathways
#'
merge_pathways <- function(switch_pw, pathways, ratio){
  switch_pw <- switch_pw[order(switch_pw$FDR),]
  switch_pw$keep <- 1
  for (i in 1:(nrow(switch_pw)-1)) {
    if (switch_pw[i,]$keep == 1) {
      for (j in (i+1):nrow(switch_pw)) {
        pw1 <- pathways[[rownames(switch_pw)[i]]]
        pw2 <- pathways[[rownames(switch_pw)[j]]]
        overlapping <- intersect(pw1, pw2)
        if (length(overlapping)/length(pw2) > ratio) {
          switch_pw[j,]$keep <- 0
        }
      }
    }
  }
  switch_pw <- switch_pw[switch_pw$keep == 1,]
  return(subset(switch_pw, select = -c(keep)))
}

#' @title Reduce redundant pathways
#'
#' @description Reduce significant pathways that are with same genes over certain ratio,
#'     do up- and down- regulated pathways separately
#'
#' @param switch_pw a data frame with significantly changed pathways
#' @param pathways a list of pathways with genes
#' @param redundant_pw_rate cutoff ratio for merging redundant pathways
#' @export
#' @return A dataframe of reduced pathways
#'
reduce_pathways <- function(switch_pw, pathways, redundant_pw_rate = 0.8){
  switch_up <- switch_pw[switch_pw$direction == "up",]
  if (nrow(switch_up) > 1) {
    switch_up <- merge_pathways(switch_up, pathways, ratio = redundant_pw_rate)
  }
  switch_down <- switch_pw[switch_pw$direction == "down",]
  if (nrow(switch_down) > 1) {
    switch_down <- merge_pathways(switch_down, pathways, ratio = redundant_pw_rate)
  }
  # combine
  switch_pw_reduce <- rbind(switch_up, switch_down)
  switch_pw_reduce <- switch_pw_reduce[order(switch_pw_reduce$FDR),]
  return(switch_pw_reduce)
}


#' @title Find significantly changed pathways and switching timepoint
#'
#' @description This function finds significantly changed pathways and determine
#'     the switching timepoint for each pathway
#'
#' @param sce SingleCellExperiment
#' @param pathways a list of pathways with genes
#' @param toplotgl_ptw swiching genes to plot
#' @param sig_FDR FDR cut off for significant pathways
#' @return A dataframe of significantly changed pathways ordered by FDR
#'
#' @export
#'
find_switch_pathway <- function(scerowdata, pathways = msigdb_h_c2_c5, toplotgl_sig,
                                sig_FDR = 0.05) {
  # make all the genes from pathways into a vector
  gps <- c()
  for (i in 1:length(pathways)) {
    gps <- c(gps, pathways[[i]])
  }
  genesINpathways <- unique(gps)
  N <- length(intersect(rownames(scerowdata), genesINpathways))

  # up-regulated pathways
  toplotgl_ptw <- toplotgl_sig[rownames(toplotgl_sig) %in% genesINpathways & toplotgl_sig$direction == "up",]
  switch_pw_sigup <- phyper_pathway(N, pathways, toplotgl_ptw, sig_FDR = 0.05, direction = "up")
  # down-regulated pathways
  toplotgl_ptw <- toplotgl_sig[rownames(toplotgl_sig) %in% genesINpathways & toplotgl_sig$direction == "down",]
  switch_pw_sigdown <- phyper_pathway(N, pathways, toplotgl_ptw, sig_FDR = 0.05, direction = "down")
  # combine
  switch_pw_sig <- rbind(switch_pw_sigup, switch_pw_sigdown)
  if (nrow(switch_pw_sig) > 0) {
    switch_pw_sig$feature_type <- "pathways"
    switch_pw_sig$feature_name <- rownames(switch_pw_sig)
  }
  switch_pw_ord <- switch_pw_sig[order(switch_pw_sig$FDR),]

  return(switch_pw_ord)
}

#' @title Pathways ridge plot
#'
#' @description This function generates pathways ridge plots
#'
#' @param switch_pw_re significant pathways to plot
#' @param toplotgl_sig switching genes
#' @param direction switching direction of the pathway, up or down
#' @param orderbytime order the pathways by switching time (mean time of switching genes) if TRUE,
#'                    order the pathways by FDR if FALSE
#' @return A ggplot object representing the pathway ridge plot
#'
#' @importFrom ggridges geom_density_ridges
#' @export
#'
plot_pathway_density <- function(switch_pw_re, toplotgl_sig, pw_direction = c("up", "down"), orderbytime = TRUE){
  switch_pw_re <- switch_pw_re[switch_pw_re$direction %in% pw_direction,]
  if (orderbytime == TRUE) {
    switch_pw_re <- switch_pw_re[order(switch_pw_re$switch_at_time, decreasing = TRUE),]
  } else {
    switch_pw_re <- switch_pw_re[order(switch_pw_re$FDR, decreasing = TRUE),]
  }
  gl <- c()
  for (i in 1:nrow(switch_pw_re)) {
    pn <- rownames(switch_pw_re)[i]
    pgl <- data.frame(msigdb_h_c2_c5[pn], pn, stringsAsFactors = FALSE)
    colnames(pgl) <- c("Genes", "Pathways")
    tgn <- length(unique(pgl$Genes))
    genestoplot <- intersect(rownames(toplotgl_sig), pgl$Genes)
    pgl <- pgl[pgl$Genes %in% genestoplot,]
    toplotgl <- toplotgl_sig[pgl$Genes, ]
    pgl <- as.data.frame(cbind(pgl, toplotgl))
    pgl <- pgl[pgl$direction == switch_pw_re[i,]$direction,]
    pgl$Pathways <- paste0(pgl$Pathways, "(",nrow(pgl),"/",tgn,")")
    gl <- rbind(gl, pgl)
  }

  gl$Pathways <- factor(gl$Pathways, levels = unique(gl$Pathways))
  p <- ggplot(gl, aes(x = switch_at_time, y = Pathways, fill = direction, col = direction))
  p <- p + theme_classic()
  p <- p + xlab("Pseudo-timeline") + geom_density_ridges(alpha = 0.5) +
    labs(fill = "Regulation", color = "Regulation")
  p <- p + theme(text = element_text(size = 12, family = "Helvetica"),
                 panel.background = element_rect(fill = "white", colour = NA),
                 axis.line = element_line(colour = "black"),
                 legend.key.size = unit(10, "pt"),
                 legend.text = element_text(size = 10,colour = "black"),
                 legend.title = element_text(size = 11, colour = "black")) +
    scale_color_manual(values=c("forestgreen", "chocolate2")) +
    scale_fill_manual(values=c("forestgreen", "chocolate2"))
  return(p)
}
