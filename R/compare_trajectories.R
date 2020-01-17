#' @title Identify distinct switching genes for each path
#'
#' @description This function identifies distinct switching genes for each path
#'
#' @param toplotgl_Rsub1 switching genes of path1
#' @param toplotgl_Rsub2 switching genes of path2
#' @param path1name name of path1 given by user
#' @param path2name name of path2 given by user
#' @param r2cutoff pseudo R^2 cutoff
#' @return
#'
#' @export
#'
distinct_genes <- function(toplotgl_Rsub1, toplotgl_Rsub2, path1name = "Path1Genes", path2name = "Path2Genes",
                           r2cutoff = 0.05, scale_timeline = FALSE, path1time = NULL, path2time = NULL, bin = 100){
  toplotgl_Rsub1$genenames <- rownames(toplotgl_Rsub1)
  toplotgl_Rsub2$genenames <- rownames(toplotgl_Rsub2)

  gl1 <- toplotgl_Rsub1$genenames
  gl2 <- toplotgl_Rsub2$genenames
  glin1 <- setdiff(gl1, gl2)
  glin2 <- setdiff(gl2, gl1)
  gs_p1 <- toplotgl_Rsub1[glin1,]
  gs_p2 <- toplotgl_Rsub2[glin2,]
  gs_p1$Paths <- path1name
  gs_p2$Paths <- path2name
  if (scale_timeline == TRUE) {
    steptime1 <- (max(path1time) - min(path1time))/bin
    gs_p1$switch_at_time <- round((gs_p1$switch_at_time - min(path1time))/steptime1)
    steptime2 <- (max(path2time) - min(path2time))/bin
    gs_p2$switch_at_time <- round((gs_p2$switch_at_time - min(path2time))/steptime2)
  }
  # combine distinct gene into one dataframe
  toplotgl <- rbind(gs_p1[,c("geneID","zerop_gene","switch_at_time","pvalues","FDR","pseudoR2s",
                             "estimates","prd_quality","direction","switch_at_timeidx",
                             "genenames","Paths")],
                    gs_p2[,c("geneID","zerop_gene","switch_at_time","pvalues","FDR","pseudoR2s",
                             "estimates","prd_quality","direction","switch_at_timeidx",
                             "genenames","Paths")])
  toplotgl <- toplotgl[toplotgl$pseudoR2s > r2cutoff,]
  return(toplotgl)
}

#' @title Identify common switching genes between paths
#'
#' @description This function identifies common switching genes between two paths
#'
#' @param toplotgl_Rsub1 switching genes of path1
#' @param toplotgl_Rsub2 switching genes of path2
#' @param path1name name of path1 given by user
#' @param path2name name of path2 given by user
#' @param r2cutoff pseudo R^2 cutoff
#' @return
#'
#' @export
#'
common_genes <- function(toplotgl_Rsub1, toplotgl_Rsub2, path1name = "Path1Genes", path2name = "Path2Genes",
                         r2cutoff = 0.05){
  toplotgl_Rsub1 <- toplotgl_Rsub1[toplotgl_Rsub1$pseudoR2s > r2cutoff,]
  toplotgl_Rsub2 <- toplotgl_Rsub2[toplotgl_Rsub2$pseudoR2s > r2cutoff,]

  toplotgl_Rsub1$genenames <- rownames(toplotgl_Rsub1)
  toplotgl_Rsub2$genenames <- rownames(toplotgl_Rsub2)
  gl1 <- toplotgl_Rsub1$genenames
  gl2 <- toplotgl_Rsub2$genenames
  comgl <- intersect(gl1, gl2)

  if (all(toplotgl_Rsub1[comgl,]$direction != toplotgl_Rsub2[comgl,]$direction)){
    print("Directions are not consistent.")
    return(NULL)
  }

  toplotgl_Rsub1$genetype <- path1name
  toplotgl_Rsub2$genetype <- path2name
  ggData <- as.data.frame(rbind(toplotgl_Rsub1[comgl, c("geneID","zerop_gene","switch_at_time","pvalues","FDR","pseudoR2s",
                                                        "estimates","prd_quality","direction","switch_at_timeidx",
                                                        "genenames","genetype")],
                                toplotgl_Rsub2[comgl, c("geneID","zerop_gene","switch_at_time","pvalues","FDR","pseudoR2s",
                                                        "estimates","prd_quality","direction","switch_at_timeidx",
                                                        "genenames","genetype")]))
  ggData$genetype <- factor(ggData$genetype, levels = c(path1name, path2name))
  return(ggData)
}

#' @title Plot common switching genes between paths
#'
#' @description This function plots common switching genes between two paths
#'
#' @param ggData data frame for common genes
#' @param timedata timedata to show on plot
#' @return
#'
#' @export
#'
common_genes_plot <- function(ggData, timedata){
  ggData$genetypenum <- as.numeric(ggData$genetype)

  comtml_plot <- ggplot(ggData, aes(switch_at_time, genetypenum, group = genenames, col=direction)) +
    geom_line() + geom_point() + ylim(0.5, 3.5)# + xlim(20,28)

  comtml_plot<-comtml_plot+theme_classic()

  # Plot horizontal black line for timeline
  comtml_plot<-comtml_plot+geom_hline(yintercept=c(1,2),#as.integer(ggData$type),
                                      color = "black", size=0.6)

  pseudotime_step <- (max(timedata) - min(timedata))/10
  pseudotime_range <- seq(min(timedata), max(timedata), by=pseudotime_step)
  pseudotime_df <- data.frame(pseudotime_range, pseudotime_format=round(pseudotime_range,1))

  pathposy <- c(1+0.07, 2-0.07)
  pathposx <- max(pseudotime_range)-0.9
  pathnames <- levels(ggData$genetype)
  pathnames_df <- data.frame(pathposx, pathposy, pathnames)
  comtml_plot <- comtml_plot + geom_text(data = pathnames_df, inherit.aes = FALSE, size=3.5,
                                         aes(x=pathposx, y=pathposy, label=pathnames))

  ##text position
  positions <- seq(0.05, 1.15, by=0.1)
  tml <- as.data.frame(ggData[ggData$genetype == levels(ggData$genetype)[2],])
  tml <- tml[with(tml, order(switch_at_time)), ]
  tml_uni <- tml[!duplicated(tml$switch_at_time),]

  line_pos <- data.frame(
    "switch_at_time"=tml_uni$switch_at_time,
    "position"=rep(positions, length.out=length(tml_uni$switch_at_time))
  )

  tml <- merge(x=tml, y=line_pos, by="switch_at_time", all = TRUE)

  text_offset <- 0.05
  tml$switch_at_time_count <- ave(tml$switch_at_time==tml$switch_at_time, tml$switch_at_time, FUN=cumsum)
  tml$text_position <- (tml$switch_at_time_count * text_offset) + tml$position

  # Plot vertical segment lines for milestones
  comtml_plot<-comtml_plot+geom_segment(data=tml[tml$switch_at_time_count == 1,], aes(y=2+position,yend=2,xend=switch_at_time),
                                        color='black', size=0.4)#size=tml$pseudoR2s)

  # Plot scatter points at zero and date
  # tml_plot<-tml_plot+geom_point(aes(y=0), size=2)

  # Don't show axes, appropriately position legend
  comtml_plot<-comtml_plot+theme(axis.line.y=element_blank(),
                                 axis.text.y=element_blank(),
                                 axis.title.x=element_blank(),
                                 axis.title.y=element_blank(),
                                 axis.ticks.y=element_blank(),
                                 axis.text.x =element_blank(),
                                 axis.ticks.x =element_blank(),
                                 axis.line.x =element_blank(),
                                 legend.position = "bottom", legend.key.size = unit(12, "pt"),
                                 legend.text = element_text(size = 12),
                                 text = element_text(size = 12, family = "Helvetica")) +
    labs(fill = "Regulation", color = "Regulation") +
    scale_color_manual(values=c("forestgreen", "chocolate2"))

  comtml_plot<-comtml_plot+geom_text(data=pseudotime_df, inherit.aes = FALSE,
                                     aes(x=pseudotime_range, y=0.9, label=pseudotime_format), size=3.8, color='black')

  # Show text for each milestone
  comtml_plot<-comtml_plot+geom_text(data=tml, inherit.aes = FALSE,
                                     aes(x=switch_at_time,y=2+text_position,label=genenames),size=3)

  return(comtml_plot)
}

