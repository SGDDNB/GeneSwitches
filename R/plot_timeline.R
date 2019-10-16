#' @title Extract switching gene list of interesting
#'
#' @description This function extract a list of significant switching genes
#'
#' @param sce switching genes
#' @param allgenes if use all genes
#' @param pathway_name a list of pathway name(s) to plot
#' @param genelists a gene list to plot
#' @param genetype specific gene type to plot c("EMT", "reprogramming", "stem", "surface", "TF")
#' @param r2cutoff pseudo R^2 cutoff
#' @param direction switching direction, up or down
#' @param topnum number of top genes ordered by pseudo R^2 value
#' @return
#'
#' @import plyr
#' @export
#'
filter_switchgenes <- function(sce, allgenes = FALSE, pathway_name = NULL, genelists = GeneSwitches:::genelists,
                               genetype = c("Surface proteins", "TFs"),
                               r2cutoff = 0.03, direction = c("up", "down"), topnum = 100000) {
  if (allgenes == TRUE) {
    toplotgl <- rowData(sce)
    toplotgl$feature_type <- "All genes"
  } else if(!is.null(pathway_name)) {
    gl <- c()
    for (pn in pathway_name) {
      pgl <- data.frame(msigdb_h_c2_c5[pn], pn, stringsAsFactors = FALSE)
      colnames(pgl) <- c("feature_name", "feature_type")
      gl <- rbind(gl, pgl)
    }
    multi <- gl$feature_name[duplicated(gl$feature_name)]
    if (length(multi) > 0) {
      gl <- ddply(gl,.(feature_name),paste)[,c(1,3)]
      rownames(gl) <- gl$feature_name
      colnames(gl)[2] <- "types"
      gl$feature_type <- gl$types
      gl[multi,]$feature_type <- "Multiple"
      genestoplot <- intersect(rownames(sce), gl$feature_name)
      toplotgl <- rowData(sce)[genestoplot, ]
      toplotgl$feature_type <- gl[genestoplot, ]$feature_type
      toplotgl$types <- gl[genestoplot, ]$types
    } else {
      rownames(gl) <- gl$feature_name
      genestoplot <- intersect(rownames(sce), gl$feature_name)
      toplotgl <- rowData(sce)[genestoplot, ]
      toplotgl$feature_type <- gl[genestoplot, ]$feature_type
    }
  } else {
    genelists_sub <- genelists[genelists$genetypes %in% genetype, ]
    genelists_sub <- genelists_sub[!duplicated(genelists_sub$genenames), ]
    rownames(genelists_sub) <- genelists_sub$genenames
    genestoplot <- intersect(rownames(sce), genelists_sub$genenames)
    toplotgl <- rowData(sce)[genestoplot, ]
    toplotgl$feature_type <- genelists_sub[genestoplot, ]$genetypes
  }

  toplotgl_sub <- toplotgl[toplotgl$prd_quality == 1 & toplotgl$pseudoR2s > r2cutoff &
                           toplotgl$direction %in% direction, ]
  if (nrow(toplotgl_sub) > topnum) {
    toplotgl_sub <- toplotgl_sub[order(toplotgl_sub$pseudoR2s, decreasing = TRUE),]
    toplotgl_sub <- toplotgl_sub[1:topnum,]
  }
  return(toplotgl_sub)
}

#' @title Plot switching genes
#'
#' @description This function plots switching genes on the pseudo-timeline
#'
#' @param tml switching genes
#' @param timedata pseudotime for cells
#' @param iffulltml if plot the full timeline
#' @param txtsize text size for gene names
#' @param color_by the cell attribute (e.g. the column of tml) to map to each cell's color
#' @return
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'
plot_timeline_ggplot <- function(tml, timedata, iffulltml = TRUE, txtsize = 3.5, color_by = "feature_type") {
  tml <- as.data.frame(tml)
  tml <- tml[order(tml$switch_at_time), ]
  tml$direction_num <- -1
  if ("up" %in% tml$direction) {
    tml[tml$direction == "up", ]$direction_num <- 1
  }
  tml$color_by <- as.factor(tml[,color_by])
  tml$feature_name <- rownames(tml)
  head(tml)

  if (iffulltml) {
    pseudotime_step <- (max(timedata) - min(timedata))/4
    pseudotime_range <- seq(min(timedata), max(timedata), by = pseudotime_step)
    pseudotime_df <- data.frame(pseudotime_range, pseudotime_format = round(pseudotime_range, 1))
  } else {
    pseudotime_step <- (max(tml$switch_at_time) - min(tml$switch_at_time))/4
    pseudotime_range <- seq(min(tml$switch_at_time), max(tml$switch_at_time), by = pseudotime_step)
    pseudotime_df <- data.frame(pseudotime_range, pseudotime_format = round(pseudotime_range, 1))
  }

  tml_plot <- ggplot(tml, aes(x = switch_at_time, y = pseudoR2s * direction_num, col = color_by, label = feature_name)) +
    geom_point(size = txtsize/3) + xlab("Pseudo-timeline") + ylab("Quality of fitting (R^2)")
  tml_plot <- tml_plot + theme_classic()
  # Plot horizontal black line for timeline
  tml_plot <- tml_plot + geom_hline(yintercept = 0, color = "black", size = 0.6)
  tml_plot <- tml_plot + geom_label(data = pseudotime_df, aes(x = pseudotime_range, y = 0, label = pseudotime_format), size = (txtsize-0.5),
                                    color = "black")

  tml_plot <- tml_plot + geom_text_repel(aes(x = switch_at_time, y = pseudoR2s * direction_num, label = feature_name),
                                         size = txtsize, fontface = "bold", show.legend = FALSE)

  tml_plot <- tml_plot + theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(10, "pt"),
                               text = element_text(size = 12, family = "Helvetica"))

  return(tml_plot)
}
