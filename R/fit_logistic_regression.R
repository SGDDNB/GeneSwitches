#' @title Random downsampling of zero expression
#'
#' @description This function does randomly downsampling of cells with zero expression for one gene
#'
#' @param glmdata binary data of one gene
#' @param ratio_ds downsampling zeros to this proportion
#' @return
#'
downsample_zeros <- function(glmdata, ratio_ds = 0.7) {
  p = as.numeric(ratio_ds)
  set.seed(42)   # Set seed for consistency
  downsample <- sample(which(glmdata$expvalue == 0), length(which(glmdata$expvalue == 0)) - round(sum(glmdata$expvalue != 0) * p/(1 - p)))
  if (length(downsample) > 0) {
    subdata <- glmdata[-downsample, ]
  } else {subdata <- glmdata}
  return(subdata)
}

#' @title Fit fast logistic regression and find switching timepoint
#'
#' @description This function fits fast logistic regression and find switching timepoint for each gene
#'
#' @param sce SingleCellExperiment
#' @param downsample Logical. if do random downsampling of zeros
#' @param ds_cutoff only do downsampling if zero percentage is over this cutoff
#' @param zero_ratio downsampling zeros to this proportion
#' @param sig_FDR FDR cut off for significant genes
#' @return
#'
#' @import fastglm
#' @export
#'
find_switch_logistic_fastglm <- function(sce, downsample = FALSE, ds_cutoff = 0.7, zero_ratio = 0.7,
                                         sig_FDR = 0.05, show_warnings = TRUE) {
  binarydata <- assays(sce)$binary
  expdata <- assays(sce)$expdata
  binarydata <- binarydata[which(rowData(sce)$passBinary == TRUE), ]
  expdata <- expdata[which(rowData(sce)$passBinary == TRUE), ]
  genes <- rowData(sce)[which(rowData(sce)$passBinary == TRUE), ]
  timedata <- sce$Pseudotime
  pvalues <- binarydata[, 1]
  pseudoR2s <- binarydata[, 1]
  estimates <- binarydata[, 1]
  switch_at_time <- binarydata[, 1]
  prd_quality <- binarydata[, 1]
  CI <- binarydata[, 1]

  for (i in 1:nrow(binarydata)) {
    glmdata <- cbind(State = as.numeric(binarydata[i, ]), expvalue = as.numeric(expdata[i, ]),
                     timedata = sce$Pseudotime)
    glmdata <- as.data.frame(glmdata)

    if (downsample == TRUE & round(genes$zerop_gene[i],3) > ds_cutoff) {
      glmdata <- downsample_zeros(glmdata, ratio_ds = zero_ratio)
    }

    if (show_warnings == TRUE) {
      glm_results <- fastglm(x = model.matrix(State ~ timedata, data = glmdata),
                             y = glmdata$State, family = binomial(link = "logit"))
    } else {
      glm_results <-suppressWarnings(fastglm(x = model.matrix(State ~ timedata, data = glmdata),
                                             y = glmdata$State, family = binomial(link = "logit")))
    }
    pvalues[i] <- coef(summary(glm_results))[, 4][2]
    ll.null <- glm_results$null.deviance/-2
    ll.proposed <- glm_results$deviance/-2
    # McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
    pseudoR2s[i] <- (ll.null - ll.proposed)/ll.null
    estimates[i] <- coef(summary(glm_results))[, 1][2]
    # p=0.5
    switch_at_time[i] <- (log(0.5/(1 - 0.5)) - coef(glm_results)[1])/coef(glm_results)[2]
    if (switch_at_time[i] >= max(glmdata$timedata)) {
      switch_at_time[i] = max(glmdata$timedata)
      prd_quality[i] = 0
    } else {
      prd_quality[i] = 1
    }
    if (switch_at_time[i] <= min(glmdata$timedata)) {
      switch_at_time[i] = min(glmdata$timedata)
      prd_quality[i] = 0
    }
    se <- summary(glm_results)$coefficients[, 2]
    CI[i] <- sqrt((se[1]*1.96/coef(glm_results)[1])^2 + (se[2]*1.96/coef(glm_results)[2])^2)*
      abs(coef(glm_results)[1]/coef(glm_results)[2])
    remove(glm_results)
  }

  result_switch <- cbind(switch_at_time, CI, pvalues, pseudoR2s, estimates, prd_quality)
  rownames(result_switch) <- rownames(binarydata)
  result_switch <- as.data.frame(result_switch)
  result_switch$direction <- "up"
  result_switch[result_switch$estimates < 0, ]$direction <- "down"
  result_switch$FDR <- p.adjust(result_switch$pvalues, method = "BH")
  steptime <- (max(timedata) - min(timedata))/100
  result_switch$switch_at_timeidx <- round((result_switch$switch_at_time - min(timedata))/steptime)

  # process_resultswitch --------------------------------------------------------------------------
  # check significance FDR < sig_FDR
  if (max(result_switch$FDR) > sig_FDR) {
    result_switch[result_switch$FDR > sig_FDR, ]$prd_quality <- 0
  }

  geneinfo <- merge(rowData(sce), result_switch, by=0, all=TRUE)[,-1] #[,1:11]
  rownames(geneinfo) <- geneinfo$geneID
  # all(rownames(geneinfo) == rownames(sce))
  geneinfo <- geneinfo[rownames(sce), ]
  rowData(sce) <- geneinfo
  return(sce)
}

