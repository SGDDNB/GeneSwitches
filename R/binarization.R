
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to GeneSwitches!")
}

#' @title Binarize gene expression
#'
#' @description This function generates on/off binarized data for gene expression
#'
#' @param sce SingleCellExperiment
#' @param fix_cutoff Logical. if use fixed global cutoff for binarization, default FALSE
#' @param binarize_cutoff fixed global cutoff for binarization, default 0.2
#' @param ncores number of cores
#' @return
#'
#' @import parallel
#' @importFrom mixtools normalmixEM
#' @export
#'
binarize_exp <- function(sce, fix_cutoff = FALSE, binarize_cutoff = 0.2, ncores = 3) {
  # calculate zero percentage
  #zerop_g <- c()
  expdata <- assays(sce)$expdata
  #for (i in 1:nrow(expdata)) {
  #  zp <- length(which(expdata[i, ] == 0))/ncol(expdata)
  #  zerop_g <- c(zerop_g, zp)
  #}
  zerop_g <- apply(expdata==0, 1, sum)
  zerop_g <- zerop_g/ncol(expdata)

  if (fix_cutoff == TRUE) {
    expdata <- assays(sce)$expdata
    is.na(expdata) <- assays(sce)$expdata == 0
    exp_reduced_binary <- as.matrix((expdata > binarize_cutoff) + 0)
    exp_reduced_binary[is.na(exp_reduced_binary)] = 0
    assays(sce)$binary <- exp_reduced_binary
    oupBinary <- data.frame(geneID = rownames(sce),
                            zerop_gene = zerop_g,
                            passBinary = TRUE)
    rowData(sce) <- oupBinary
  } else {
    expdata <- assays(sce)$expdata
    # Add gaussian noise to gene expression matrix
    # Here we use a sd of 0.1
    LogCountsadd = expdata + matrix(rnorm(nrow(expdata)*ncol(expdata),
                                          mean = 0, sd = 0.1),
                                    nrow(expdata), ncol(expdata))
    # Start fitting mixture models for each gene
    oupBinary = do.call(
      rbind, mclapply(rownames(LogCountsadd), function(iGene){
        set.seed(42)   # Set seed for consistency
        tmpMix = normalmixEM(LogCountsadd[iGene, ], k = 2)
        if (tmpMix$mu[1] < tmpMix$mu[2]) {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[1],
                              mu2 = tmpMix$mu[2],
                              sigma1 = tmpMix$sigma[1],
                              sigma2 = tmpMix$sigma[2],
                              lambda1 = tmpMix$lambda[1],
                              lambda2 = tmpMix$lambda[2],
                              loglik = tmpMix$loglik)
        } else {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[2],
                              mu2 = tmpMix$mu[1],
                              sigma1 = tmpMix$sigma[2],
                              sigma2 = tmpMix$sigma[1],
                              lambda1 = tmpMix$lambda[2],
                              lambda2 = tmpMix$lambda[1],
                              loglik = tmpMix$loglik)
        }
        return(tmpOup)
      }, mc.cores = ncores))

    # Check if non-bimodal genes
    oupBinary$passBinary = TRUE
    oupBinary[oupBinary$lambda1 < 0.1, ]$passBinary = FALSE
    oupBinary[oupBinary$lambda2 < 0.1, ]$passBinary = FALSE
    oupBinary[(oupBinary$mu2 - oupBinary$mu1) < (oupBinary$sigma1 + oupBinary$sigma2), ]$passBinary = FALSE
    # table(oupBinary$passBinary)

    # Solve for intersection for remaining genes
    oupBinary$root = -1
    for(iGene in oupBinary[oupBinary$passBinary == TRUE, ]$geneID){
      tmpMix = oupBinary[oupBinary$geneID == iGene, ]
      tmpInt = uniroot(function(x, l1, l2, mu1, mu2, sd1, sd2) {
        dnorm(x, m = mu1, sd = sd1) * l1 -
          dnorm(x, m = mu2, sd = sd2) * l2},
        interval = c(tmpMix$mu1,tmpMix$mu2),
        l1 = tmpMix$lambda1, mu1 = tmpMix$mu1, sd1 = tmpMix$sigma1,
        l2 = tmpMix$lambda2, mu2 = tmpMix$mu2, sd2 = tmpMix$sigma2)
      oupBinary[oupBinary$geneID == iGene, ]$root = tmpInt$root
    }
    # Binarize expression
    binLogCounts = expdata[oupBinary$geneID,]
    binLogCounts = t(scale(t(binLogCounts), scale = FALSE,
                           center = oupBinary$root))
    binLogCounts[binLogCounts >= 0] = 1
    binLogCounts[binLogCounts < 0] = 0
    assays(sce)$binary <- binLogCounts

    oupBinary$zerop_gene <- zerop_g
    rowData(sce) <- oupBinary
  }
  return(sce)
}
