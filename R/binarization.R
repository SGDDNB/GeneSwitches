
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to GeneSwitches!")
}

#' @title Binarize gene expression
#'
#' @description This function generates on/off binarized data for gene expression
#'
#' @param sce SingleCellExperiment
#' @return
#'
#' @export
#'
binarize_exp <- function(sce, binarize_cutoff = 0.2) {
  exp_reduced <- assays(sce)$expdata
  is.na(exp_reduced) <- assays(sce)$expdata == 0
  exp_reduced_binary <- as.matrix((exp_reduced > binarize_cutoff) + 0)
  exp_reduced_binary[is.na(exp_reduced_binary)] = 0
  assays(sce)$binary <- exp_reduced_binary
  return(sce)
}
