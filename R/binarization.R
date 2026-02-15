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
#' @param gaussian_weight_cutoff float if the second gaussian has a weight (lambda) less than this cutoff assume the gene is not expressed.
#' @return A SingleCellExperiment object with an added binary assay in `assays(sce)$binary` and updated `rowData(sce)` containing binarization metadata.
#'
#' @import parallel
#' @importFrom pbapply pblapply
#' @importFrom mixtools normalmixEM
#' @importFrom Matrix rowSums colSums t Matrix
#' @export
#'
binarize_exp <- function(sce, fix_cutoff = FALSE, binarize_cutoff = 0.2, ncores = 3, gaussian_weight_cutoff = 0.03) {
  
  # Check if binarize_cutoff is valid -must be greater than 0 and numeric
  if (!is.numeric(binarize_cutoff) || binarize_cutoff < 0) {
    stop("The 'binarize_cutoff' must be a numeric value greater than or equal to 0.")
  }

  # Check if expdata exists
  if (is.null(assays(sce)$expdata)) {
    stop("Expression data 'expdata' not found in the SingleCellExperiment object assays.")
  }

  # load in the expression data from the sce object
  expdata <- assays(sce)$expdata

  # ensure expdata is a sparse matrix for downstream efficiency
  if (!inherits(expdata, "dgCMatrix")) {
    message("assays(sce)$expdata is not sparse (dgCMatrix).\nCoercing expression data to dgCMatrix for efficiency...")
    expdata <- as(expdata, "dgCMatrix")
  }

  # calculate the percentage of zeros for each gene
  # this variable is not used in the current binarization funciton but used downstream?
  # assumes expdata is a sparse matrix
  # uses Matrix package to calculate the percentage of zeros for each gene
  zerop_g <- 1 - (Matrix::rowSums(expdata != 0) / ncol(expdata))

  if (fix_cutoff == TRUE) {
    exp_reduced_binary <- (expdata > binarize_cutoff) + 0 # Remains sparse 'dgCMatrix'
    assays(sce)$binary <- exp_reduced_binary
    oupBinary <- data.frame(geneID = rownames(sce),
                            zerop_gene = zerop_g,
                            passBinary = TRUE) # this may be missleading as no checks were done - MTN 04/02/26
    rowData(sce) <- oupBinary
  } else {
    # transpose the expdata to have genes as columns and cells as rows for easier processing in the parallel loop
    expdata_t <- Matrix::t(expdata)

    # before fitting the mixture models, filter out genes that are expressed in less than 5 cells,
    # as these are likely to fail the bimodality tests and cause long runtimes due to non-convergence of the normalmixEM algorithm.
    # these will later be binarized to all 0's, which is likely more accurate than trying to fit a bimodal distribution to a gene that is only expressed in a few cells. 
    n_cells_expressed <- Matrix::colSums(expdata_t > 0)
    low5_exp_genes <- names(n_cells_expressed)[n_cells_expressed <= 5]
    if (length(low5_exp_genes) > 0) {
        message(paste(length(low5_exp_genes), "genes are expressed in 5 or fewer cells and will be binarized to 0's."))
    }

    # filter out these genes from expdata_t for the mixture model fitting
    expdata_t <- expdata_t[, !(colnames(expdata_t) %in% low5_exp_genes)]

    # --- MIXTURE MODEL FITTING ---
    message("Fitting mixture models for each gene... (use multiple cores for faster processing!)")
    # Start fitting mixture models for each gene
    # If on Windows, use 1 core (safe). If on Linux, use 'ncores' (fast).
    use_cores <- if(.Platform$OS.type == "windows") 1 else ncores
    oupBinary = do.call(
        rbind, pbapply::pblapply(colnames(expdata_t), function(iGene){ # itterate through each gene (cols of expdata_t)

        # START TRYCATCH to avoid errors from the normalmixEM function which can fail to converge for some genes.
        tryCatch({

        # Extract sparse vector for gene
        raw_counts <- as.numeric(expdata_t[, iGene])

        # Add gaussian noise to the gene expression to smooth the distribution for better fitting of mixture models.
        # # Here we use a sd of 0.1
        set.seed(42)   # Set seed for consistency
        noisy_counts <- raw_counts + rnorm(length(raw_counts), mean = 0, sd = 0.1)

        
        # fit mixture model with two components per gene
        # stricter max itterations and restarts to avoid long runtimes on genes that fail to converge.
        # similar logic for the higher epsilon value.
        tmpMix = mixtools::normalmixEM(noisy_counts,
                                            k = 2,
                                            maxit = 400,
                                            maxrestarts = 15,
                                            epsilon = 1e-6,
                                            verb = FALSE)


        if (tmpMix$mu[1] < tmpMix$mu[2]) {
            tmpOup = data.frame(geneID = iGene,
                                mu1 = tmpMix$mu[1],
                                mu2 = tmpMix$mu[2],
                                sigma1 = tmpMix$sigma[1],
                                sigma2 = tmpMix$sigma[2],
                                lambda1 = tmpMix$lambda[1],
                                lambda2 = tmpMix$lambda[2],
                                loglik = tmpMix$loglik,
                                passBinary = TRUE)
        } else { # correct the order of the two gausians. 
            tmpOup = data.frame(geneID = iGene,
                                mu1 = tmpMix$mu[2],
                                mu2 = tmpMix$mu[1],
                                sigma1 = tmpMix$sigma[2],
                                sigma2 = tmpMix$sigma[1],
                                lambda1 = tmpMix$lambda[2],
                                lambda2 = tmpMix$lambda[1],
                                loglik = tmpMix$loglik,
                                passBinary = TRUE)
        }
        return(tmpOup)
        }, error = function(e) {
                # ERROR HANDLER: If normalmixEM fails, return NA values for this gene
                message("Error fitting gene: ", iGene)
                # likley due to non-convergence
                # Return a row with NA stats and passBinary = FALSE
                data.frame(geneID = iGene,
                    mu1 = NA, mu2 = NA,
                    sigma1 = NA, sigma2 = NA,
                    lambda1 = NA, lambda2 = NA,
                    loglik = NA,
                    passBinary = FALSE)
            })
            # END TRYCATCH

    }, cl = use_cores))

    # report number of genes which normalmixEM failed to fit a model 
    if (sum(oupBinary$passBinary == FALSE) > 0) {
        message(paste(sum(oupBinary$passBinary == FALSE), "genes failed to fit a mixture model and will removed."))
    }


    # --- Biomodality Tests ---

    # Check for genes that do not pass the bimodality tests

    # check for overlap between the two gausians,
    # if the distance between the means is less than the sum of the SDs, mark as non-bimodal
    # assume these are continuous genes that do not have a clear switch point and are not good candidates for binarization.

    # using which(oupBinary$passBinary) to ensure we don't calculate on NAs
    pass_idx <- which(oupBinary$passBinary) # Only look at currently passing genes
    diff_mu <- oupBinary$mu2[pass_idx] - oupBinary$mu1[pass_idx]
    sum_sigma <- oupBinary$sigma1[pass_idx] + oupBinary$sigma2[pass_idx]
    fail_sep <- pass_idx[diff_mu < sum_sigma]
    oupBinary$passBinary[fail_sep] <- FALSE
    # table(oupBinary$passBinary)
    # print the number of genes that fail the separation test
    if (length(fail_sep) > 0) {
        message(paste(length(fail_sep), "genes failed the separation test and will be removed."))
    }

    # We must use 'which' to avoid crashing on the genes that failed above (they introduce NA values)
    # Check if one gaussian dominates the other, if one gaussian has less than x% weight mark as non-bimodal
    # changing this to a variable, default should be more like ~3% based on my testing. - MTN 05/02/26

    # report the number of genes that fail the gaussian weight cutoff for each gaussian.
    idx_l2_too_small <- which(oupBinary$passBinary & oupBinary$lambda2 < gaussian_weight_cutoff)
    if (length(idx_l2_too_small) > 0) {
        message(paste(length(idx_l2_too_small), "genes have lambda2 <", gaussian_weight_cutoff, "and will be binarized to 0's."))
    }
    # Note: the logic above relies on the testing of distance between the means being done first.

    # filter for l1 being too small removed as it would remove housekeeping genes,
    # I cant find many/any scenarios where a gene with a small l1 would be a problem for binarization

    # These genes should be marked as non-bimodal and set to 0 expression in the binarization
    # initialize a new column in oupBinary to store the roots,
    # doing this early alows me to add Inf values for the genes that fail the lambda2 bimodality test.
    oupBinary$root <- NA
    # Inf values mean these genes will be binarized to all 0's
    # only do this for genes where oupBinary$passBinary is TRUE to avoid adding Inf values to genes that already failed the separation test and have passBinary = FALSE.
    # Use the integer indices (idx_l2_too_small) directly; they already filtered for passBinary == TRUE
    oupBinary$root[idx_l2_too_small] <- Inf

    # --- CALCULATE Roots/intersections ---

    # Identify the intersection between the two gausians 
    # aka the switching point of a switching gene

    # Identify the indicies of the genes which passed the tests for bimodality
    pass_indices = which(oupBinary$passBinary == TRUE & is.na(oupBinary$root))

    # Itterate through each gene and find the root/intersection
    for(i in pass_indices){
        # select the values for the current gene
        mu1 = oupBinary$mu1[i]
        mu2 = oupBinary$mu2[i]
        sd1 = oupBinary$sigma1[i]
        sd2 = oupBinary$sigma2[i]
        l1 = oupBinary$lambda1[i]
        l2 = oupBinary$lambda2[i]

        # Added tryCatch here because in rare cases the intersection is not between the two means, which causes uniroot to fail.
        tryCatch({
            # calculate the root/intersection between the two gaussians using uniroot
            tmpInt = uniroot(function(x) {
                dnorm(x, mean = mu1, sd = sd1) * l1 -
                dnorm(x, mean = mu2, sd = sd2) * l2},
                interval = c(mu1, mu2))
            # store the root in the dataframe
            oupBinary$root[i] = tmpInt$root
            
        }, error = function(e) {
            # If uniroot fails, we mark the gene as failed
            oupBinary$passBinary[i] <<- FALSE # Use double arrow to update outside function
        })
        
    }


        # --- BINARIZATION ---

    # match the order of genes in expdata_t
    # (because oupBinary might be in a different order after the parallel loop)
    rownames(oupBinary) <- oupBinary$geneID
    oupBinary <- oupBinary[colnames(expdata_t), ]

    # remove the genes that failed the bimodality tests and have NA roots, as these will cause errors in the binarization step.
    # these genes will be removed.
    genes_to_remove <- oupBinary$geneID[is.na(oupBinary$root) & oupBinary$passBinary == FALSE]
      if (length(genes_to_remove) > 0) {
        warning(paste("A total of:", length(genes_to_remove), "genes will be removed from sce."))
        expdata_t <- expdata_t[, !(colnames(expdata_t) %in% genes_to_remove)]
        oupBinary <- oupBinary[!(oupBinary$geneID %in% genes_to_remove), ]
        sce <- sce[!(rownames(sce) %in% genes_to_remove), ]
      }


    # Perform Sparse Binarization efficiently
    # binarize based on the root value for each gene, 
    # if the expression is greater than the root, it is 1, otherwise 0.

    # Transpose back to Genes x Cells immediately.
    # We do this because we need to compare each row (gene) to a specific threshold.
    mat_genes <- Matrix::t(expdata_t)

    # Force to dgCMatrix to ensure @i and @x slots are accessible
    # (t() can sometimes return dgTMatrix or other types)
    if (!inherits(mat_genes, "dgCMatrix")) {
        mat_genes <- as(mat_genes, "dgCMatrix") 
    }
    
    # Extract thresholds aligned with the rows of mat_genes
    gene_roots <- oupBinary$root
    
    # We leverage the internal structure of dgCMatrix:
    # @x contains non-zero values, @i contains their 0-based row indices.
    # We only need to update the non-zero values. 
    # (Assumption: gene_roots >= 0. Therefore 0 > root is always FALSE, so zeros stay zeros).
    # aka they are already binarized to 0's, so we only need to update the non-zero values to 1's or 0's based on the root comparison.
    
    # Get 1-based row indices for every non-zero element
    nonZero_row_indices <- mat_genes@i + 1
    # Compare each non-zero value to its specific gene threshold
    # If val > root -> 1. If val <= root -> 0. 
    # Update the matrix values
    mat_genes@x  <- (mat_genes@x > gene_roots[nonZero_row_indices]) + 0

    # Drop explicit zeros (values that were <= root became 0) to keep matrix sparse
    binLogCounts <- Matrix::drop0(mat_genes)

    # add back the genes that were filtered out for being expressed in less than 5 cells, and binarize them to all 0's.
    if (length(low5_exp_genes) > 0) {
        # Create sparse matrix of 0s:
        # with the same number of columns as low5_exp_genes and the same number of rows as expdata_t.
        low5_binary <- Matrix::Matrix(0, nrow = length(low5_exp_genes), ncol = ncol(binLogCounts),
                                      dimnames = list(low5_exp_genes, colnames(binLogCounts)))
        
        # also add them to the metadata dataframe oupBinary
        # add the skipped genes back to oupBinary 
        # with passBinary = TURE (as they will be binarized to all 0's, and should be used for downstream analysis)
        # set the other parameters asside from geneID and passBinary to NA as they were not tested. 
        oupSkipped_genes <- data.frame(geneID = low5_exp_genes,
                                       mu1 = NA, mu2 = NA,
                                       sigma1 = NA, sigma2 = NA,
                                       lambda1 = NA, lambda2 = NA,
                                       loglik = NA,
                                       passBinary = TRUE,
                                       root = NA)
        rownames(oupSkipped_genes) <- oupSkipped_genes$geneID
        oupBinary <- rbind(oupBinary, oupSkipped_genes)
        
        # Bind this with the binarized matrix
        # (Genes as Rows as we have transposed back)
        binLogCounts <- rbind(binLogCounts, low5_binary)
    }

    # ensure the order of genes in binLogCounts matches the order of genes in sce
    binLogCounts <- binLogCounts[rownames(sce), ]
    # Store in SingleCellExperiment
    assays(sce)$binary <- binLogCounts

    # Add metdata 
    # ensure the order of genes in oupBinary matches the order of genes in sce
    oupBinary <- oupBinary[rownames(sce), ]
    # zerop_g to oupBinary
    oupBinary$zerop_gene <- zerop_g[oupBinary$geneID]
    # oupBinary to rowData of sce
    rowData(sce) <- oupBinary

    # print the number of genes that passed the binarization
    # best print this to inform the user about the stringency of the tests and potential need to adjust parameters.
    message(paste("Number of binarized genes:", sum(oupBinary$passBinary)))

    message("Done! Binary assay added to 'assays(sce)$binary'.")
    #add a warning if binary includes NA's # They should not include NA's as we have removed the genes that fail the bimodality tests, but this is just a safety check.
    if (any(is.na(assays(sce)$binary))) {
        warning("The binary assay contains NA values, which may cause issues in downstream analysis if not removed.")
    }
  }
  return(sce)
}
