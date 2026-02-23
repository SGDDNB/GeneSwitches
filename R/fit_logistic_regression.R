#' @title Fit fast logistic regression and find switching timepoint
#'
#' @description This function fits fast logistic regression and find switching timepoint for each gene
#'
#' @param sce SingleCellExperiment
#' @param downsample Logical. if do random downsampling of zeros
#' @param ds_cutoff only do downsampling if zero percentage is over this cutoff
#' @param zero_ratio downsampling zeros to this proportion
#' @param sig_FDR FDR cut off for significant genes
#' @return A SingleCellExperiment object with the regression results added to the rowData.
#'
#' @import fastglm
#' @export
#'
find_switch_logistic_fastglm <- function(sce,
                                         downsample = TRUE,
                                         ds_cutoff = 0.7,
                                         zero_ratio = ds_cutoff,
                                         sig_FDR = 0.05) {
  # CHECK: Ensure the necessary assays are present in the SingleCellExperiment object
  if (!("binary" %in% names(assays(sce))) || !("expdata" %in% names(assays(sce)))) {
    stop("The SingleCellExperiment object must contain 'binary' and 'expdata' assays. Please run 'binarize_exp()' first to create these assays.")
  }                                           

  # load the data binarized by binarize_exp(), and the original expression data
  binarydata <- assays(sce)$binary
  expdata <- assays(sce)$expdata 

  # this function assumes the matrices to be sparse
  # ensure expdata is a sparse matrix for downstream efficiency
  if (!inherits(expdata, "dgCMatrix")) {
    message("assays(sce)$expdata is not sparse (dgCMatrix).", 
            "Coercing expression data to dgCMatrix for efficiency...")
    expdata <- as(expdata, "dgCMatrix")
  }
  # ensure binarydata is a sparse matrix for downstream efficiency
  if (!inherits(binarydata, "dgCMatrix")) {
    message("assays(sce)$binary is not sparse (dgCMatrix).",
            "Coercing expression data to dgCMatrix for efficiency...")
    binarydata <- as(binarydata, "dgCMatrix")
  }

  # --- TRANSPOSE ---
  # looping through a sparse matrix's columns is faster than looping through rows.
  binary_t <- Matrix::t(binarydata)
  exp_t <- Matrix::t(expdata)

  # --- Load/PRE-CALCULATE CONSTANTS ---
  # initialize vectors to store the results of the regression for each gene
  # these vectors have the same length as the number of genes 
  # intialize as NA, for the genes which fail regression.
  # and contan the gene names as rownames
  n_genes <- ncol(binary_t)
  n_cells <- nrow(binary_t)
  pvalues         <- rep(NA_real_, n_genes)
  pseudoR2s       <- rep(NA_real_, n_genes)
  estimates       <- rep(NA_real_, n_genes)
  switch_at_time  <- rep(NA_real_, n_genes)
  prd_quality     <- rep(0, n_genes) # Keep quality as 0 (numeric)
  CI              <- rep(NA_real_, n_genes)

  # Pre-calculate Time Bounds
  # extract pseudotime data from the sce object
  timedata <- sce$Pseudotime 
  # The loop calls min() and max() on 'timedata' repeatedly. 
  # Since pseudotime is the same for all genes, calculate this once.
  t_min <- min(timedata)
  t_max <- max(timedata)

  # Pre-calculate the Design Matrix
  # Instead of building this matrix ~20,000 times inside the loop, build it once.
  # This matrix contains the Intercept (col 1) and Pseudotime (col 2).
  design_mat_full <- model.matrix(~ timedata)

  # Pre-calculate Downsampling Trigger
  # Instead of doing the 'round()' and boolean checks inside the loop,
  # create a simple TRUE/FALSE vector for all genes at once.
  if (downsample) {
    should_downsample <- round(rowData(sce)$zerop_gene, 3) > ds_cutoff  # not sure why we need to round here, but keeping consistent with original code -MTN 08/02/26
  } else {
    should_downsample <- logical(n_genes) # Initialize a vector of FALSE values for all genes
  }

  # SPARSE CHECK - Skip genes where the regression is unlikely to work due to too few 1s or too few 0s.
  # Calculate the number of cells expressing each gene (Sum of the column in binary_t)
  gene_counts <- Matrix::colSums(binary_t)
  # Define the threshold
  min_cells_required <- 10
  # Identify genes to skip:
  # 1. Expressed in < 10 cells (Too few 1s)
  # 2. Expressed in > (Total - 10) cells (Too few 0s)
  skip_gene <- (gene_counts < min_cells_required) | (gene_counts > (n_cells - min_cells_required))
  # report number of skipped genes
  message(sum(skip_gene), " genes will be skipped due to insufficient expression.")
          # "If a gene is expressed in fewer than 10 cells regression results may be unreliable"

  message("Fitting logistic regression models...")

  for (i in 1:ncol(binary_t)) {
    # Skip genes that are too sparse for reliable regression
    if (skip_gene[i]) {
        next
    }

    # extract the binary state for the current gene
    gene_exp_vec <- as.numeric(binary_t[, i]) 

    #
    glm_response_var <- gene_exp_vec
    glm_predictor_var <- design_mat_full

    # --- DOWNSAMPLING ZEROS---
    # downsample zeros if the percentage of zeros is above the cutoff and downsampling is enabled
    # Downsample inside the loop, but only for genes that meet the downsampling criteria.
    # This methods avoids building a new dataframe for every gene, and instead just slices the existing vectors/matrices.
    if (should_downsample[i]) {
    # We only fetch the raw expression data IF we actually need to downsample
    exp_vec <- as.numeric(exp_t[, i]) 
    
    # Identify zeros
    zero_indices <- which(exp_vec == 0)
    n_non_zero   <- length(exp_vec) - length(zero_indices)
    
    # Calculate number of zeros to keep based on the desired ratio
    n_keep <- round(n_non_zero * (zero_ratio / (1 - zero_ratio)))
    
    # Sample if necessary
    set.seed(42) # 
    keep_zeros <- sample(zero_indices, n_keep)
    
    # Combine indices: (All non-zeros) + (Selected zeros)
    keep_idx <- c(which(exp_vec != 0), keep_zeros)
    
    # Update our working variables by slicing
    glm_response_var <- gene_exp_vec[keep_idx]
    glm_predictor_var <- design_mat_full[keep_idx, , drop = FALSE]
    }

    # fit a logistic regression model using fastglm, 
    # with the binary expression as the response variable and pseudotime as the predictor variable
    # Fit the model using the prepared matrices
    # use trycatch here if errors occur due to non-convergence or other issues with the regression.
    fit <- fastglm(x = glm_predictor_var, y = glm_response_var, family = binomial(link = "logit"))

    # If fit failed (NULL), skip to next gene
    if (is.null(fit)) {
        next
    }

    # --- EXTRACT REGRESSION RESULTS ---
    # # extract information from the regression results.
    # # calculate the time at which the predicted probability of expression is 0.5 using the coefficients from the logistic regression model.
    # # check if the switching time is within the range of the pseudotime data, and adjust it if it is outside the range.
    # # this could be optional, as in some (rare) downstream applications, it might be useful to know  which genes are predicted to switch before or after the observed pseudotime range. 
    # calculate the confidence interval for the switching time using the standard errors of the coefficients from the logistic regression model.

    # Extract Summary ONCE
    sum_fit <- summary(fit)
    coefs   <- coef(sum_fit)
    
    # Check for valid output (needs at least 2 rows: Intercept + Time)
    if (nrow(coefs) < 2) { 
        pvalues[i] <- NA; next 
    }

    #Store Basic Stats
    pvalues[i]   <- coefs[2, 4] # P-value (Row 2, Col 4)
    estimates[i] <- coefs[2, 1] # Slope (Row 2, Col 1)

    #Pseudo R2
    ll_null     <- fit$null.deviance / -2
    ll_proposed <- fit$deviance / -2
    pseudoR2s[i] <- (ll_null - ll_proposed) / ll_null

    #Switch Time Calculation
    # Formula: Time = (TargetLogit - Intercept) / Slope
    # Since TargetLogit for 0.5 is 0, Time = -Intercept / Slope
    intercept_val <- coefs[1, 1]
    slope_val     <- coefs[2, 1]
    
    calc_time <- -intercept_val / slope_val
    
    # Quality Control (Using pre-calculated t_min/t_max)
    # if the calculated switching time is outside the observed pseudotime range,
    #   we set the switching time to the nearest bound (t_min or t_max) and mark the prediction quality as 0 (low).
    if (calc_time >= t_max) {
        switch_at_time[i] <- t_max; prd_quality[i] <- 0
    } else if (calc_time <= t_min) {
        switch_at_time[i] <- t_min; prd_quality[i] <- 0
    } else {
        switch_at_time[i] <- calc_time; prd_quality[i] <- 1
    }

    # Confidence Interval (Delta Method)
    # combines the standard errors of the intercept and slope to estimate the uncertainty in the switching time.
    # if statement to check if calc_time is NA (which can happen if slope_val is 0) 
    # or if the intercept is very close to zero (which can lead to errors in the CI calculation). 
    # If either of these conditions is true, we set the CI to NA to avoid misleading results.
    if (!is.na(calc_time) && abs(intercept_val) > 1e-10) {
        se_intercept <- coefs[1, 2]
        se_slope     <- coefs[2, 2]
        term1 <- (se_intercept * 1.96 / intercept_val) ^ 2
        term2 <- (se_slope * 1.96 / slope_val) ^ 2
        CI[i] <- sqrt(term1 + term2) * abs(intercept_val / slope_val)
    } else {
        CI[i] <- NA
    }
  }

  message("Finished fitting models. Processing results...")
  
# --- POST-PROCESSING ---
# Create Dataframe of Results
# We use 'colnames(binary_t)' because our data was transposed
result_switch <- data.frame(
  switch_at_time = switch_at_time,
  CI             = CI,
  pvalues        = pvalues,
  pseudoR2s      = pseudoR2s,
  estimates      = estimates,
  prd_quality    = prd_quality,
  row.names      = colnames(binary_t) # Ensures gene names match results
)

# Add Derived Columns (Vectorized Operations)
# Use ifelse for "up/down" logic—it's cleaner and faster
result_switch$direction <- ifelse(result_switch$estimates > 0, "up", "down")
result_switch$FDR       <- p.adjust(result_switch$pvalues, method = "BH")

# Use our pre-calculated t_max and t_min constants
steptime <- (t_max - t_min) / 100
result_switch$switch_at_timeidx <- round((result_switch$switch_at_time - t_min) / steptime)

# Quality Filter 
# Identify genes with high FDR. 
# Use which() because it handles NAs gracefully.
# (If FDR is NA, which() simply ignores it, preventing the 'if' crash)
high_fdr_indices <- which(result_switch$FDR > sig_FDR)
if (length(high_fdr_indices) > 0) {
  result_switch$prd_quality[high_fdr_indices] <- 0
}


# --- INTEGRATION ---
# Instead of a slow 'merge', we assign values directly to the SCE object.
# This preserves row order and handles the fact that we filtered genes earlier.

# Extract rowData to a temp variable (modifying S4 objects repeatedly is slow)
r_data <- rowData(sce)

# Initialize columns with NA (or 0) for ALL genes
r_data$switch_at_time     <- NA_real_
r_data$CI                 <- NA_real_
r_data$pvalues            <- NA_real_
r_data$pseudoR2s          <- NA_real_
r_data$estimates          <- NA_real_
r_data$prd_quality        <- 0  # Default quality is 0 (low) until proven otherwise
r_data$direction          <- NA_character_
r_data$FDR                <- NA_real_
r_data$switch_at_timeidx  <- NA_real_


# Inject results into the processed genes only
r_data$switch_at_time    <- result_switch$switch_at_time
r_data$CI                <- result_switch$CI
r_data$pvalues           <- result_switch$pvalues
r_data$pseudoR2s         <- result_switch$pseudoR2s
r_data$estimates         <- result_switch$estimates
r_data$prd_quality       <- result_switch$prd_quality
r_data$direction         <- result_switch$direction
r_data$FDR               <- result_switch$FDR
r_data$switch_at_timeidx <- result_switch$switch_at_timeidx

# Save back to SCE
rowData(sce) <- r_data

return(sce)
}

