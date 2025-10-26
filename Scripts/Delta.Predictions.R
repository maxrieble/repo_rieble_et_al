################################################################################
# Delta_Predictions.R (Documentation)
#
# DESCRIPTION:
# Calculates the full posterior predicted delta trajectories ('delta.3.pred')
# from the raw JAGS output ('out' object). Saves the result to a file.
#
# INPUTS (expected in global environment when sourced):
# - `out`: Raw JAGS output object (from out*.rdata).
# - `TABLE_DATA_DELTA`: Needed for Standardize value.
# - `ind`: Number of fish/capture events (integer).
# - `time.after`: Time vector, typically 1:10 (numeric vector).
# - `Standardize_val`: The value used to standardize DELTA (numeric).
# - `sp`: Species name ("Carp" or "Perch") (optional, for filename).
#
# OUTPUTS:
# - Creates 'delta.3.pred' data frame in the global environment.
# - Saves 'delta.3.pred.*.rdata' file to Results/ folder.
################################################################################
# [# NOTE:] This script is typically called by Supp_A... script if the
# pre-computed delta.3.pred.*.rdata file is missing.

message("Starting Delta_Predictions.R script...")

# --- Input Validation ---
required_vars <- c("out", "TABLE_DATA_DELTA", "ind", "time.after", "Standardize_val")
if (!all(required_vars %in% ls(.GlobalEnv))) {
  missing_vars <- required_vars[!required_vars %in% ls(.GlobalEnv)]
  stop("Delta_Predictions.R requires variables to be defined in the global environment: ",
       paste(missing_vars, collapse=", "))
}
if (!"BUGSoutput" %in% names(out) || !"sims.list" %in% names(out$BUGSoutput)) {
  stop("'out' object does not have the expected JAGS output structure.")
}

# --- Initialization ---
delta.3.pred <- data.frame() # Initialize empty data frame
n_iter_total <- out$BUGSoutput$n.keep # n.keep already incorporates chains*thin*keep
n_fish <- ind                       # Use 'ind' passed from calling script

message(paste("... Calculating predictions for", n_fish, "individuals across", n_iter_total, "MCMC samples..."))

# Check dimensions of MCMC samples match n_fish
mcmc_params <- c("kone.i", "rr.i", "alpha.i", "beta.i", "ktwo.i")
for(p in mcmc_params) {
  if(!p %in% names(out$BUGSoutput$sims.list) || ncol(out$BUGSoutput$sims.list[[p]]) != n_fish) {
    stop(paste("Mismatch between n_fish (", n_fish, ") and columns in out$BUGSoutput$sims.list$", p))
  }
}

# --- Prediction Loop ---
# Iterate through each saved MCMC sample
for (PredIndex in 1:n_iter_total) {
  # Matrix to store predictions for this iteration (rows=time, cols=fish)
  delta.2.pred_matrix <- matrix(NA, nrow = length(time.after), ncol = n_fish)
  
  # Iterate through each fish
  for (IdIndex in 1:n_fish) {
    Pred <- numeric(length(time.after)) # Vector for this fish's trajectory
    
    # Safely get parameters for this iteration and fish
    # [# NOTE:] Parameter names match those in the BACIPS.R JAGS model
    kone_i  = out$BUGSoutput$sims.list$kone.i[PredIndex, IdIndex]
    rr_i    = out$BUGSoutput$sims.list$rr.i[PredIndex, IdIndex]
    alpha_i = out$BUGSoutput$sims.list$alpha.i[PredIndex, IdIndex]
    beta_i  = out$BUGSoutput$sims.list$beta.i[PredIndex, IdIndex]
    ktwo_i  = out$BUGSoutput$sims.list$ktwo.i[PredIndex, IdIndex]
    
    # Skip calculation if any parameter is missing/NA for this iteration/fish
    if(any(is.na(c(kone_i, rr_i, alpha_i, beta_i, ktwo_i)))) {
      # Assign NA to the entire trajectory for this fish in this iteration
      delta.2.pred_matrix[, IdIndex] <- NA
      next # Skip to the next fish
    }
    
    # Calculate trajectory step-by-step using the BACIPS recurrence relation
    for (StepIndex in 1:length(time.after)) {
      current_time <- time.after[StepIndex] # Use actual time value (1 to 10)
      
      if (current_time == 1) {
        # Formula for the first step after impact
        Pred[StepIndex] <- kone_i + rr_i * kone_i * (1 - (kone_i / ktwo_i)) -
          (alpha_i / (1 + beta_i * current_time)) * kone_i
      } else {
        # Formula for subsequent steps, depends on the previous step's prediction
        Pred[StepIndex] <- Pred[StepIndex - 1] + rr_i * Pred[StepIndex - 1] * (1 - (Pred[StepIndex - 1] / ktwo_i)) -
          (alpha_i / (1 + beta_i * current_time)) * Pred[StepIndex - 1]
      }
      # Handle potential non-finite results (Inf, -Inf, NaN)
      if (!is.finite(Pred[StepIndex])) {
        Pred[StepIndex] <- NA
      }
    } # End time step loop
    # Store the calculated trajectory for this fish
    delta.2.pred_matrix[, IdIndex] <- Pred
  } # End fish loop
  
  # Undo the standardization applied before the model run
  # [# NOTE:] The predictions are generated on the TDELTA2 scale and need
  # to be converted back to the original DELTA scale.
  delta.2.pred_matrix <- delta.2.pred_matrix - Standardize_val
  
  # Convert matrix to data frame for easier handling
  delta.2.pred_df <- as.data.frame(delta.2.pred_matrix)
  # Name columns by fish ID (1 to n_fish)
  colnames(delta.2.pred_df) <- as.character(1:n_fish)
  
  # Add timestep and iteration information
  delta.2.pred_df$timestep <- time.after
  delta.2.pred_df$iteration <- PredIndex
  
  # Append the results for this iteration to the main data frame
  delta.3.pred <- rbind(delta.3.pred, delta.2.pred_df)
  
  # Optional progress update
  if (PredIndex %% 1000 == 0) { message(paste("... completed iteration", PredIndex, "of", n_iter_total)) }
  
} # End MCMC iteration loop

# --- Save Output ---
# Determine filename based on species variable if it exists
prefix <- "Unknown"
if(exists("sp", inherits = FALSE)) {
  prefix <- substr(sp, 1, 1)
} else if (exists("species_name", inherits=FALSE)) {
  prefix <- substr(species_name, 1, 1)
}

output_filename <- file.path("Results", paste0("delta.3.pred.", prefix, ".rdata"))
save(delta.3.pred, file = output_filename)

message(paste("... delta.3.pred calculation complete. Saved to:", output_filename))

# Clean up large object if running within another script
# rm(delta.2.pred_matrix, delta.2.pred_df)