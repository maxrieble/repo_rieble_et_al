################################################################################
# SCRIPT Supp_A: Cluster Analysis and Data Merging (Documentation)
#
# DESCRIPTION:
# Documents the Random Forest clustering process, LD/MD calculation,
# merging of parameters with metadata, and meta-model analysis used to
# create the final 'All*.rdata' files. This script shows how the pre-computed
# final 'All' files (used by 02_Generate_Plots.R) were generated.
#
# INPUTS:
# - Results/All_raw_C.rdata, Results/All_raw_P.rdata (RAW JAGS summary)
# - Input/CarpCatch.rdata, Input/PerchCatch.rdata (must contain env. data?)
# - Results/delta.3.pred.C, Results/delta.3.pred.P (or Scripts/Delta_Predictions.R)
# - Results/TAB_DEL_C.rdata, Results/TAB_DEL_P.rdata
#
# OUTPUTS (if run):
# - Results/Clustnames.rdata (Cluster assignments)
# - Results/AllC.rdata, Results/AllP.rdata (FINAL combined table)
# - Console output (Meta-model summaries).
################################################################################

# --- 1. SETUP ---
# setwd("D:/Daten/Paper/Repository") # Set repo root if running standalone
library(dplyr)
library(randomForest)
library(cluster)
library(nnet) # For multinom()

message("Starting Supplementary Cluster Analysis and Merging Script (Supp_A)")

# --- 2. CONFIGURATION ---
species_list <- c("Perch", "Carp")
set.seed(90) # Set seed for reproducibility of RF clustering

# --- 3. PROCESSING LOOP ---
Clustnames <- list() # Initialize list to store cluster results for saving

for (sp in species_list) {
  message(paste("\n", paste(rep("-", 40), collapse = "")))
  message(paste(">>> Clustering and Merging for species:", sp))
  message(paste(rep("-", 40), collapse = ""))
  
  prefix <- substr(sp, 1, 1)
  
  # --- 3.1. LOAD REQUIRED DATA ---
  message("... Loading required data.")
  
  # Load Raw JAGS summary (output from Script 01)
  all_raw_path <- file.path("Results", paste0("All_raw_", prefix, ".rdata"))
  if (!file.exists(all_raw_path)) {
    stop(paste("FATAL: Raw 'All' summary not found:", basename(all_raw_path), ". Run Script 01 first."))
  }
  load(all_raw_path); all_summary_raw <- All; rm(All) # Loads raw 'All' summary
  
  # Load CATCH data (check for required columns later)
  catch_path <- file.path("Input", paste0(sp, "Catch.rdata"))
  if (!file.exists(catch_path)) { stop(paste("FATAL: CATCH file not found:", basename(catch_path))) }
  load(catch_path); CATCH <- if (sp == "Carp") CarpCatch2 else PerchCatch2; n_fish <- nrow(CATCH)
  
  # Load or Generate predicted delta values ('delta.3.pred')
  pred_file <- file.path("Results", paste0("delta.3.pred.", prefix, ".rdata"))
  if (!file.exists(pred_file)) {
    delta_script_path <- file.path("Scripts","Delta_Predictions.R")
    if(file.exists(delta_script_path)){
      message("... delta.3.pred file not found. Running Delta_Predictions.R...")
      # Load 'out' object needed by Delta_Predictions.R
      out_path <- file.path("Results", paste0("out", prefix, ".rdata"))
      if(!file.exists(out_path)) stop(paste("Cannot run Delta_Predictions.R: out file missing:", basename(out_path)))
      load(out_path) # Loads 'out' or 'outC'/'outP'
      out_obj_name <- if(exists("out", inherits = FALSE)) "out" else paste0("out", prefix)
      if(!exists(out_obj_name, inherits = FALSE)) stop("Could not find 'out' object in loaded file.")
      assign("out", get(out_obj_name), envir = .GlobalEnv) # Ensure it's named 'out' globally
      
      # Load TABLE_DATA_DELTA for Standardize value
      delta_path <- file.path("Results", paste0("TAB_DEL_", prefix, ".rdata"))
      if (!file.exists(delta_path)) stop("Cannot run Delta_Predictions.R: TAB_DEL file missing.")
      load(delta_path) # Loads TABLE_DATA_DELTA
      # Define variables needed by Delta_Predictions.R in the global environment
      assign("ind", n_fish, envir = .GlobalEnv)
      assign("time.after", 1:10, envir = .GlobalEnv)
      assign("Standardize_val", -min(TABLE_DATA_DELTA$DELTA, na.rm = TRUE) + 1000, envir = .GlobalEnv)
      
      # Source the script
      source(delta_script_path) # Creates delta.3.pred globally
      if(!exists("delta.3.pred")) stop("Delta_Predictions.R did not create delta.3.pred object.")
      message("... delta.3.pred generated.")
      # Clean up global variables set for the sourced script
      rm(out, ind, time.after, Standardize_val, TABLE_DATA_DELTA, envir = .GlobalEnv)
    } else {
      stop(paste("FATAL: Cannot find required prediction file", basename(pred_file), "or script", basename(delta_script_path)))
    }
  } else {
    load(pred_file, .GlobalEnv) # Loads delta.3.pred object
    message("... Loaded pre-computed delta.3.pred.")
  }
  
  # Load TABLE_DATA_DELTA again (if cleaned up) just for Standardize value
  delta_path <- file.path("Results", paste0("TAB_DEL_", prefix, ".rdata"))
  if (!exists("TABLE_DATA_DELTA")) load(delta_path)
  Standardize_val <- -min(TABLE_DATA_DELTA$DELTA, na.rm = TRUE) + 1000
  
  # --- 3.2. EXTRACT PARAMETERS & CALCULATE LD/MD ---
  message("... Extracting median parameters and calculating LD/MD.")
  
  # Extract median parameter estimates from raw summary
  param_medians <- data.frame(ID = 1:n_fish) # BUGS model index (1 to n_fish)
  params_to_extract <- c("alpha.i", "beta.i", "rr.i", "kone.i", "ktwo.i")
  for(param in params_to_extract) {
    param_medians[[param]] <- sapply(1:n_fish, function(i) {
      rowname <- paste0(param, "[", i, "]")
      # Check if the parameter exists in the raw summary
      if(rowname %in% rownames(all_summary_raw)) {
        return(all_summary_raw[rowname, "50%"]) # Get the median estimate
      } else {
        warning(paste("Parameter", rowname, "not found in raw All summary for", sp))
        return(NA)
      }
    })
  }
  
  # Calculate mean predicted trajectory (V1-V10) from delta.3.pred
  PredMeans <- matrix(NA, nrow = n_fish, ncol = 10); colnames(PredMeans) <- paste0("V", 1:10)
  expected_cols <- as.character(1:n_fish)
  if(!all(expected_cols %in% colnames(delta.3.pred))) stop("Column names in delta.3.pred do not match expected fish IDs (1:n_fish).")
  for (ID in 1:n_fish) {
    for (Step in 1:10) {
      # Select column by character name "1", "2", etc.
      pred_vals <- delta.3.pred[delta.3.pred$timestep == Step, as.character(ID)]
      if(length(pred_vals) > 0) PredMeans[ID, Step] <- mean(pred_vals, na.rm = TRUE)
    }
  }
  
  # Calculate LD/MD (difference on the standardized TDelta2 scale)
  # LD: Lowest point of mean trajectory relative to baseline (kone.i)
  # MD: Mean point of mean trajectory relative to baseline (kone.i)
  LD <- numeric(n_fish); MD <- numeric(n_fish)
  for (ID in 1:n_fish) {
    kone_val <- param_medians$kone.i[ID] # Baseline on standardized scale
    pred_vals <- PredMeans[ID, ]        # Mean trajectory on standardized scale
    LD[ID] <- tryCatch(kone_val - min(pred_vals, na.rm=TRUE), error=function(e) NA) # kone.i is higher -> positive LD
    MD[ID] <- tryCatch(kone_val - mean(pred_vals, na.rm=TRUE), error=function(e) NA) # kone.i is higher -> positive MD
  }
  
  # Assemble data frame for clustering input
  DataForClustering <- param_medians %>%
    mutate(LD = LD, MD = MD)
  
  # --- 3.3. PERFORM CLUSTER ANALYSIS ---
  message("... Performing Random Forest clustering.")
  
  # Select predictors for RF clustering
  # [# NOTE:] Predictors chosen include BACIPS parameters (alpha.i, beta.i, rr.i, kone.i)
  # and derived summary metrics (LD, MD). ktwo.i was excluded in some original scripts.
  predictors_for_rf <- c("alpha.i", "beta.i", "rr.i", "kone.i", "LD", "MD")
  ScaledData <- DataForClustering
  
  # Check for NA values before scaling
  na_check <- sapply(ScaledData[predictors_for_rf], function(x) sum(is.na(x)))
  if(any(na_check > 0)) {
    warning(paste("NA values found in RF predictors for", sp, ":",
                  paste(names(na_check[na_check>0]), collapse=", "),
                  ". RF will omit rows with NAs."))
  }
  # Scale predictors used in RF
  ScaledData[predictors_for_rf] <- lapply(ScaledData[predictors_for_rf], function(x) as.numeric(scale(x)))
  
  # Apply Carp exclusion TEMPORARILY before clustering
  original_indices <- 1:n_fish # Track original BUGS index
  indices_to_cluster <- original_indices
  scaled_data_for_rf <- ScaledData # Work with a copy
  
  if (sp == "Carp") {
    exclude_index <- which(param_medians$ID == 10) # Find row index for BUGS ID 10
    if(length(exclude_index) > 0) {
      scaled_data_for_rf <- ScaledData[-exclude_index, ] # Exclude row from RF input
      indices_to_cluster <- original_indices[-exclude_index] # Keep track of original IDs included
      message("... Temporarily removing Carp ID 10 for clustering.")
    }
  }
  
  # Run Random Forest to get proximity matrix
  rf.fit <- randomForest(
    x = scaled_data_for_rf[, predictors_for_rf], # Use scaled data without excluded fish
    y = NULL,               # Unsupervised clustering
    ntree = 10000,          # Number of trees
    proximity = TRUE,       # Calculate proximity matrix
    oob.prox = TRUE,        # Use out-of-bag samples for proximity
    na.action = na.omit     # Omit rows with NA predictors
  )
  
  # Hierarchical clustering on the RF proximity matrix
  hclust.rf <- hclust(as.dist(1 - rf.fit$proximity), method = "ward.D2")
  
  # Cut the dendrogram into clusters
  # [# NOTE:] k=4 clusters were chosen based on visual inspection of the dendrogram
  # and subsequent biological interpretation.
  k_clusters <- 4
  rf_cluster_results <- cutree(hclust.rf, k = k_clusters)
  
  # Store cluster assignments, mapping back to original fish IDs
  cluster_vector <- rep(NA, n_fish) # Full vector including potentially excluded fish
  # Ensure na.action did not remove unexpected rows
  if(nrow(na.omit(scaled_data_for_rf[, predictors_for_rf])) != length(rf_cluster_results)){
    stop(paste("FATAL: Mismatch after na.omit in RF for", sp, ". Check NA handling."))
  }
  # Find which original indices correspond to the rows kept by RF
  kept_indices <- indices_to_cluster[!apply(ScaledData[indices_to_cluster, predictors_for_rf], 1, anyNA)]
  if(length(rf_cluster_results) == length(kept_indices)) {
    cluster_vector[kept_indices] <- rf_cluster_results
  } else {
    stop(paste("FATAL: Mismatch between RF cluster results and non-NA fish indices for", sp))
  }
  Clustnames[[sp]] <- cluster_vector # Store results for this species
  
  # --- 3.4. ASSEMBLE FINAL 'All' OBJECT ---
  message("... Assembling final 'All' object with merged data.")
  
  # Start with median parameters
  All_final <- param_medians
  
  # Add identifiers and TL from CATCH (assumes CATCH rows match ID 1:n_fish order)
  All_final$Lotek.ID <- CATCH$Lotek.ID
  All_final$TL <- CATCH$TL # Assumes TL column exists
  All_final$CatchNr <- CATCH$CatchNr
  
  # Add predicted mean trajectory (V1-V10)
  All_final <- cbind(All_final, as.data.frame(PredMeans))
  
  # Add Cluster number (vector saved above)
  All_final$Cluster <- Clustnames[[sp]]
  
  # Add RankCluster (Manual ranking - requires adjustment based on inspection)
  # [# NOTE:] Define the mapping from numeric Cluster (1-4) to RankCluster labels.
  # This mapping MUST be confirmed by inspecting cluster profiles (e.g., mean trajectories).
  if (sp == "Perch") {
    # Example based on previous script: Cluster 1->low, 2->mid, 3->high1, 4->high2
    rank_map <- c("low", "mid", "high1", "high2")
    rank_levels <- c("low", "mid", "high1", "high2")
  } else { # Carp
    # Example based on previous script: Cluster 1->low, 2->mid, 3->high1, 4->high2
    # [# NOTE:] ** ADJUST THIS MAPPING BASED ON YOUR ACTUAL CARP RESULTS **
    rank_map <- c("low", "mid", "high1", "high2")
    rank_levels <- c("low", "mid", "high1", "high2")
  }
  # Initialize with NAs, then fill for non-NA clusters
  All_final$RankCluster <- factor(NA, levels = rank_levels)
  valid_clusters <- !is.na(All_final$Cluster)
  if(any(All_final$Cluster[valid_clusters] > length(rank_map) | All_final$Cluster[valid_clusters] < 1)){
    warning(paste("Invalid cluster numbers detected for", sp,". Check Clustnames generation."))
  } else {
    All_final$RankCluster[valid_clusters] <- factor(rank_map[All_final$Cluster[valid_clusters]], levels = rank_levels)
  }
  
  
  # Add Environmental Data from CATCH
  # [# NOTE:] This assumes relevant columns (Temperature, Turbidity, ODO, ODOsat)
  # exist in the CATCH object. Add warnings if they are missing.
  env_cols <- c("Temperature", "Turbidity", "ODO", "ODOsat")
  for(col in env_cols) {
    if(col %in% names(CATCH)) {
      All_final[[col]] <- CATCH[[col]]
    } else {
      All_final[[col]] <- NA # Add NA column if data is missing
      warning(paste("Environmental variable '", col, "' not found in CATCH data for ", sp))
    }
  }
  
  # Add Condition Factor
  # [# NOTE:] Calculates from weight/length if available in CATCH, otherwise adds NA.
  if("CondFactor" %in% names(CATCH)) {
    All_final$CondFactor <- CATCH$CondFactor
  } else if (all(c("Weight..g.", "Total.Length..mm.") %in% names(CATCH))) {
    # Use the specific column names found in CarpCatch2 printout
    All_final$CondFactor <- 100000 * (CATCH$`Weight..g.` / CATCH$`Total.Length..mm.` ^ 3)
  } else {
    All_final$CondFactor <- NA
    warning(paste("Could not calculate or find CondFactor for ", sp))
  }
  
  # Add LD and MD metrics
  All_final$LD <- LD
  All_final$MD <- MD
  
  # --- 3.5. FINAL CARP EXCLUSION ---
  # Apply exclusion permanently to the final merged table
  if (sp == "Carp") {
    exclude_index <- which(All_final$ID == 10) # Find row by BUGS ID
    if(length(exclude_index) > 0) {
      All_final <- All_final[-exclude_index, , drop = FALSE]
      message("... Permanently removed Carp ID 10 from final 'All' object.")
    } else {
      # This might happen if ID 10 had NAs and was omitted by RF na.action
      warning("Carp exclusion: Fish ID 10 not found in the final 'All' object (possibly removed due to NAs).")
    }
  }
  
  # --- 3.6. SAVE FINAL 'All' OBJECT ---
  # This is the object loaded by 02_Generate_Plots.R
  All <- All_final # Assign to standard 'All' variable name for saving
  all_path <- file.path("Results", paste0("All", prefix, ".rdata"))
  save(All, file = all_path)
  message("... Saved FINAL combined 'All' object to: ", basename(all_path))
  
  # Clean up large objects from global environment before next loop iteration
  if(exists("delta.3.pred")) rm(delta.3.pred, envir = .GlobalEnv)
  if(exists("All")) rm(All) # Remove the species-specific All from global env
  
} # End species loop

# --- 4. SAVE CLUSTNAMES ---
# Save the list containing cluster assignments for both species
clust_path <- file.path("Results", "Clustnames.rdata")
save(Clustnames, file = clust_path)
message("\n... Saved cluster assignments for all species to: ", basename(clust_path))

# --- 5. META-MODEL ANALYSIS (Optional, Documentation) ---
message("\nRunning Meta-Model logic (GLM, Multinomial) for documentation...")
# [# NOTE:] This section documents the meta-models run on the final 'All' data.
# Results depend on the availability of environmental predictors in the final 'All' files.
meta_results <- list() # Store results

for (sp in species_list) {
  message(paste("-- Meta-models for:", sp, "--"))
  prefix <- substr(sp, 1, 1)
  
  # Load the FINAL merged 'All' file just created
  load(file.path("Results", paste0("All", prefix, ".rdata"))) # Loads 'All'
  
  # Prepare scaled data frame for models
  AllScaled <- All
  # Define potential predictors present in the final 'All' file
  predictors_for_meta <- c("TL", "Temperature", "Turbidity", "ODO") # Add CondFactor if desired/available
  # Check which predictors are actually available (not all NA)
  predictors_present <- intersect(predictors_for_meta, names(AllScaled))
  predictors_present <- predictors_present[sapply(AllScaled[predictors_present], function(x) !all(is.na(x)))]
  
  if (length(predictors_present) > 0) {
    message(paste("... Using predictors:", paste(predictors_present, collapse=", ")))
    # Scale the available predictors
    AllScaled[predictors_present] <- lapply(AllScaled[predictors_present], function(x) as.numeric(scale(x)))
    
    # Define model formulas dynamically based on available predictors
    formula_rhs <- paste(predictors_present, collapse=" + ")
    formula_glm_md <- as.formula(paste("MD ~", formula_rhs))
    formula_glm_ld <- as.formula(paste("LD ~", formula_rhs))
    formula_multinom <- as.formula(paste("RankCluster ~", formula_rhs))
    
    # Run GLM models for MD and LD
    MetaFit_MD <- tryCatch(glm(formula_glm_md, data = AllScaled, na.action = "na.exclude"), error=function(e) {message("MD GLM failed:", e$message); NULL})
    MetaFit_LD <- tryCatch(glm(formula_glm_ld, data = AllScaled, na.action = "na.exclude"), error=function(e) {message("LD GLM failed:", e$message); NULL})
    
    # Run Multinomial model for RankCluster
    ClusterRankmodel <- tryCatch(nnet::multinom(formula_multinom, data = AllScaled, na.action = "na.exclude", maxit=1000, trace=FALSE), error=function(e) {message("Multinomial model failed:", e$message); NULL})
    
    # Store and print results
    meta_results[[sp]] <- list(MD_GLM=MetaFit_MD, LD_GLM=MetaFit_LD, Cluster_Multinom=ClusterRankmodel)
    
    if(!is.null(MetaFit_MD)) { message("-- MD GLM Summary --"); print(summary(MetaFit_MD)) }
    if(!is.null(MetaFit_LD)) { message("-- LD GLM Summary --"); print(summary(MetaFit_LD)) }
    if(!is.null(ClusterRankmodel)) {
      message("-- Multinomial Model Summary --"); print(summary(ClusterRankmodel))
      # Calculate p-values using Wald z-tests
      z <- summary(ClusterRankmodel)$coefficients / summary(ClusterRankmodel)$standard.errors
      p <- (1 - pnorm(abs(z), 0, 1)) * 2
      message("-- p-values for Multinomial Model Coefficients --"); print(p)
    }
  } else {
    message("... Skipping meta-models: No non-NA environmental predictors found in final 'All' file.")
    meta_results[[sp]] <- list(MD_GLM=NULL, LD_GLM=NULL, Cluster_Multinom=NULL)
  }
} # End meta-model loop

message("\n✅ Supp_A script finished.")