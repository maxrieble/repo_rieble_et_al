################################################################################
# SCRIPT 01: Run Hierarchical BACIPS Model (Core Method)
#
# DESCRIPTION:
# Executes the core Bayesian analysis using 'BACIPS.R', taking
# 'TAB_DEL_*.rdata' as input. Saves raw model outputs ('out*.rdata') and the
# raw JAGS summary ('All_raw_*.rdata').
# WARNING: Takes several hours.
#
# INPUTS: Results/TAB_DEL_*.rdata, Scripts/BACIPS.R
# OUTPUTS: Results/out*.rdata, Results/All_raw_*.rdata
################################################################################

# --- SETUP ---
message("Starting Script 01: Run BACIPS Model (Core Method)")
# setwd("D:/Daten/Paper/Repository") # Set repo root if running standalone

# --- CONFIGURATION ---
species_list <- c("Perch", "Carp")

# --- PROCESSING LOOP ---
for (sp in species_list) {
  message(paste("\n>>> Running BACIPS Model for species:", sp))
  prefix <- substr(sp, 1, 1)
  
  # Load Model Input
  delta_path <- file.path("Results", paste0("TAB_DEL_", prefix, ".rdata"))
  if (!file.exists(delta_path)) {
    warning("FATAL: Delta table not found for ", sp, ". Cannot proceed.")
    next
  }
  load(delta_path) # Loads 'TABLE_DATA_DELTA'
  message("...Loaded delta table.")
  
  # Run JAGS Model
  species_num <- switch(sp, "Perch" = 1, "Carp" = 2)
  message("...Calling BACIPS.R (JAGS model). This takes several hours...")
  source(file.path("Scripts", "BACIPS.R")) # Creates 'out' and raw 'All' summary
  message("...JAGS run complete.")
  
  # Save Raw Outputs
  out_path <- file.path("Results", paste0("out", prefix, ".rdata"))
  all_raw_path <- file.path("Results", paste0("All_raw_", prefix, ".rdata"))
  
  save(out, file = out_path)
  save(All, file = all_raw_path) # 'All' here is just out$BUGSoutput$summary
  message("...Saved raw 'out' object to: ", basename(out_path))
  message("...Saved raw 'All' summary object to: ", basename(all_raw_path))
}
message("\n✅ Script 01 finished.")