# =========================================
# SCRIPT 02: Generate Manuscript Figures
#
# DESCRIPTION:
# Loads FINAL pre-computed model results ('All*.rdata', 'out*.rdata') and
# supplementary data to generate the final manuscript figures, including
# Population Activity, Cluster Plots (using base R), and LD/MD Boxplots.
# This is the main script reviewers need to run to reproduce figures.
#
# INPUTS:
# - Results/AllC.rdata, Results/AllP.rdata (FINAL merged tables)
# - Results/outC.rdata, Results/outP.rdata (Raw JAGS output)
# - Results/TAB_DEL_C.rdata, Results/TAB_DEL_P.rdata
# - Input/CarpCatch.rdata, Input/PerchCatch.rdata
# - Input/Cbehaviors2.rdata, Input/Pbehaviors2.rdata
#
# OUTPUTS: Figures saved to /Figures directory.

# !!! IMPORTANT !!!
# This script requires the large raw JAGS output files ('outC.rdata' and
# 'outP.rdata') which are not stored directly in this GitHub repository due
# to size limits. Please download them using the links provided in the
# main README.md file and place them in the 'Results/' folder before running
# this script.
# =========================================

# ---- 1. Setup and Directory Configuration ----
setwd("D:/Daten/Paper/Repository") # Set to repository root

OUTPUT_FOLDER <- "Figures/" # Use "Figures" folder
if (!dir.exists(OUTPUT_FOLDER)) {
  message(paste("Creating output directory:", OUTPUT_FOLDER))
  dir.create(OUTPUT_FOLDER, recursive = TRUE)
} else {
  message(paste("Output directory already exists:", OUTPUT_FOLDER))
}

# Load necessary libraries
library(broom); library(mice); library(miceadds); library(png)
library(fasttime); library(ggplot2); library(gridExtra); library(rjags)
library(dplyr); library(scales); library(stringr) # Added stringr if needed

# ---- 2. Helper: Load species-specific data ----
# Loads FINAL All*.rdata (merged), raw out*, TABLE_DATA_DELTA, CATCH, Behaviors
load_species_data <- function(species_name) {
  out <- list()
  load_env <- new.env() # Load into temporary environment
  
  if(species_name == "Carp") {
    prefix <- "C"; catch_obj_name <- "CarpCatch2"; beh_obj_name <- "Cbehaviors2"
  } else if(species_name == "Perch") {
    prefix <- "P"; catch_obj_name <- "PerchCatch2"; beh_obj_name <- "Pbehaviours2"
  } else stop("Unknown species")
  
  # Load required files
  load(file.path("Results", paste0("All", prefix, ".rdata")), envir = load_env)
  out$All <- load_env$All
  
  load(file.path("Results", paste0("TAB_DEL_", prefix, ".rdata")), envir = load_env)
  out$TABLE_DATA_DELTA <- load_env$TABLE_DATA_DELTA
  
  # Handle potential outC/outP naming in the .rdata file
  out_file_path <- file.path("Results", paste0("out", prefix, ".rdata"))
  out_loaded_names <- load(out_file_path, envir = load_env)
  out_obj_name <- out_loaded_names[1] # Assume first object is the one we want
  out$out <- load_env[[out_obj_name]]
  
  load(file.path("Input", paste0(species_name, "Catch.rdata")), envir = load_env)
  if(!catch_obj_name %in% ls(envir=load_env)){ stop(paste("Object", catch_obj_name, "not found in Catch file."))}
  out$CATCH <- load_env[[catch_obj_name]]
  
  beh_file_path <- file.path("Input", paste0(beh_obj_name, ".rdata")) # Load behaviors file by its direct name
  if(!file.exists(beh_file_path)){stop(paste("Behaviors file not found:", beh_file_path))}
  load(beh_file_path, envir = load_env)
  if(!beh_obj_name %in% ls(envir=load_env)){ stop(paste("Object", beh_obj_name, "not found in Behaviors file."))}
  out$Behaviors <- load_env[[beh_obj_name]]
  
  rm(load_env) # Clean up
  
  # Validation
  expected_cols <- c("ID", "Lotek.ID", "Cluster", "CatchNr", "RankCluster", "TL", "LD", "MD")
  if (!all(expected_cols %in% names(out$All))) {
    warning(paste("Loaded 'All' for", species_name, "missing expected metadata."))
  }
  return(out)
}

# ---- 3. Load species data ----
message("Loading all data for plotting...")
carp_data    <- load_species_data("Carp")
perch_data <- load_species_data("Perch")
species_data <- list(Carp = carp_data, Perch = perch_data)
message("...Data loaded.")

# Define common time variables
time.after  <- seq(1, 10, 1)
time.before <- seq(-5, -1, 1)
time.true   <- c(time.before, time.after)

# ---- 4. PLOT 1: Population Activity ----
message("Generating Plot 1: Population Activity...")

PopAct_plot_list <- list()
species_order <- c("Carp", "Perch")
panel_labels <- c("A", "B")

for (i in seq_along(species_order)) {
  species_name <- species_order[i]
  data <- species_data[[species_name]]
  Beh <- data$Behaviors
  CATCH <- data$CATCH
  
  # Ensure Behaviors data has necessary columns
  if (!all(c("Day2", "Daynum", "distance") %in% names(Beh))) {
    warning(paste("Skipping Pop Activity for", species_name, "- Behaviors missing Day2/Daynum/distance."))
    next
  }
  # Ensure CATCH data has necessary columns
  if (!"Time" %in% names(CATCH) || !inherits(CATCH$Time, "POSIXct")) {
    warning(paste("Skipping Pop Activity for", species_name, "- CATCH missing valid 'Time' column."))
    next
  }
  
  
  # Calculate daily mean distance
  # Convert Day2 to Date if it's not already
  if (!inherits(Beh$Day2, "Date")) {
    Beh$Day2 <- tryCatch(as.Date(fastPOSIXct(Beh$Day2, tz="GMT")), error = function(e) as.Date(Beh$Day2))
    if (!inherits(Beh$Day2, "Date")) {
      warning(paste("Could not convert Day2 to Date for", species_name))
      next
    }
  }
  
  PopAct <- Beh %>%
    filter(!is.na(distance)) %>%
    group_by(Day2) %>%
    summarise(dists = mean(distance, na.rm = TRUE), .groups = 'drop') %>%
    rename(days = Day2) # Rename Date column
  
  # Get catch days (as Date objects)
  capdays <- as.Date(CATCH$Time)
  capdays_min <- min(capdays, na.rm=TRUE)
  capdays_max <- max(capdays, na.rm=TRUE)
  
  
  p <- ggplot(PopAct, aes(x = days, y = dists)) +
    geom_point(size = 1.5, colour = "black") +
    stat_smooth(method = "loess", colour = "blue", se = FALSE) +
    # Add rect only if capdays are valid
    { if (!is.na(capdays_min) && !is.na(capdays_max)) {
      annotate("rect", xmin = capdays_min, xmax = capdays_max, ymin = 0, ymax = Inf, alpha = 0.2, fill = "red")
    }
    } +
    theme_bw(base_size = 20) + # Using specific theme for this plot
    labs(x = "Time", y = "Mean daily swimming distance (m)") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      axis.text = element_text(size = 16),
      axis.title.x = element_text(size = 18, face = "plain"),
      axis.title.y = element_text(size = 17, face = "plain"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    annotate("text", x = max(PopAct$days, na.rm=T), y = max(PopAct$dists, na.rm = TRUE),
             label = panel_labels[i], hjust = 1.1, vjust = 1.1, size = 8, fontface = "bold")
  
  PopAct_plot_list[[species_name]] <- p
}

# Save the combined Population Activity plot
if(length(PopAct_plot_list) > 0) {
  out_path_popact <- file.path(OUTPUT_FOLDER, "Population_Activity_Combined.png")
  message(paste("Saving population activity plot to:", out_path_popact))
  png(out_path_popact, width = 50, height = 14.5, units = "cm", res = 300, type = "cairo")
  grid.arrange(grobs = PopAct_plot_list, ncol = length(PopAct_plot_list))
  dev.off()
} else {
  message("Skipped saving Population Activity plot.")
}

# ---- 5. PLOT 2: Cluster Plots (Base R) ----
# [# NOTE:] This uses the base R plotting code from your script.
message("Generating Plot 2: Cluster Plots (Base R)...")

out_path_clusters <- file.path(OUTPUT_FOLDER, "Clusters_Combined.png")
png(out_path_clusters, width=55, height=29, units="cm", res=300, pointsize=12, type="cairo")

par(mfrow=c(2,4), mar=c(5,5,4,2), cex.lab=1.7, cex.axis=1.5)

species_order <- c("Carp","Perch")
cluster_sizes <- list(
  Carp = c("low n = 3", "mid n = 2", "high1 n = 3", "high2 n = 2"), # Adjust n based on final All file
  Perch = c("low n = 4", "mid n = 4", "high1 n = 6", "high2 n = 6") # Adjust n based on final All file
)

for(species_name in species_order){
  data <- species_data[[species_name]]
  out_sp <- data$out
  All_sp <- data$All # This is the final merged 'All'
  TABLE_sp <- data$TABLE_DATA_DELTA
  
  # Check if 'All' has cluster info
  if (!"Cluster" %in% names(All_sp)) {
    warning(paste("Skipping Base R Cluster plot for", species_name, "- 'All' object missing 'Cluster' column."))
    # Need to draw empty plots to maintain grid structure
    for(k in 1:4) plot(1, type="n", axes=FALSE, xlab="", ylab="", main=paste(species_name, "- Cluster", k, "(Data Missing)"))
    next
  }
  
  Standardize <- -min(TABLE_sp$DELTA, na.rm=TRUE) + 1000
  # Robust calculation of y-limits
  y_min <- suppressWarnings(min(TABLE_sp$DELTA, na.rm=TRUE))
  y_max <- suppressWarnings(max(TABLE_sp$DELTA, na.rm=TRUE))
  if(!is.finite(y_min)) y_min <- -1000 # Default if no data
  if(!is.finite(y_max)) y_max <- 1000  # Default if no data
  y_min <- y_min - 300
  y_max <- y_max + 300
  
  for(k in 1:4){
    idx <- which(All_sp$Cluster == k)
    
    if (length(idx) == 0) {
      # Draw empty plot if cluster is empty
      plot(1, type="n", axes=TRUE, xlab="Days since capture", ylab=expression(paste(Delta,"(m)")),
           ylim=c(y_min, y_max), main = paste(species_name, "- Cluster", k, "(Empty)"))
      mtext(species_name, side=3, line=1, adj=0, cex=1.3, font=3)
      text(x=10, y=y_max, labels=paste("Cluster",k, "n = 0"), adj=c(1,1), cex=1.6, font=2)
      next
    }
    
    # Adjust cluster size label based on actual data
    current_cluster_size_label <- paste0("Cl ", k, " n = ", length(idx)) # Generic label
    
    n_iter <- out_sp$BUGSoutput$n.keep # Number of samples per chain * chains = n.keep already
    ClusterEst <- matrix(NA, nrow=n_iter, ncol=10)
    kones <- numeric(n_iter)
    
    for(i in 1:n_iter){
      # Extract parameters for fish in cluster 'k' for MCMC sample 'i'
      # Need to handle case where idx might point outside bounds if All_sp was filtered (e.g., Carp 10)
      # We need to map the row index 'i' to the actual fish IDs included in the MCMC output
      
      # Use the BUGS index from All_sp for fish in this cluster
      bugs_indices_in_cluster <- All_sp$ID[idx] 
      
      # Filter bugs_indices to only those present in the MCMC output columns
      valid_bugs_indices <- bugs_indices_in_cluster[bugs_indices_in_cluster <= ncol(out_sp$BUGSoutput$sims.list$kone.i)]
      
      if(length(valid_bugs_indices) == 0) next # Skip if no valid fish found (shouldn't happen here)
      
      # Calculate mean parameters across fish in the cluster for this iteration
      kone  <- mean(out_sp$BUGSoutput$sims.list$kone.i[i, valid_bugs_indices], na.rm=T)
      alpha <- mean(out_sp$BUGSoutput$sims.list$alpha.i[i, valid_bugs_indices], na.rm=T)
      beta  <- mean(out_sp$BUGSoutput$sims.list$beta.i[i, valid_bugs_indices], na.rm=T)
      rr    <- mean(out_sp$BUGSoutput$sims.list$rr.i[i, valid_bugs_indices], na.rm=T)
      ktwo  <- mean(out_sp$BUGSoutput$sims.list$ktwo.i[i, valid_bugs_indices], na.rm=T)
      
      if(any(is.na(c(kone, alpha, beta, rr, ktwo)))) next # Skip if mean params are NA
      
      delta_clust <- numeric(10)
      for(j in 1:10){
        if(j==1){ delta_clust[j] <- kone + rr*kone*(1 - kone/ktwo) - (alpha/(1+beta*j))*kone
        } else { delta_clust[j] <- delta_clust[j-1] + rr*delta_clust[j-1]*(1 - delta_clust[j-1]/ktwo) - (alpha/(1+beta*j))*delta_clust[j-1] }
        if(!is.finite(delta_clust[j])) delta_clust[j] <- NA # Handle potential Inf/-Inf
      }
      ClusterEst[i,] <- delta_clust
      kones[i] <- kone
    }
    
    # Calculate Quantiles (robust to NAs)
    Q50 <- apply(ClusterEst, 2, median, na.rm=TRUE) - Standardize
    Q25 <- apply(ClusterEst, 2, quantile, probs=0.25, na.rm=TRUE) - Standardize
    Q75 <- apply(ClusterEst, 2, quantile, probs=0.75, na.rm=TRUE) - Standardize
    Q05 <- apply(ClusterEst, 2, quantile, probs=0.025, na.rm=TRUE) - Standardize
    Q95 <- apply(ClusterEst, 2, quantile, probs=0.975, na.rm=TRUE) - Standardize
    k50 <- median(kones, na.rm=TRUE) - Standardize
    k25 <- quantile(kones, 0.25, na.rm=TRUE) - Standardize
    k75 <- quantile(kones, 0.75, na.rm=TRUE) - Standardize
    k05 <- quantile(kones, 0.025, na.rm=TRUE) - Standardize
    k95 <- quantile(kones, 0.975, na.rm=TRUE) - Standardize
    
    # Check if quantiles are valid
    if(any(!is.finite(c(Q50, k50)))) {
      plot(1, type="n", axes=TRUE, xlab="Days since capture", ylab=expression(paste(Delta,"(m)")),
           ylim=c(y_min, y_max), main = paste(species_name, "- Cluster", k, "(NA Quantiles)"))
      next
    }
    
    x_vals <- -4:10
    Q50_full <- c(rep(k50,5), Q50)
    Q25_full <- c(rep(k25,5), Q25)
    Q75_full <- c(rep(k75,5), Q75)
    Q05_full <- c(rep(k05,5), Q05)
    Q95_full <- c(rep(k95,5), Q95)
    
    plot(x_vals, Q50_full, type="l", lwd=2, col="black",
         ylim=c(y_min, y_max),
         xlab="Days since capture", ylab=expression(paste(Delta,"(m)")),
         cex.lab=1.7, cex.axis=1.5)
    
    lines(x_vals, Q25_full, col="red", lwd=2)
    lines(x_vals, Q75_full, col="red", lwd=2)
    lines(x_vals, Q05_full, col="blue", lwd=2)
    lines(x_vals, Q95_full, col="blue", lwd=2)
    
    abline(h=k50, lty=3, col="darkgrey", lwd=2)
    abline(v=0, lty=2)
    
    mtext(species_name, side=3, line=1, adj=0, cex=1.3, font=3)
    # Use dynamically calculated N
    text(x=max(x_vals), y=y_max, labels=current_cluster_size_label, adj=c(1,1), cex=1.6, font=2)
  }
}
dev.off() # Close the PNG device


# ---- 6. PLOT 3: LD/MD Cluster Boxplots (ggplot2) ----
# [# NOTE:] Code added from older plotting script
message("Generating Plot 3: LD/MD Cluster Boxplots...")

LDMD_plot_list <- list()
plot_counter <- 1
panel_labels_ldmd <- c("(A)", "(B)", "(C)", "(D)") # Renamed
cluster_labels <- c("low", "mid", "high1", "high2")
species_order <- c("Carp", "Perch")
metrics <- c("LD", "MD")

for (species_name in species_order) {
  All_df <- species_data[[species_name]]$All
  
  if(!all(c("LD", "MD", "Cluster", "RankCluster") %in% names(All_df))) {
    warning(paste("Skipping LD/MD plot for", species_name, "- 'All' missing columns."))
    # Add placeholders to keep grid structure
    LDMD_plot_list[[paste0(species_name, "_LD")]] <- ggplot() + theme_void() + ggtitle(paste(species_name, "LD - Data Missing"))
    LDMD_plot_list[[paste0(species_name, "_MD")]] <- ggplot() + theme_void() + ggtitle(paste(species_name, "MD - Data Missing"))
    plot_counter <- plot_counter + 1 # Increment twice
    next
  }
  
  if (!is.factor(All_df$RankCluster)) {
    All_df$RankCluster <- factor(All_df$Cluster, levels = 1:4, labels = cluster_labels)
  } else { levels(All_df$RankCluster) <- cluster_labels }
  
  for (metric in metrics) {
    p <- ggplot(All_df, aes(x = RankCluster, y = .data[[metric]])) +
      geom_boxplot(fill = ifelse(metric == "LD", "lightblue", "lightgreen"), na.rm = TRUE) +
      labs(x = "Cluster", y = paste(metric, "Response Metric")) +
      theme_bw(base_size = 16) +
      theme( axis.text = element_text(size = 14), axis.title = element_text(size = 15), plot.margin = unit(c(1.2, 1.2, 1.2, 1.2), "cm") ) +
      annotate("text", x = 1, y = max(All_df[[metric]], na.rm = TRUE) * 1.05,
               label = paste0(panel_labels_ldmd[plot_counter], " ", species_name),
               hjust = 0, vjust = 1, size = 6, fontface = "bold")
    
    # Use unique names for the list elements
    LDMD_plot_list[[paste0(species_name, "_", metric)]] <- p
    plot_counter <- plot_counter + 1
  }
}

# Arrange and save the 2x2 panel plot
if(length(LDMD_plot_list) >= 4) { # Check if enough plots were generated
  # Ensure correct order for grid.arrange
  plot_order <- c("Carp_LD", "Perch_LD", "Carp_MD", "Perch_MD")
  valid_plots <- LDMD_plot_list[plot_order[plot_order %in% names(LDMD_plot_list)]]
  
  if(length(valid_plots) == 4) {
    out_path_ldmd <- file.path(OUTPUT_FOLDER, "LDMD_Clusters_Combined.png")
    message(paste("Saving LD/MD cluster plot to:", out_path_ldmd))
    png(out_path_ldmd, width = 40, height = 30, units = "cm", res = 300)
    grid.arrange(grobs = valid_plots, nrow = 2, ncol = 2)
    dev.off()
  } else {
    message("Skipped saving LD/MD plot - not all panels could be generated.")
  }
} else {
  message("Skipped saving LD/MD plot due to errors or missing data.")
}


message("\n✅ Script 02 finished successfully. All plots saved to /Figures.")