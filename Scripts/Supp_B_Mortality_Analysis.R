################################################################################
# SCRIPT Supp_B: Mortality Determination Process and Rate Comparison (Documentation)
#
# DESCRIPTION:
# This supplementary script documents the process used to determine fish
# mortality status and performs the statistical comparison of mortality rates
# between impact (caught) and control groups using the final Mort_All dataset.
# Mortality status was primarily determined through visual inspection of daily
# fish tracks plotted against the lake shoreline, aided by automated flagging
# of low-movement days. GLMs predicting mortality based on metadata were NOT
# performed as part of the final analysis presented.
#
# INPUTS:
# - Input/CarpPos.rdata, Input/PerchPos.rdata (Raw positions for demo)
# - Input/CarpCatch.rdata, Input/PerchCatch.rdata (Needed for context)
# - Input/uferlinie.shp (and associated files, for plotting demo)
# - Results/Mort_All.rdata (Contains FINAL mortality status for statistics)
#
# OUTPUTS:
# - Console output (statistical summaries).
# - Plots (if plotting examples are uncommented and run manually).
#
################################################################################

# --- 1. SETUP ---
# setwd("D:/Daten/Paper/Repository") # Set repo root if running standalone
library(fasttime)
library(dplyr)
library(ggplot2) # Although not used in this version, kept for consistency
library(rgdal)   # For reading shapefiles (legacy)
library(sf)      # For reading shapefiles (modern)

message("Starting Supplementary Mortality Analysis Script (Supp_B)")

# --- 2. LOAD CORE DATA ---
message("Loading Position, Catch, and Shoreline data...")
# [# NOTE:] Load necessary raw data files from the Input folder.
tryCatch({ load("Input/CarpPos.rdata", .GlobalEnv) }, error = function(e){ warning("CarpPos.rdata not found.") })
tryCatch({ load("Input/PerchPos.rdata", .GlobalEnv) }, error = function(e){ warning("PerchPos.rdata not found.") })
tryCatch({ load("Input/CarpCatch.rdata", .GlobalEnv) }, error = function(e){ warning("CarpCatch.rdata not found.") }) # Loads CarpCatch2
tryCatch({ load("Input/PerchCatch.rdata", .GlobalEnv) }, error = function(e){ warning("PerchCatch.rdata not found.") }) # Loads PerchCatch2

# Load shoreline shapefile robustly
shoreline <- NULL
tryCatch({
  shoreline <- rgdal::readOGR(dsn = "Input", layer = "uferlinie", verbose = FALSE)
  message("Shoreline loaded using rgdal.")
}, error = function(e_rgdal) {
  tryCatch({
    shoreline <- sf::st_read(dsn = "Input/uferlinie.shp", quiet = TRUE)
    message("Shoreline loaded using sf.")
  }, error = function(e_sf) {
    warning("Could not load shoreline shapefile (Input/uferlinie.shp). Plotting functions disabled.")
  })
})

# --- 3. DEMONSTRATION OF MORTALITY DETERMINATION PROCESS ---
# [# NOTE:] This section recalculates the 'alive' stats from raw position data
# to demonstrate the automated flagging method used prior to visual inspection.
# It uses Carp data as an example.

message("\nDemonstrating mortality determination process (using Carp data)...")

if (exists("CarpPos")) {
  # Select Carp data and preprocess time
  DATA_demo <- CarpPos
  DATA_demo$time <- as.POSIXct(paste(DATA_demo$date, DATA_demo$time), format = "%d.%m.%Y %H:%M:%S", tz = "GMT")
  DATA_demo$days <- fastPOSIXct(cut(DATA_demo$time, "days"), tz = "GMT")
  IDlist_demo <- sort(unique(DATA_demo$ID))
  
  # 3.1 Calculate Automated Flags: Identify potential inactivity periods
  message("... Calculating daily movement stats for automated flagging demonstration...")
  alive_stats_demo <- DATA_demo %>%
    group_by(ID, days) %>%
    summarise(
      sdx = sd(x, na.rm = TRUE),
      sdy = sd(y, na.rm = TRUE),
      # Calculate depth range robustly
      rangez = {
        depth_vals <- cur_data()[["sensor.depth"]] # Access column safely
        if(!is.null(depth_vals) && sum(!is.na(depth_vals)) > 1) {
          max(depth_vals, na.rm = TRUE) - min(depth_vals, na.rm = TRUE)
        } else { NA_real_ }
      },
      n = n(),
      .groups = 'drop'
    ) %>%
    # Add IDnum and daynum for easier referencing
    mutate(
      IDnum = match(ID, IDlist_demo),
      daynum = match(days, sort(unique(DATA_demo$days))),
      # Flag days with minimal depth range or missing depth data
      maybedead = case_when(
        is.na(rangez) ~ 1,         # Flag if depth data is missing
        rangez < 0.45 ~ 1,         # Flag if range is very small
        TRUE ~ 0                   # Otherwise, assume alive
      )
    )
  message("... Example 'alive' stats calculation complete.")
  # Example: print(head(alive_stats_demo))
  
  # 3.2 Visual Inspection Process: Plotting function definition
  message("... Defining plotting function used for visual inspection...")
  plot_pos_demo <- function(IDnum, daynum, data_source = DATA_demo, shoreline_obj = shoreline,
                            alive_df = alive_stats_demo, id_list = IDlist_demo) {
    # [# NOTE:] This function plots daily fish positions against the shoreline.
    # It was used iteratively to visually assess movement patterns, especially
    # for fish/days flagged by 'maybedead'.
    if (is.null(shoreline_obj)) { message("Shoreline not loaded, skipping plot."); return() }
    
    # Get ID and Date
    ID <- id_list[IDnum]
    if (length(ID) == 0 || is.na(ID)) { print(paste("Invalid IDnum", IDnum)); return() }
    day_data <- alive_df %>% filter(IDnum == !!IDnum) %>% arrange(daynum)
    if (daynum > nrow(day_data) || daynum <= 0) { print(paste("daynum", daynum, "out of range.")); return()}
    day <- day_data$days[daynum]
    day_str <- format(day, format = "%Y-%m-%d")
    print(paste("Plotting:", ID, "(IDnum", IDnum, ") on", day_str, "(daynum", daynum, ")"))
    
    # Plot shoreline
    plot_args <- list(lwd = 1, main = "")
    if (inherits(shoreline_obj, "Spatial")) { do.call(plot, c(list(x = shoreline_obj), plot_args)) }
    else if (inherits(shoreline_obj, "sf")) { do.call(plot, c(list(x = sf::st_geometry(shoreline_obj)), plot_args)) }
    
    # Add fish positions
    fish_day_data <- data_source[which(data_source$ID == ID & data_source$days == day), ]
    if(nrow(fish_day_data) > 0) {
      points(fish_day_data$x, fish_day_data$y, pch = 16, cex = 0.5)
    }
    
    # Add label
    usr_coords <- par("usr")
    text(usr_coords[1] + (usr_coords[2]-usr_coords[1])*0.05,
         usr_coords[4] - (usr_coords[4]-usr_coords[3])*0.1,
         paste(IDnum, daynum, day_str), col = "red", font = 2, adj=0)
  }
  # Example Manual Call:
  # plot_pos_demo(IDnum = 4, daynum = 10) # Plots day 10 for fish IDnum 4
  
  # [# NOTE:] The final mortality status (0 or 1) recorded in 'Mort_All.rdata'
  # was determined by visually inspecting plots like these, especially for days
  # flagged by 'maybedead'. This involved checking if the fish remained stationary
  # over consecutive days, potentially near the shore. Lists like 'mortalityC_m2'
  # in older scripts represent the output of this manual process.
  
} else {
  message("... Skipping demonstration: CarpPos.rdata not loaded.")
}


# --- 4. STATISTICAL COMPARISON OF MORTALITY RATES ---
message("\nCalculating Mortality Rate Comparison using Mort_All.rdata...")

# Load the final dataset containing definitive mortality status
mort_all_path <- "Results/Mort_All.rdata"
if (!file.exists(mort_all_path)) {
  stop("FATAL: Mort_All.rdata not found in Results/. Cannot perform statistical comparison.")
}
load(mort_all_path, .GlobalEnv) # Loads Mort_All object

# Initialize results table
mortality_results <- data.frame(
  Species = character(),
  N_Impact = integer(), N_Impact_Dead = integer(),
  N_Control = integer(), N_Control_Dead = integer(),
  MT_percent = numeric(), LowerCI_percent = numeric(), UpperCI_percent = numeric()
)
species_names_map <- c("Perch", "Carp") # Assumes Species column uses 1 and 2

for (species_code in c(1, 2)) {
  species_name <- species_names_map[species_code]
  
  # Check if Mort_All contains data for this species
  if(!species_code %in% Mort_All$Species) {
    message(paste("Skipping", species_name, "- No data found in Mort_All."))
    next
  }
  data_subset <- Mort_All[Mort_All$Species == species_code, ]
  
  # Calculate counts directly from Mort_All
  # Group=1 is Impact, Group=0 is Control
  # Mortality=1 is Dead, Mortality=0 is Alive
  nf  <- sum(data_subset$group == 1, na.rm=TRUE) # N Impact
  nfd <- sum(data_subset$group == 1 & data_subset$Mortality == 1, na.rm=TRUE) # N Impact Dead
  nc  <- sum(data_subset$group == 0, na.rm=TRUE) # N Control
  ncd <- sum(data_subset$group == 0 & data_subset$Mortality == 1, na.rm=TRUE) # N Control Dead
  
  # Proceed only if both groups have fish
  if(nf == 0 || nc == 0) {
    message(paste("Skipping", species_name, "- zero fish in impact or control group."))
    next
  }
  
  # Calculate rates and difference
  MR <- nfd / nf
  MC <- ncd / nc
  MT <- MR - MC
  
  # Calculate variance and confidence interval (Wilde et al. 2004)
  varMC <- (MC * (1 - MC)) / nc
  varMR <- (MR * (1 - MR)) / nf
  varMC[is.nan(varMC)] <- 0 # Handle cases where rates are 0 or 1
  varMR[is.nan(varMR)] <- 0
  varMT <- varMC + varMR
  
  # Use Z=1.96 for 95% CI based on normal approximation
  CI_half_width <- 1.96 * sqrt(varMT)
  
  # Convert to percentages
  upper_ci = (MT + CI_half_width) * 100
  lower_ci = (MT - CI_half_width) * 100
  mt_percent = MT * 100
  
  # Store results
  result_row <- data.frame(
    Species = species_name,
    N_Impact = nf, N_Impact_Dead = nfd,
    N_Control = nc, N_Control_Dead = ncd,
    MT_percent = mt_percent,
    LowerCI_percent = lower_ci,
    UpperCI_percent = upper_ci
  )
  mortality_results <- rbind(mortality_results, result_row)
} # End species loop

print("--- Mortality Rate Difference (Impact - Control) based on Mort_All.rdata ---")
print(mortality_results)

# [# NOTE:] The calculation follows Wilde et al. (2004) N Am J Fish Manage 23:779-786
# for estimating the variance of the difference between two proportions.

message("\n✅ Supp_B script finished.")
# Clean up demonstration objects and loaded data (optional)
rm(list=ls(pattern="_demo$|^shoreline$|Mort_All$|CarpPos$|PerchPos$|CarpCatch2$|PerchCatch2$"))