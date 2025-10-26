# Fish Movement Analysis (Carp & Perch) - BACIPS Modeling

**Repository for the manuscript cjfas-2025-0348: "Hooking mortality, behavioural impact and recovery from catch-and-release angling in a piscivorous (Eurasian perch, Perca fluviatilis) and a non-piscivorous fish (common carp, Cyprasinus carpio) assessed using high-resolution acoustic telemetry in a shallow natural lake" submitted to Aquatic Sciences.**

This repository contains the data and R scripts necessary to reproduce the final figures for the study using pre-computed model results. It also includes scripts documenting the full analysis pipeline used to generate those results.

---

## Important Note on Large Files

The raw JAGS output files (`Results/outC.rdata` - ~112MB, `Results/outP.rdata` - ~197MB) are required to run `Scripts/02_Generate_Plots.R` but exceed GitHub's file size limits.

**Please download both files from the Google Drive folder linked below and place them in the `Results/` folder before running the plotting script:**

* **Download Folder:** [https://drive.google.com/drive/folders/1Q7OVPNtXZy4qSJrn0ExPLqP9TjlQPX_9?usp=drive_link]

These `out*` files will be included in the final data archive (e.g., Dryad) upon publication. The very large raw position files (`Input/*Pos.rdata` - >400MB each) are also omitted from this repository for size reasons and will be included in the final archive.

---

## How to Reproduce Figures (Recommended)

The primary script for reproducing the manuscript figures uses pre-computed results and runs quickly. Please ensure your R working directory is set to the repository root (`D:/Daten/Paper/Repository`).

1.  **Download `outC.rdata` and `outP.rdata`** using the link above and place them in the `/Results` folder.
2.  **Run `Scripts/02_Generate_Plots.R`**:
    * This script loads the **final, pre-computed `All*.rdata` files** (containing merged parameters and metadata) and the downloaded raw JAGS output (`out*.rdata`).
    * It generates all manuscript figures and saves them into the `/Figures` directory (which will be created if it doesn't exist).

---

## Repository Structure

* `/Input`: Contains supplementary data such as catch records (`*Catch.rdata`), behavior metrics (`*behaviors2.rdata`), raw position data (`*Pos.rdata` - used only for documentation scripts), and the lake shoreline file (`uferlinie.*`).
* `/Results`: Contains model-ready delta tables (`TAB_DEL_*.rdata`), the pre-computed final model results including metadata (`All*.rdata`), raw JAGS summaries (`All_raw_*`), mortality analysis data (`Mort*`), and cluster assignments (`Clustnames.rdata`). **Note:** Large raw JAGS output (`out*.rdata`) must be downloaded separately (see link above). Predicted trajectory files (`delta.3.pred.*`) are also not included but can be generated using `Scripts/Delta_Predictions.R`.
* `/Scripts`: Contains the main plotting script (`02_Generate_Plots.R`), the core model execution script (`01_Run_BACIPS_Model.R`), the model definition (`BACIPS.R`), and supplementary documentation scripts (`Supp_*`, `Delta_Predictions.R`).
* `/Figures`: Destination for all generated plots (created automatically by `02_Generate_Plots.R`).

---

## Full Workflow Documentation (Optional, Time-Intensive)

The following scripts document the complete workflow used to generate the pre-computed results found in the `/Results` folder. They are provided for transparency and full reproducibility but are not required just to generate the final figures.

1.  **`Scripts/01_Run_BACIPS_Model.R` (Core Method):**
    * Executes the main Bayesian analysis using `BACIPS.R`, taking `Results/TAB_DEL_*.rdata` as input.
    * Saves the raw model outputs (`Results/out*.rdata`, `Results/All_raw_*.rdata`).
    * **Warning:** Running this script takes **several hours**. ⏳

2.  **(Process documented in `Scripts/Supp_A_...`)**:
    * Details the post-processing steps: Random Forest clustering, LD/MD metric calculation, merging parameters with metadata (from `Input/*Catch.rdata`), and running meta-models.
    * This process uses `Results/All_raw_*.rdata` and `Results/delta.3.pred.*` (generated via `Delta_Predictions.R`) as input and creates the final `Results/All*.rdata` and `Results/Clustnames.rdata`.

---

## Supplementary Scripts (Documentation Only)

These scripts provide details on specific analysis steps:

* **`Scripts/Supp_A_Cluster_Analysis_and_Data_Merging.R`**: Shows the Random Forest clustering, LD/MD calculation, merging process to create the final `All*.rdata`, and meta-model analysis. **Note:** Requires `delta.3.pred.*` files, which must be generated first by running `Scripts/Delta_Predictions.R` (takes significant time).
* **`Scripts/Supp_B_Mortality_Analysis.R`**: Documents the process for determining mortality status (flagging using position data, visual inspection examples) and compares mortality rates between groups using `Results/Mort_All.rdata`.
* **`Scripts/Delta_Predictions.R`**: Calculates the full predicted delta trajectories (`delta.3.pred.*` files) from raw JAGS output (`Results/out*.rdata`). **Warning:** Takes significant time to run.
* **`Scripts/BACIPS.R`**: Defines the JAGS model structure used by `01_Run_BACIPS_Model.R`.