# Restore4Cs_GHG_chamber_fluxes
This repository contains R scripts used to calculate **instantaneous CO₂ and CH₄ fluxes from static chamber incubations** conducted within the **RESTORE4Cs project**. 
The workflow computes multiple flux estimates from greenhouse gas concentration time-series and selects the **most appropriate flux model for each incubation**. The repository implements a workflow separated into two scripts: 

## `Fluxes_from_RData.R`

This scripts produces three instantaneous flux estimates for each static chamber CO₂ and CH₄ incubation conducted. The script:

- Loads raw GHG concentration time-series data and required metadata for every RESTORE4Cs incubation (available at Zenodo)
- Calculates **LM** and **HM** fluxes using the `goFlux` package 
- Calculates an additional **total.flux** estimate (robust to ebullition), based on net concentration change using initial and final 10s-windows
- Optionally produces diagnostic plots of each incubation
- Outputs tables containing all calculated flux estimates and associated model performance metrics.

## `Best_model_selection.R`

This script selects the **best flux estimate** for each CO₂ and CH₄ time-series. Model selection follows a sequential set of criteria to determine whether **LM**, **HM**, or **total.flux** provides the most appropriate estimate:
  1. Incubations flagged due to artifacts or manipulations precluding flux calculation are excluded.
  2. For CH4 time-series with ebullitive dynamics (visually assessed), the **total.flux** is chosen over **LM** unless LM.r2 is above 0.99 (HM is never considered for ebullitive time-series).
  3. For the rest of time-series (CO₂ and CH₄, non-ebullitive), the best-model is chosen beween **LM** and **HM**. The non-linear **HM** model is chosen only when all the following criteria are met (defaulting to the **LM** as best model when one or more are violated):
     - HM model produces a valid flux estimate (HM.flux ≠ NA)
     - LM flux estimate is above the minimal detectable flux (abs(LM.flux) >= MDF.lim)
     - HM curvature parameter Kappa is below the theoretical maximum (HM.k < k.max)
     - The ratio between the non-linear (HM) flux estimate and the linear (LM) estimate, the g-fact is below the gas species-specific threshold (CO₂ g-fact < 4; CH₄ g-fact < 3)
     - Akaike Information Criterion corrected for small sample size (AICc) of the HM model is lower or equal than that of the LM model (HM.AICc <= LM.AICc)
     - Mean absolute error (MAE) of HM model is at least 5% lower than that of LM model (HM.MAE <= 0.95 * LM.MAE)

# Packages

The scripts are rely on the following packages: `goFlux`, `dplyr`, `tidyr`, `ggplot2`, `purrr` and `stringr`.
Package versions are managed using **`renv`**, ensuring a reproducible environment.

# Input
Input data is stored at Zenodo (https://doi.org/10.5281/zenodo.18803756), where they are described in detail. Briefly, the input data consist on:
-  RData files: contain raw CO₂ and CH₄ concentration time-series
-  co2_auxfile.csv and ch4_auxfile.csv: contain details required to calculate instantaneous GHG flux of each incubation
-  ebullition_inspection.csv: contains the recording of vissually apparent ebullitive dynamics in each CH₄ time series

# Output

The workflow produces:

- Tables containing **all calculated flux estimates**
- Optionally, **diagnostic plots** of incubation time-series and model fits
- Final tables containing **best-flux selections for CO₂ and CH₄** (stored also at LifeWatch ERIC repo: https://doi.org/10.48372/C29B-QW38)

# Authors

Miguel Cabrera-Brufau  
Camille Minaudo  

**RESTORE4Cs Project**
