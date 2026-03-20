# Restore4Cs_GHG_chamber_fluxes
R workflow to compute CO₂ and CH₄ fluxes from static chamber incubations (RESTORE4Cs project)

This repository contains R scripts used to calculate **instantaneous CO₂ and CH₄ fluxes from static chamber incubations** conducted within the **RESTORE4Cs project** (https://www.restore4cs.eu/).  
The workflow computes flux estimates from greenhouse gas concentration time-series using three different models and automatically selects the **best model for each incubation**.
The repository implements a workflow separated into two scripts. 

# `Fluxes_from_RData.R`

This script produces **three instantaneous flux estimates** for each static chamber CO₂ and CH₄ incubation.  

The script:

- Loads raw GHG concentration time-series data and required metadata for each RESTORE4Cs incubation (available at Zenodo)
- Calculates **LM** and **HM** fluxes using the `goFlux` package
- Calculates an additional **total.flux** estimate (robust to ebullition), based on net concentration change using the first and last **10-second windows**
- Optionally produces diagnostic plots of each incubation
- Outputs tables containing all calculated flux estimates and associated model performance metrics


# `Best_model_selection.R`

This script selects the **best flux estimate** for each CO₂ and CH₄ time-series.
Model selection follows a sequential set of criteria to determine whether **LM**, **HM**, or **total.flux** provides the most appropriate estimate.

1. Incubations flagged due to artifacts or manipulations that preclude flux calculation are excluded.

2. For **CH₄ time-series with ebullitive dynamics** (visually assessed), the **total.flux** estimate is selected instead of **LM**, unless **LM.r² ≥ 0.99**.  
   The HM model is never considered for ebullitive time-series.

3. For the remaining time-series (CO₂ and CH₄ without ebullition), the best model is selected between **LM** and **HM**.

   The non-linear **HM** model is selected only when **all** the following conditions are met (otherwise **LM** is selected):

   - HM produces a valid flux estimate (`HM.flux ≠ NA`)
   - The LM flux estimate is above the minimal detectable flux (`|LM.flux| ≥ MDF.lim`)
   - The HM curvature parameter **κ** is below the theoretical maximum (`HM.k < k.max`)
   - The ratio between HM and LM flux estimates (**g-factor**) is below the gas-specific threshold:  
     - CO₂ g-factor threshold = 4 (`g.fact < 4`)
     - CH₄ g-factor threshold = 3 (`g.fact < 3`)
   - The Akaike Information Criterion corrected for small sample size (**AICc**) of the HM model is lower than or equal to that of the LM model (`HM.AICc ≤ LM.AICc`)
   - The mean absolute error (**MAE**) of the HM model is at least **5% lower** than that of the LM model (`HM.MAE ≤ 0.95 × LM.MAE`)



# Packages

The scripts rely on the following R packages: `goFlux`, `dplyr`, `tidyr`, `ggplot2`, `purrr`, and `stringr`.
Package versions are managed using **`renv`**, ensuring a reproducible computational environment.



# Input

Input data are stored in Zenodo (https://doi.org/10.5281/zenodo.18803756). It consits of: 

- **RData files** containing raw CO₂ and CH₄ concentration time-series
- **co2_auxfile.csv** and **ch4_auxfile.csv**, containing metadata required for flux calculation
- **ebullition_inspection.csv**, recording visually identified ebullitive dynamics in each CH₄ time-series



# Output

The workflow produces:

- Tables containing **all calculated flux estimates**
- Optional **diagnostic plots** of incubation time-series and model fits
- Final tables containing **best-flux selections for CO₂ and CH₄**

Best-flux datasets are also archived in the **LifeWatch ERIC repository**:  
https://doi.org/10.48372/C29B-QW38



# Authors

Miguel Cabrera-Brufau  
Camille Minaudo  

**RESTORE4Cs Project**
