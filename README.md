# Accounting for Human Movement to Improve Exposure-Health Models

**Hadiqa Tahir, Simon Smart, Samuel Cai, André Ng, Joshua Vande Hey, Tim CD Lucas**

School of Medical Sciences, University of Leicester, Leicester, UK

> Preprint forthcoming on medRxiv. This repository will be updated with a DOI when available.

------------------------------------------------------------------------

## Overview

This repository contains the code for the simulation study described in the paper. The study develops and evaluates statistical models that explicitly account for human movement when linking disaggregated environmental exposures to individual-level health outcomes.

Three models are implemented and compared:

- **Fixed-Weight (FW)**: inverse-distance weights with α = 1 fixed a priori
- **Learned-Weight (LW)**: extends FW by estimating the distance-decay exponent α from the data
- **Restricted Cubic Spline (RCS)**: applies a nonlinear (spline) transformation to NDVI before passing to the FW model

All three share the same core structure:

```         
logit P(ARI_j) = Σ_i [ w_i × f(NDVI_i) ]
```

where `w_i = (1/d_i^α) / Σ(1/d_i^α)` and `f(·)` is linear (FW/LW) or a restricted cubic spline (RCS).

Model performance was assessed across three sample sizes (N = 1,114; 50,000; 100,000) using bias, empirical standard error (EmpSE), mean square error (MSE), and credible interval coverage, each with Monte Carlo errors.

------------------------------------------------------------------------

## Data

### Real DHS data (not included)

The paper uses the Albania 2017–18 Demographic and Health Surveys (DHS) dataset, linking acute respiratory infection (ARI) outcomes in children under five to pixel-level NDVI within a 3 km buffer around DHS cluster centroids.

**The original DHS data cannot be shared publicly.** DHS data are subject to data use agreements and require an approved application. To reproduce the paper's exact results:

1.  Apply for access at <https://dhsprogram.com>
2.  Request the Albania 2017–18 Individual Recode (IR) and Geographic datasets
3.  Replace `data/df_1114_synthetic.parquet` with the processed dataset

The pipeline expects a Parquet file with the following columns:

| Column | Description |
|------------------------------------|------------------------------------|
| `HOUSE` | Individual ID (integer, one per child) |
| `ARI` | Binary outcome: 1 = ARI, 0 = no ARI |
| `NDVI` | NDVI value for a pixel within the individual's 3 km buffer |
| `DISTANCES` | Euclidean distance (metres) from pixel centroid to cluster centroid |
| `WEIGHT` | Inverse-distance weight for this pixel: `(1/d) / Σ(1/d)` |
| `X_cluster` | Easting of cluster centroid (UTM Zone 34N, EPSG:32634) |
| `Y_cluster` | Northing of cluster centroid (UTM Zone 34N, EPSG:32634) |
| `x_ndvi` | Easting of NDVI pixel centroid |
| `y_ndvi` | Northing of NDVI pixel centroid |

In this data, '`HOUSE'` is the same as 'individual' or 'child'. The data are in long format: one row per pixel per individual.

### Synthetic data (included)

Because the real data cannot be shared, a synthetic dataset is provided that mimics the structure of the original data without containing any real survey records.

To regenerate the synthetic dataset, run:

``` r
source("data/generate_synthetic_data.R")
```

This requires the Albania NDVI raster (`data/albania_ndvi.tif`) and boundary shapefile (`data/df_albania_boundary_shp/`), which are included in the repository.

In order to test the sample size, this observed DHS data were sampled with replacement to sample sizes: 50,000 and 100,000. These data are too large to provide, therefore a script is available to generate resampled data.

``` r
source("data/generate_resampled_data.R")
```

------------------------------------------------------------------------

## Repository Structure

```         
.
├── data/
│   ├── generate_synthetic_data.R          # Generates the synthetic dataset
│   ├── generate_resampled_data.R          # Generates the resampled datasets
│   └── df_1114_synthetic.parquet  # Synthetic dataset (N = 1,114)
│   └── albania_ndvi.tif                   # NDVI Raster file
│   └── df_albania_boundary.shp            # folder containing Albania's shapefile.
│
├── fixed_weight_model/
│   ├── Fixed_weight_model.cpp             # TMB C++ model (shared with RCS)
│   ├── FW_Source_Functions.R              # Helper functions
│   └── FW_Model_Simulation_Study.R        # Main simulation study script
│
├── rcs_model/
│   ├── Fixed_weight_model.cpp             # TMB C++ model (shared with RCS)
│   ├── RCS_Source_Functions.R             # Helper functions
│   └── RCS_Model_Simulation_Study.R       # Main simulation study script
│
└── learned_weight_model/
    ├── Learned_weights_model.cpp          # TMB C++ model
    ├── LW_Source_Functions.R              # Helper functions
    └── LW_Model_Simulation_Study.R        # Main simulation study script
```

The RCS model reuses `Fixed_weight_model.cpp` — the spline transformation is applied to the exposure matrix in R before passing it to TMB.

------------------------------------------------------------------------

## Dependencies

### R packages

``` r
install.packages(c(
  "dplyr", "ggplot2", "TMB", "data.table",
  "arrow", "readr", "here", "gridExtra",
  "splines", "terra", "sf"
))
```

### TMB

TMB (Template Model Builder) interfaces C++ models with R. Install from CRAN:

``` r
install.packages("TMB")
```

A working C++ compiler is required.

------------------------------------------------------------------------

## How to Run

Each model has a single self-contained script that runs all six steps in sequence: fit on observed data → simulate outcomes → refit on simulations → extract results → compute performance measures → produce plots.

Before running, open the script and set:

``` r
sample_size <- '1114'   # or '50K' or '100K' (resampled datasets not provided due to memory, script available to generate.)
nsims       <- 10       # increase for the full simulation study
```

Then run the scripts to conduct the simulation study:

``` r
source("fixed_weight_model/FW_Model_Simulation_Study.R")
source("rcs_model/RCS_Model_Simulation_Study.R")
source("learned_weight_model/LW_Model_Simulation_Study.R")
```

Outputs are saved to `results_<sample_size>/` within each model folder, including RDS files for estimates, covariance matrices, and performance measures, and PDF diagnostic plots.

------------------------------------------------------------------------

## Correspondence

Hadiqa Tahir — [ht233\@leicester.ac.uk](mailto:ht233@leicester.ac.uk) Tim CD Lucas — [tim.lucas\@leicester.ac.uk](mailto:tim.lucas@leicester.ac.uk)
