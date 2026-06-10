# =============================================================================#
# LW_Model_Simulation_Study.R
# Author:  Hadiqa Tahir
# Date:    June 2026
# =============================================================================#
# Research Questions:
#       - Can the model recover parameter estimates?
#       - How does sample size affect model performance?
#
# Learned-Weight (LW) model details:
#   logit P(ARI_j) = sum_i [ w_i * (b0 + b1 * NDVI_i) ]
#   where w_i = (1/d_i^alpha) / sum(1/d_i^alpha)  (alpha is estimated, not fixed)
#   Parameters: b0, b1, log_alpha
#   Priors: b0, b1 ~ N(0, 10^2);  log_alpha ~ N(-1, 1^2)
#
#   Unlike the FW model, raw distances are passed to TMB and weights are
#   computed internally using the learned alpha parameter.
#
# This file covers:
#   1. Fitting the LW model to 'observed' data to obtain "true" parameters.
#   2. Simulating ARI outcomes from the fitted model (nsims replicates).
#   3. Refitting the model on each simulated dataset.
#   4. Extracting MAP estimates, SEs, and covariance matrices.
#   5. Computing performance measures with Monte Carlo errors.
#   6. Producing diagnostic plots (zipper, histogram, scatter, SE distribution).
#
# Dependencies:
#   - LW_Source_Functions.R                : helper functions
#   - Learned_weights_model.cpp            : TMB C++ model
#   - df_<sample_size>_synthetic.parquet   : synthetic dataset mimicking real DHS data
#
# Outputs (written to learned_weight_model/results_<sample_size>/):
#   - results_<sample_size>.rds               : per-simulation MAP estimates and SEs
#   - covariance_matrix_<sample_size>.rds     : per-simulation covariance matrices
#   - performance_measures_<sample_size>.rds
#   - plots/                                  : PDF diagnostic plots
#
# Real DHS data cannot be shared publicly. To reproduce the paper's results,
# apply for the Albania 2017-18 DHS data at https://dhsprogram.com
# =============================================================================#

library(dplyr)
library(ggplot2)
library(TMB)
library(data.table)
library(arrow)
library(readr)
library(here)
library(gridExtra)

# Set working directory
setwd(here::here("learned_weight_model"))

# Functions
source('LW_Source_Functions.R')

#==============================================================================#
# 1. Run model on observed data                                               ----
#==============================================================================#

sample_size <- '1114'         # Change to '50K' or '100K' for resampled datasets.
# Resampled data not provided on github due to large memory.

# Load the dataset
data_path <- here::here(paste0("df_", sample_size, "_synthetic.parquet"))
df <- open_dataset(data_path) %>% collect()

# Outcome data as a vector
ARI <- df %>%
  select(HOUSE, ARI) %>%
  group_by(HOUSE) %>%
  summarise(ARI = first(ARI)) %>%
  collect() %>%
  pull(ARI)

# Load in the indices
n_pixels_vec <- df %>% count(HOUSE) %>% pull(n)      # how many pixels per person?
startindex   <- cumulative_sum(n_pixels_vec)          # start/end index for each pixel

# Exposure data as a matrix
ndvi_and_cov <- as.matrix(df$NDVI)

# Raw distances — passed directly to TMB; weights are computed inside the model
# using the learned alpha: w_i = (1/d_i^alpha) / sum(1/d_i^alpha)
distances <- df$DISTANCES

# Priors for log_alpha
priormean_log_alpha <- -1
priorsd_log_alpha   <- 1

# Compile and load model
compile("Learned_weights_model.cpp", "&> logfile.log")
dyn.load(dynlib("Learned_weights_model"))

# Run model
results <- TMB_model(ARI, ndvi_and_cov, distances, startindex, n_pixels_vec,
                     priormean_log_alpha, priorsd_log_alpha)
print(results)

# Save the observed model fit
path <- paste0('observed_fit_', sample_size, '.rds')
saveRDS(results, file = path)

rm(ARI, data_path, path)

#==============================================================================#
# 2. Simulate data from observed model fit                                    ----
#==============================================================================#

nsims <- 10
dir.create(paste0("simulated_outcomes_", sample_size), showWarnings = TRUE)
dir.create(paste0("results_", sample_size),            showWarnings = TRUE)

# Extract true parameters from observed fit
true_b0        <- results$par.fixed[['intercept']]
true_b1        <- results$par.fixed[['slope']]
true_log_alpha <- results$par.fixed[['log_alpha']]
alpha          <- exp(true_log_alpha)

# Compute weights using learned alpha, then probabilities
# w_i = (1/d_i^alpha) / sum(1/d_i^alpha)
df$WEIGHT <- (1 / df$DISTANCES^alpha) / ave((1 / df$DISTANCES^alpha), df$HOUSE, FUN = sum)

n_j   <- with(df, stats::ave(WEIGHT * (true_b0 + (true_b1 * NDVI)), HOUSE, FUN = sum, na.rm = TRUE))
p_ari <- 1 / (1 + exp(-n_j))
df_nj_pj <- data.frame(HOUSE = df$HOUSE, WEIGHT = df$WEIGHT, n_j = n_j, p_ari = p_ari)

# Generate the simulated outcomes
output_folder <- paste0("simulated_outcomes_", sample_size)
for (sim_index in 1:nsims) {
  simulation_result           <- simulate_disease_data(df_nj_pj, seed = 4000, i = sim_index)
  simulation_result           <- data.frame(ARI = simulation_result)
  colnames(simulation_result) <- paste0("Sim_", sim_index)
  file_name <- file.path(output_folder, paste0("Sim_", sim_index, ".parquet"))
  write_parquet(simulation_result, file_name)
  message("Simulation ", sim_index, " saved.")
}

rm(df_nj_pj, n_j, p_ari, alpha, output_folder, sim_index, simulation_result, file_name)

#==============================================================================#
# 3. Run model on simulated data                                              ----
#==============================================================================#

# All inputs are the same as section 1; only the outcome vector changes.
simulation_file <- file.path(paste0("simulated_outcomes_", sample_size), paste0("Sim_", 1:nsims, ".parquet"))
results_folder  <- paste0("results_", sample_size)

for (i in 1:nsims) {
  ARI     <- arrow::open_dataset(simulation_file[i]) %>% collect()
  cases   <- as.vector(as.integer(ARI[[1]]))
  results <- TMB_model(cases, ndvi_and_cov, distances, startindex, n_pixels_vec,
                       priormean_log_alpha, priorsd_log_alpha)
  name        <- gsub(".parquet", "_results.rds", basename(simulation_file[i]))
  output_file <- file.path(results_folder, name)
  saveRDS(results, output_file)
  message("Saved ", output_file)
}

rm(ndvi_and_cov, distances, startindex, n_pixels_vec,
   priormean_log_alpha, priorsd_log_alpha,
   simulation_file, ARI, cases, results, name, output_file, i)

#==============================================================================#
# 4. Extract results from simulation model runs                               ----
#==============================================================================#

files          <- list.files(path = results_folder, pattern = "Sim_.*_results.rds", full.names = TRUE)
estimates_data <- data.frame()

for (i in seq_along(files)) {
  r <- readRDS(files[i])
  estimates_data <- rbind(estimates_data, data.frame(
    Simulation    = i,
    b0            = r$par.fixed[['intercept']],
    b1            = r$par.fixed[['slope']],
    log_alpha     = r$par.fixed[['log_alpha']],
    se_b0         = sqrt(diag(solve(r$hess)))[[1]],
    se_b1         = sqrt(diag(solve(r$hess)))[[2]],
    se_log_alpha  = sqrt(diag(solve(r$hess)))[[3]],
    nllprior      = r$nll_priors,
    nllfinal      = r$nll_final
  ))
  rm(r)
}

saveRDS(estimates_data, file.path(results_folder, paste0("results_", sample_size, ".rds")))

# Extract covariance matrices
covariance_matrix <- lapply(files, function(f) readRDS(f)$cov.fixed)
saveRDS(covariance_matrix, file.path(results_folder, paste0("covariance_matrix_", sample_size, ".rds")))

rm(i, cumulative_sum, simulate_disease_data, TMB_model, files)

#==============================================================================#
# 5. Assess Performance Measures: Bias, EmpSE, MSE, Coverage                 ----
#==============================================================================#

# Credible intervals for each parameter
estimates_data$ci_lower_b0        <- ci_lower(estimates_data$b0, estimates_data$se_b0, 0.025)
estimates_data$ci_upper_b0        <- ci_upper(estimates_data$b0, estimates_data$se_b0, 0.025)
estimates_data$ci_lower_b1        <- ci_lower(estimates_data$b1, estimates_data$se_b1, 0.025)
estimates_data$ci_upper_b1        <- ci_upper(estimates_data$b1, estimates_data$se_b1, 0.025)
estimates_data$ci_lower_log_alpha <- ci_lower(estimates_data$log_alpha, estimates_data$se_log_alpha, 0.025)
estimates_data$ci_upper_log_alpha <- ci_upper(estimates_data$log_alpha, estimates_data$se_log_alpha, 0.025)

# Variance of each parameter from covariance matrix
for (i in 1:nsims) {
  estimates_data$b0_var[i]        <- covariance_matrix[[i]]["intercept",  "intercept"]
  estimates_data$b1_var[i]        <- covariance_matrix[[i]]["slope",      "slope"]
  estimates_data$log_alpha_var[i] <- covariance_matrix[[i]]["log_alpha",  "log_alpha"]
}
rm(i, covariance_matrix)

# Contraction and z-scores
# Note: log_alpha has a tighter prior (sd = 1) than b0/b1 (sd = 10)
sd_prior_b0        <- 10
sd_prior_b1        <- 10
sd_prior_log_alpha <- 1

estimates_data <- estimates_data %>%
  mutate(
    contraction_b0        = 1 - (b0_var / sd_prior_b0^2),
    contraction_b1        = 1 - (b1_var / sd_prior_b1^2),
    contraction_log_alpha = 1 - (log_alpha_var / sd_prior_log_alpha^2),
    z_score_b0            = (b0 - true_b0) / sqrt(b0_var),
    z_score_b1            = (b1 - true_b1) / sqrt(b1_var),
    z_score_log_alpha     = (log_alpha - true_log_alpha) / sqrt(log_alpha_var)
  )

rm(sd_prior_b0, sd_prior_b1, sd_prior_log_alpha)

# Bias, Empirical SE, MSE, Coverage
performance_measures <- data.frame(
  Parameter  = c("B0", "B1", "log_alpha"),
  True_Value = c(true_b0, true_b1, true_log_alpha),
  Bias = c(
    Bias(estimates_data$b0, true_b0, nsims),
    Bias(estimates_data$b1, true_b1, nsims),
    Bias(estimates_data$log_alpha, true_log_alpha, nsims)
  ),
  MCE_Bias = c(
    MCE_bias(nsims, estimates_data$b0),
    MCE_bias(nsims, estimates_data$b1),
    MCE_bias(nsims, estimates_data$log_alpha)
  ),
  Empirical_SE = c(
    EmpiricalSE(nsims, estimates_data$b0),
    EmpiricalSE(nsims, estimates_data$b1),
    EmpiricalSE(nsims, estimates_data$log_alpha)
  ),
  MCE_Empirical_SE = c(
    MCE_empirical(nsims, estimates_data$b0),
    MCE_empirical(nsims, estimates_data$b1),
    MCE_empirical(nsims, estimates_data$log_alpha)
  ),
  MSE = c(
    MSE(nsims, estimates_data$b0, true_b0),
    MSE(nsims, estimates_data$b1, true_b1),
    MSE(nsims, estimates_data$log_alpha, true_log_alpha)
  ),
  MCE_MSE = c(
    MCE_mse(nsims, estimates_data$b0, true_b0),
    MCE_mse(nsims, estimates_data$b1, true_b1),
    MCE_mse(nsims, estimates_data$log_alpha, true_log_alpha)
  ),
  Coverage = c(
    coverage(ci_lower(estimates_data$b0, estimates_data$se_b0, 0.025),
             ci_upper(estimates_data$b0, estimates_data$se_b0, 0.025), true_b0, nsims),
    coverage(ci_lower(estimates_data$b1, estimates_data$se_b1, 0.025),
             ci_upper(estimates_data$b1, estimates_data$se_b1, 0.025), true_b1, nsims),
    coverage(ci_lower(estimates_data$log_alpha, estimates_data$se_log_alpha, 0.025),
             ci_upper(estimates_data$log_alpha, estimates_data$se_log_alpha, 0.025), true_log_alpha, nsims)
  ),
  MCE_Coverage = c(
    MCE_coverage(coverage(ci_lower(estimates_data$b0, estimates_data$se_b0, 0.025),
                          ci_upper(estimates_data$b0, estimates_data$se_b0, 0.025), true_b0, nsims), nsims),
    MCE_coverage(coverage(ci_lower(estimates_data$b1, estimates_data$se_b1, 0.025),
                          ci_upper(estimates_data$b1, estimates_data$se_b1, 0.025), true_b1, nsims), nsims),
    MCE_coverage(coverage(ci_lower(estimates_data$log_alpha, estimates_data$se_log_alpha, 0.025),
                          ci_upper(estimates_data$log_alpha, estimates_data$se_log_alpha, 0.025), true_log_alpha, nsims), nsims)
  )
)

print(performance_measures, digits = 2)
saveRDS(performance_measures, file.path(results_folder, paste0('performance_measures_', sample_size, '.rds')))

#==============================================================================#
# 6. Visualisation of results                                                 ----
#==============================================================================#

dir.create(file.path(results_folder, "plots"), recursive = TRUE)
plots_dir <- file.path(results_folder, "plots")

# Contraction vs Z-Score
z1 <- plot_contraction_vs_zscore(estimates_data, "contraction_b0",  "z_score_b0",  expression(beta[0]),    0,     "Contraction_vs_Zscore_B0.pdf")
z2 <- plot_contraction_vs_zscore(estimates_data, "contraction_b1", "z_score_b1", expression(beta[1]),    1,    "Contraction_vs_Zscore_B1.pdf")
z3 <- plot_contraction_vs_zscore(estimates_data, "contraction_log_alpha", "z_score_log_alpha", expression(log(alpha)), "log_alpha", "Contraction_vs_Zscore_log_alpha.pdf")
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Contraction_vs_Zscore.pdf")

# Zipper Plots (log scale + original scale for alpha)
z1 <- zipper_plot("b0", true_b0, estimates_data, nsims, 0) + theme(legend.position = "none")
z2 <- zipper_plot("b1", true_b1, estimates_data, nsims, 1) + theme(legend.position = "none")
z3 <- zipper_plot("log_alpha", true_log_alpha, estimates_data, nsims, "log_alpha") + theme(legend.position = "none")
z4 <- zipper_plot_original_scale("log_alpha", true_log_alpha, estimates_data, nsims) + theme(legend.position = "none")
save_grid_plot(grid.arrange(z1, z2, z3, z4, ncol = 4), "Zipper_Plot.pdf", width = 16, height = 6)

# Histogram of MAP Estimates
z1 <- plot_histogram(estimates_data, "b0", true_b0, expression(beta[0]), "Hist_MAP_B0.pdf", 0, nsims)
z2 <- plot_histogram(estimates_data, "b1", true_b1, expression(beta[1]), "Hist_MAP_B1.pdf", 1, nsims)
z3 <- plot_histogram(estimates_data, "log_alpha", true_log_alpha, expression(log(alpha)), "Hist_MAP_log_alpha.pdf", "log_alpha", nsims)
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Hist_MAP_Estimates.pdf")

# Scatter Plot of MAP Estimates vs SEs
z1 <- plot_scatter_map_vs_se(estimates_data, true_b0, "b0", "se_b0", expression(beta[0]), "Scatter_MAP_vs_SE_B0.pdf", 0)
z2 <- plot_scatter_map_vs_se(estimates_data, true_b1, "b1", "se_b1", expression(beta[1]), "Scatter_MAP_vs_SE_B1.pdf", 1)
z3 <- plot_scatter_map_vs_se(estimates_data, true_log_alpha, "log_alpha", "se_log_alpha", expression(log(alpha)), "Scatter_MAP_vs_SE_log_alpha.pdf", "log_alpha")
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Scatter_MAP_vs_SE.pdf")

# Histogram of SEs
z1 <- plot_histogram_se(estimates_data, "se_b0", expression(beta[0]), "Hist_SE_B0.pdf", 0)
z2 <- plot_histogram_se(estimates_data, "se_b1", expression(beta[1]), "Hist_SE_B1.pdf", 1)
z3 <- plot_histogram_se(estimates_data, "se_log_alpha", expression(log(alpha)), "Hist_SE_log_alpha.pdf", "log_alpha")
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Hist_SEs.pdf")

rm(z1, z2, z3, z4)

