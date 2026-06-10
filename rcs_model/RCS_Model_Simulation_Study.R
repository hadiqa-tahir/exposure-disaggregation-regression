# =============================================================================#
# RCS_Model_Simulation_Study.R
# Author:  Hadiqa Tahir
# Date:    June 2026
# =============================================================================#
# Research Questions:
#       - Can the model recover parameter estimates?
#       - How does sample size affect model performance?
#
# Restricted Cubic Spline (RCS) model details:
#   logit P(ARI_j) = sum_i [ w_i * (b0 + b1*basis1_i + b2*basis2_i) ]
#   where basis columns are splines::ns(NDVI, knots = median, Boundary.knots = c(0.1, 0.99))
#   and   w_i = (1/d_i) / sum(1/d_i)  (inverse-distance weights, alpha = 1 fixed)
#   Priors: b0, b1, b2 ~ N(0, 10^2)
#
#   Note: Uses the same Fixed_weight_model.cpp as the FW model — the RCS
#   transformation is applied to the exposure matrix before passing to TMB.
#
# This file covers:
#   1. Fitting the RCS model to 'observed' data to obtain "true" parameters.
#   2. Simulating ARI outcomes from the fitted model (nsims replicates).
#   3. Refitting the model on each simulated dataset.
#   4. Extracting MAP estimates, SEs, and covariance matrices.
#   5. Computing performance measures with Monte Carlo errors.
#   6. Producing diagnostic plots (zipper, histogram, scatter, SE distribution).
#
# Dependencies:
#   - RCS_Source_Functions.R               : helper functions
#   - Fixed_weight_model.cpp               : TMB C++ model (shared with FW model)
#   - df_<sample_size>_synthetic.parquet   : synthetic dataset mimicking real DHS data
#
# Outputs (written to rcs_model/results_<sample_size>/):
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
library(splines)

# Set working directory
setwd(here::here("rcs_model"))

# Functions
source('RCS_Source_Functions.R')

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

# Exposure data: transform NDVI to RCS basis columns
# Knots are computed once from observed data and reused in all subsequent sections
# to ensure a consistent basis across the observed fit and all simulation refits.
ndvi_raw <- df %>% select(NDVI) %>% collect() %>% as.matrix()
knots  <- quantile(ndvi_raw, probs = c(0.5),        names = FALSE)   # 1 internal knot
bknots <- quantile(ndvi_raw, probs = c(0.1, 0.99),  names = FALSE)   # boundary knots
ndvi_and_cov <- splines::ns(ndvi_raw, knots = knots, intercept = FALSE, Boundary.knots = bknots)
rm(ndvi_raw)

# Weights as a vector
weights <- df %>% select(WEIGHT) %>% collect() %>% pull()
# check: sum(df$WEIGHT[df$HOUSE == 4]) should equal 1

# Compile and load model
orig_wd <- getwd()
setwd(here::here("fixed_weight_model"))
compile("Fixed_weight_model.cpp", "&> logfile.log")
dyn.load(dynlib("Fixed_weight_model"))
setwd(orig_wd)
rm(orig_wd)

# Run model
results <- TMB_model(ARI, ndvi_and_cov, weights, startindex, n_pixels_vec)
print(results)

# Save the observed model fit
path <- paste0('observed_fit_', sample_size, '.rds')
saveRDS(results, file = path)

rm(ARI, data_path, path, ndvi_raw)

#==============================================================================#
# 2. Simulate data from observed model fit                                    ----
#==============================================================================#

nsims <- 10
dir.create(paste0("simulated_outcomes_", sample_size), showWarnings = TRUE)
dir.create(paste0("results_", sample_size),            showWarnings = TRUE)

# Extract true parameters from observed fit
# logit P(ARI_j) = sum_i [ w_i * (b0 + b1*basis1_i + b2*basis2_i) ]
true_b0 <- results$par.fixed[['intercept']]
true_b1 <- results$par.fixed[[2]]
true_b2 <- results$par.fixed[[3]]

spline_coefs <- c(true_b1, true_b2)
Restricted_Cubic_Spline <- true_b0 + ndvi_and_cov %*% spline_coefs
n_j   <- with(df, stats::ave(WEIGHT * as.numeric(Restricted_Cubic_Spline), HOUSE, FUN = sum, na.rm = TRUE))
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

rm(df_nj_pj, n_j, p_ari, Restricted_Cubic_Spline, spline_coefs,
   output_folder, sim_index, simulation_result, file_name)

#==============================================================================#
# 3. Run model on simulated data                                              ----
#==============================================================================#

# All inputs are the same as section 1; only the outcome vector changes.
simulation_file <- file.path(paste0("simulated_outcomes_", sample_size), paste0("Sim_", 1:nsims, ".parquet"))
results_folder  <- paste0("results_", sample_size)

for (i in 1:nsims) {
  ARI    <- arrow::open_dataset(simulation_file[i]) %>% collect()
  cases  <- as.vector(as.integer(ARI[[1]]))
  results <- TMB_model(cases, ndvi_and_cov, weights, startindex, n_pixels_vec)
  name   <- gsub(".parquet", "_results.rds", basename(simulation_file[i]))
  output_file <- file.path(results_folder, name)
  saveRDS(results, output_file)
  message("Saved ", output_file)
}

rm(ndvi_and_cov, weights, startindex, n_pixels_vec,
   simulation_file, ARI, cases, results, name, output_file, i)

#==============================================================================#
# 4. Extract results from simulation model runs                               ----
#==============================================================================#

files          <- list.files(path = results_folder, pattern = "Sim_.*_results.rds", full.names = TRUE)
estimates_data <- data.frame()
for (i in seq_along(files)) {
  r <- readRDS(files[i])
  estimates_data <- rbind(estimates_data, data.frame(
    Simulation = i,
    b0         = r$par.fixed[[1]],
    b1         = r$par.fixed[[2]],
    b2         = r$par.fixed[[3]],
    se_b0      = sqrt(diag(solve(r$hess)))[[1]],
    se_b1      = sqrt(diag(solve(r$hess)))[[2]],
    se_b2      = sqrt(diag(solve(r$hess)))[[3]],
    nllprior   = r$nll_priors,
    nllfinal   = r$nll_final
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
estimates_data$ci_lower_b0 <- ci_lower(estimates_data$b0, estimates_data$se_b0, 0.025)
estimates_data$ci_upper_b0 <- ci_upper(estimates_data$b0, estimates_data$se_b0, 0.025)
estimates_data$ci_lower_b1 <- ci_lower(estimates_data$b1, estimates_data$se_b1, 0.025)
estimates_data$ci_upper_b1 <- ci_upper(estimates_data$b1, estimates_data$se_b1, 0.025)
estimates_data$ci_lower_b2 <- ci_lower(estimates_data$b2, estimates_data$se_b2, 0.025)
estimates_data$ci_upper_b2 <- ci_upper(estimates_data$b2, estimates_data$se_b2, 0.025)

# Variance of each parameter from covariance matrix
for (i in 1:nsims) {
  estimates_data$b0_var[i] <- covariance_matrix[[i]][1, 1]
  estimates_data$b1_var[i] <- covariance_matrix[[i]][2, 2]
  estimates_data$b2_var[i] <- covariance_matrix[[i]][3, 3]
}
rm(i, covariance_matrix)

# Contraction and z-scores
sd_prior <- 10      # Priors for b0, b1, b2 were N(0, 10^2)
estimates_data <- estimates_data %>%
  mutate(
    contraction_b0 = 1 - (b0_var / sd_prior^2),
    contraction_b1 = 1 - (b1_var / sd_prior^2),
    contraction_b2 = 1 - (b2_var / sd_prior^2),
    z_score_b0     = (b0 - true_b0) / sqrt(b0_var),
    z_score_b1     = (b1 - true_b1) / sqrt(b1_var),
    z_score_b2     = (b2 - true_b2) / sqrt(b2_var)
  )
rm(sd_prior)

# Bias, Empirical SE, MSE, Coverage
performance_measures <- data.frame(
  Parameter = c("B0", "B1", "B2"),
  True_Value = c(true_b0, true_b1, true_b2),
  Bias = c(
    Bias(estimates_data$b0, true_b0, nsims),
    Bias(estimates_data$b1, true_b1, nsims),
    Bias(estimates_data$b2, true_b2, nsims)
  ),
  MCE_Bias = c(
    MCE_bias(nsims, estimates_data$b0),
    MCE_bias(nsims, estimates_data$b1),
    MCE_bias(nsims, estimates_data$b2)
  ),
  Empirical_SE = c(
    EmpiricalSE(nsims, estimates_data$b0),
    EmpiricalSE(nsims, estimates_data$b1),
    EmpiricalSE(nsims, estimates_data$b2)
  ),
  MCE_Empirical_SE = c(
    MCE_empirical(nsims, estimates_data$b0),
    MCE_empirical(nsims, estimates_data$b1),
    MCE_empirical(nsims, estimates_data$b2)
  ),
  MSE = c(
    MSE(nsims, estimates_data$b0, true_b0),
    MSE(nsims, estimates_data$b1, true_b1),
    MSE(nsims, estimates_data$b2, true_b2)
  ),
  MCE_MSE = c(
    MCE_mse(nsims, estimates_data$b0, true_b0),
    MCE_mse(nsims, estimates_data$b1, true_b1),
    MCE_mse(nsims, estimates_data$b2, true_b2)
  ),
  Coverage = c(
    coverage(ci_lower(estimates_data$b0, estimates_data$se_b0, 0.025),
             ci_upper(estimates_data$b0, estimates_data$se_b0, 0.025), true_b0, nsims),
    coverage(ci_lower(estimates_data$b1, estimates_data$se_b1, 0.025),
             ci_upper(estimates_data$b1, estimates_data$se_b1, 0.025), true_b1, nsims),
    coverage(ci_lower(estimates_data$b2, estimates_data$se_b2, 0.025),
             ci_upper(estimates_data$b2, estimates_data$se_b2, 0.025), true_b2, nsims)
  ),
  MCE_Coverage = c(
    MCE_coverage(coverage(ci_lower(estimates_data$b0, estimates_data$se_b0, 0.025),
                          ci_upper(estimates_data$b0, estimates_data$se_b0, 0.025), true_b0, nsims), nsims),
    MCE_coverage(coverage(ci_lower(estimates_data$b1, estimates_data$se_b1, 0.025),
                          ci_upper(estimates_data$b1, estimates_data$se_b1, 0.025), true_b1, nsims), nsims),
    MCE_coverage(coverage(ci_lower(estimates_data$b2, estimates_data$se_b2, 0.025),
                          ci_upper(estimates_data$b2, estimates_data$se_b2, 0.025), true_b2, nsims), nsims)
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
z1 <- plot_contraction_vs_zscore(estimates_data, "contraction_b0", "z_score_b0", expression(beta[0]), 0, "Contraction_vs_Zscore_B0.pdf")
z2 <- plot_contraction_vs_zscore(estimates_data, "contraction_b1", "z_score_b1", expression(beta[1]), 1, "Contraction_vs_Zscore_B1.pdf")
z3 <- plot_contraction_vs_zscore(estimates_data, "contraction_b2", "z_score_b2", expression(beta[2]), 2, "Contraction_vs_Zscore_B2.pdf")
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Contraction_vs_Zscore.pdf")

# Zipper Plots
z1 <- zipper_plot("b0", true_b0, estimates_data, nsims, 0) + theme(legend.position = "none")
z2 <- zipper_plot("b1", true_b1, estimates_data, nsims, 1) + theme(legend.position = "none")
z3 <- zipper_plot("b2", true_b2, estimates_data, nsims, 2) + theme(legend.position = "none")
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Zipper_Plot.pdf")

# Histogram of MAP Estimates
z1 <- plot_histogram(estimates_data, "b0", true_b0, expression(beta[0]), "Hist_MAP_B0.pdf", 0, nsims)
z2 <- plot_histogram(estimates_data, "b1", true_b1, expression(beta[1]), "Hist_MAP_B1.pdf", 1, nsims)
z3 <- plot_histogram(estimates_data, "b2", true_b2, expression(beta[2]), "Hist_MAP_B2.pdf", 2, nsims)
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Hist_MAP_Estimates.pdf")

# Scatter Plot of MAP Estimates vs SEs
z1 <- plot_scatter_map_vs_se(estimates_data, true_b0, "b0", "se_b0", expression(beta[0]), "Scatter_MAP_vs_SE_B0.pdf", 0)
z2 <- plot_scatter_map_vs_se(estimates_data, true_b1, "b1", "se_b1", expression(beta[1]), "Scatter_MAP_vs_SE_B1.pdf", 1)
z3 <- plot_scatter_map_vs_se(estimates_data, true_b2, "b2", "se_b2", expression(beta[2]), "Scatter_MAP_vs_SE_B2.pdf", 2)
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Scatter_MAP_vs_SE.pdf")

# Histogram of SEs
z1 <- plot_histogram_se(estimates_data, "se_b0", expression(beta[0]), "Hist_SE_B0.pdf", 0)
z2 <- plot_histogram_se(estimates_data, "se_b1", expression(beta[1]), "Hist_SE_B1.pdf", 1)
z3 <- plot_histogram_se(estimates_data, "se_b2", expression(beta[2]), "Hist_SE_B2.pdf", 2)
save_grid_plot(grid.arrange(z1, z2, z3, ncol = 3), "Hist_SEs.pdf")

rm(z1, z2, z3)

# Spline Effect Plot (observed fit)
ndvi_seq   <- seq(min(df$NDVI), max(df$NDVI), length.out = 200)
X_pred     <- splines::ns(ndvi_seq, knots = knots, Boundary.knots = bknots, intercept = FALSE)
obs        <- readRDS(paste0('observed_fit_', sample_size, '.rds'))
beta_spline <- obs$par.fixed[2:length(obs$par.fixed)]
intercept   <- obs$par.fixed[['intercept']]
pred       <- as.numeric(intercept + X_pred %*% beta_spline)
X_pred_full <- cbind(1, X_pred)
pred_se    <- sqrt(diag(X_pred_full %*% obs$cov.fixed %*% t(X_pred_full)))

spline_plot_data <- data.frame(
  NDVI   = ndvi_seq,
  Effect = pred,
  lower  = pred - 1.96 * pred_se,
  upper  = pred + 1.96 * pred_se
)

spline_plot <- ggplot(spline_plot_data, aes(x = NDVI, y = Effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#72874EFF", alpha = 0.3) +
  geom_line(color = "#023743FF", linewidth = 1.2) +
  geom_vline(xintercept = knots,  linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = bknots, linetype = "dotted", color = "#E69F00") +
  labs(
    title    = "Estimated Spline Effect of NDVI on ARI",
    subtitle = "Dashed = internal knot, Dotted = boundary knots",
    x        = "NDVI",
    y        = "Linear predictor"
  ) +
  theme_minimal()

spline_plot
save_grid_plot(spline_plot, "Spline_Effect_Observed.pdf", width = 8, height = 6)

# Spline Effect Plot on Simulated Fits 
X_pred   <- splines::ns(ndvi_seq, knots = knots, Boundary.knots = bknots, intercept = FALSE)
sim_files <- list.files(path = results_folder, pattern = "Sim_.*_results.rds", full.names = TRUE)
sim_curves <- lapply(seq_along(sim_files), function(i) {
  r            <- readRDS(sim_files[[i]])
  beta_spline  <- r$par.fixed[2:length(r$par.fixed)]
  intercept    <- r$par.fixed[['intercept']]
  data.frame(
    NDVI       = ndvi_seq,
    Effect     = as.numeric(intercept + X_pred %*% beta_spline),
    Simulation = i
  )
})

sim_curves <- do.call(rbind, sim_curves)
spline_sim_plot <- ggplot(sim_curves, aes(x = NDVI, y = Effect)) +
  geom_line(data = spline_plot_data, aes(x = NDVI, y = Effect),
            color = "#023743FF", linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE) +
  geom_line(color = "#72874EFF", linewidth = 0.8) +
  geom_vline(xintercept = knots,  linetype = "dashed", color = "gray50",  linewidth = 0.4) +
  geom_vline(xintercept = bknots, linetype = "dotted", color = "#E69F00", linewidth = 0.4) +
  facet_wrap(~ Simulation, ncol = 5, labeller = label_both) +
  labs(
    title    = "Spline Effect of NDVI Across Simulations",
    subtitle = "Green = simulated fit, Blue dashed = observed fit",
    x        = "NDVI",
    y        = "Linear predictor"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 8))
spline_sim_plot
save_grid_plot(spline_sim_plot, "Spline_Effect_Simulations.pdf", width = 8, height = 6)
rm(ndvi_seq, X_pred, sim_files, sim_curves, spline_sim_plot)