# =============================================================================
# RCS_Source_Functions.R
# Author:  Hadiqa Tahir
# Date:    June 2026
# =============================================================================
# Helper functions for RCS_Model_Simulation_Study.R
#   - TMB model wrapper (cumulative_sum, TMB_model)
#   - Simulation helper (simulate_disease_data)
#   - Performance measures (Bias, EmpiricalSE, MSE, Coverage + MCEs)
#   - Plot functions (zipper, contraction, beta estimates, histogram, scatter, SE)
# =============================================================================

#-----------------------------------------------------------------------#
#               FOR RUNNING MODEL                                 ----
#-----------------------------------------------------------------------#

# Calculates the cumulative sum and adjusts indexing to .cpp style
cumulative_sum <- function(vec) {
  cumul_sum <- cumsum(c(0, vec))  
  cumul_sum <- head(cumul_sum, -1)  
  return(cumul_sum)
}

# Wrapper function to run model
TMB_model <- function(cases, ndvi_and_cov, weights, startindex, n_pixels_vec) {
  
  datalist <- list(
    
    # Data 
    cases        = cases,
    x            = as.matrix(ndvi_and_cov),
    weights      = weights,
    startindex   = startindex,
    n_pixels_vec = n_pixels_vec,
    
    # Priors for b0, b1, b2 ~ N(0, 10^2)
    priormean_intercept = 0,
    priorsd_intercept   = 10,
    
    priormean_slope = 0,
    priorsd_slope   = 10
  )
  print('Data Loaded.')
  print(names(datalist))
  
  # slope is a vector of length ncol(ndvi_and_cov) = 2 for RCS (two spline basis columns)
  parameters <- list(intercept = 0, slope = rep(0, ncol(ndvi_and_cov)))
  print('Parameters Loaded.')
  
  obj <- MakeADFun(datalist, parameters, DLL = "Fixed_weight_model")
  print('Objective Function Loaded.')
  
  opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, 
               method = 'BFGS',
               control = list(maxit = 200))
  print('Optimised.')
  
  hess <- stats::optimHess(par = opt$par, fn = obj$fn, gr = obj$gr, maxit = 300)  
  print('Calculated Hessian Matrix.')
  
  unc <- TMB::sdreport(obj, getJointPrecision = TRUE, hessian.fixed = hess)
  print('Calculated Uncertainty.')
  
  # Extract REPORT() values from the model
  nll_priors        <- obj$report()$nll_priors
  nll_final         <- obj$report()$nll_final
  pixel_linear_pred <- obj$report()$pixel_linear_pred
  p                 <- obj$report()$p
  
  results_list <- list(
    obj               = obj,                  # comment out to reduce memory if needed
    opt               = opt,                  #
    par.fixed         = unc$par.fixed,
    cov.fixed         = unc$cov.fixed,
    hess              = hess,
    unc               = unc,                  #
    nll_priors        = nll_priors,
    nll_final         = nll_final,
    pixel_linear_pred = pixel_linear_pred,    #
    p                 = p                     #
  )
  
  rm(datalist, parameters, obj, opt, hess, unc)
  return(results_list)
}

#-----------------------------------------------------------------------#
#               FOR SIMULATING DATA                                 ----
#-----------------------------------------------------------------------#

simulate_disease_data <- function(df, seed = 4000, i) {
  set.seed(seed + i)
  ARI_simulated <- as.logical(rbinom(n    = length(unique(df$HOUSE)),
                                     size = 1,
                                     prob = df$p_ari[match(unique(df$HOUSE), df$HOUSE)]))
  return(ARI_simulated)
}

#-----------------------------------------------------------------------#
#               FOR PERFORMANCE MEASURES                            ----
#-----------------------------------------------------------------------#

# Bias 
Bias <- function(b.hat, b.true, n) {
  bias <- (1/n)*sum(b.hat - b.true)
  return(bias)
}

# Monte Carlo Estimate for Bias 
MCE_bias <- function(n, b.i){
  VB  <- (sum((b.i - mean(b.i))^2))
  MCE <- sqrt(VB/(n*(n-1)))
  return(MCE)
}

# Empirical SE
EmpiricalSE <- function(n, b.i){
  VB  <- sum((b.i - mean(b.i))^2)
  ESE <- sqrt((1/ (n-1))* VB)
  return(ESE)
}

# Monte Carlo Estimate for Empirical SE 
MCE_empirical <- function(n, b.i){
  ESE  <- EmpiricalSE(n, b.i)
  MCSE <- ESE/sqrt(2 * (n-1))
  return(MCSE)
}

# Mean Square Error
MSE <- function(n, b.i, true.b) {
  MSE <- (1/n)*sum((b.i - true.b)^2)
  return(MSE)
}

# Monte Carlo Estimate for Mean Square Error
MCE_mse <- function(n, b.i, true.b){
  VB  <- sum((((b.i - true.b)^2) - MSE(n, b.i, true.b))^2)
  MCE <- sqrt(VB/(n*(n-1))) 
  return(MCE)
}

# Credible Intervals 
ci_lower <- function(b, se.b, alpha){
  z_quantile <- qnorm(1 - alpha)
  ci <- b - z_quantile*se.b
  return(ci)
}

ci_upper <- function(b, se.b, alpha){
  z_quantile <- qnorm(1 - alpha)
  ci <- b + z_quantile*se.b
  return(ci)
}

# Coverage
coverage <- function(theta_low_hat, theta_upp_hat, theta_true, n) {
  cov <- (1/n)*sum(theta_low_hat <= theta_true & theta_true <= theta_upp_hat)
  return(cov)
}

# Monte Carlo Error of Coverage 
MCE_coverage <- function(coverage, n){
  MCE <- sqrt( (coverage*(1-coverage))/ n)
  return(MCE)
}

#-----------------------------------------------------------------------#
#               TO GENERATE PLOTS                                  ----
#-----------------------------------------------------------------------#

zipper_plot <- function(parameter_name, true_value, data, n, index) {
  
  est <- data[, c(parameter_name, paste0("ci_lower_", parameter_name), 
                  paste0("ci_upper_", parameter_name), paste0("se_", parameter_name))]
  est$Simulation <- 1:n
  
  est$within_CI <- ifelse(est[[paste0("ci_lower_", parameter_name)]] <= true_value & 
                            est[[paste0("ci_upper_", parameter_name)]] >= true_value, TRUE, FALSE)
  
  est$zi <- (est[[parameter_name]] - true_value) / est[[paste0("se_", parameter_name)]]
  est    <- est[order(abs(est$zi)), ]
  est$rank <- 1:nrow(est)
  
  title_text <- bquote("Ranked Credible Intervals for " * beta[.(index)]) 
  
  plot <- ggplot(est, aes(x = .data[[parameter_name]], y = rank)) +
    geom_pointrange(aes(xmin = .data[[paste0("ci_lower_", parameter_name)]], 
                        xmax = .data[[paste0("ci_upper_", parameter_name)]], color = within_CI),
                    size = 0.5, fatten = 1) + 
    geom_vline(xintercept = true_value, linetype = "dashed", color = "black") +
    annotate("text", x = true_value, y = max(est$rank),
             label = paste0("\u03B2", index, " = ", format(true_value, digits = 3, nsmall = 2)),
             colour = "black", hjust = 1.1, vjust = 0) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"), 
                       name = "Within CI", 
                       breaks = c(TRUE, FALSE),
                       labels = c("True Value Within CI", "True Value Outside CI")) +
    scale_y_continuous(limits = c(1, max(est$rank))) +
    theme_minimal() +
    labs(title = title_text, x = "MAP Estimate",
         y = expression("Rank Based on "*abs(z[i])))
  
  return(plot)
}

#-----------------------------------------------------------------------#
#               CONTRACTION VS Z-SCORE PLOT 
#-----------------------------------------------------------------------#

plot_contraction_vs_zscore <- function(estimates_data, contraction_value, z_score_value, beta_label, index, save_name) {
  
  l_xaxis <- bquote("Posterior Contraction")
  l_yaxis <- bquote("Posterior Z-score")
  l_title <- bquote("Contraction vs Z-score for " * beta[.(index)])
  col_plot = "#72874EFF"
  
  contraction_zscore_plot <- ggplot(estimates_data, aes_string(x = contraction_value, y = z_score_value)) +
    geom_point(alpha = 0.6, col = col_plot) +
    labs(x = l_xaxis, y = l_yaxis, title = l_title) +
    theme_minimal() +
    xlim(0, 1)
}

#-----------------------------------------------------------------------#
#          BETA ESTIMATE + 95% CREDIBLE INTERVAL PLOT 
#-----------------------------------------------------------------------#

plot_beta_estimates <- function(estimates_data, beta_value, true_beta, ci_lower, ci_upper, simulation_column, beta_label, save_name, index) {
  
  y_range <- range(estimates_data[[ci_lower]], estimates_data[[ci_upper]])
  offset  <- (y_range[2] - y_range[1]) * 0.04
  
  l_title  <- bquote(beta[.(index)] ~ "MAP Estimate and 95% Credible Intervals")
  l_yaxis  <- bquote(beta[.(index)] ~ "MAP Estimate")
  col_plot = "#72874EFF"
  col_line = "black"
  
  error_plot <- ggplot() + 
    geom_errorbar(data = estimates_data,
                  aes_string(x = simulation_column, ymin = ci_lower, ymax = ci_upper), 
                  color = col_plot, alpha = 0.7) + 
    geom_point(data = estimates_data, 
               aes_string(x = simulation_column, y = beta_value),
               colour = col_plot) + 
    geom_hline(yintercept = true_beta, lty = "dashed", colour = col_line) +
    annotate("text", x = 0, y = true_beta + offset,
             label = paste0("True \u03B2", index, " = ", format(true_beta, digits = 3, nsmall = 2)),
             colour = col_line) + 
    theme_minimal() + 
    labs(title = l_title, y = l_yaxis)
}

#-----------------------------------------------------------------------#
#               HISTOGRAM OF MAP ESTIMATES
#-----------------------------------------------------------------------#

plot_histogram <- function(estimates_data, beta_value, true_beta, beta_label, save_name, index, n_sims) {
  
  x_range <- range(estimates_data[[beta_value]])
  offset  <- (x_range[2] - x_range[1]) * 0.12
  
  l_xaxis  <- bquote(beta[.(index)])
  l_title  <- bquote("Histogram of " * beta[.(index)] * " MAP estimates over " * .(n_sims) * " simulations")
  col_line = "#023743FF"
  col_fill = "#72874EFF"
  col_plot = "#453947FF"
  
  hist_plot <- ggplot(data = estimates_data, aes_string(beta_value)) + 
    geom_histogram(binwidth = 0.1, fill = col_fill, color = col_plot, alpha = 0.7) + 
    theme_minimal() + 
    labs(x = l_xaxis, y = "Frequency", title = l_title) + 
    geom_vline(xintercept = true_beta, linetype = "dashed", color = col_line) + 
    annotate("text", x = true_beta + offset, y = Inf,
             label = paste0("True \u03B2", index, " = ", format(true_beta, digits = 3, nsmall = 2)),
             colour = col_line, vjust = 2)
}

#-----------------------------------------------------------------------#
#               SCATTER PLOT OF MAP VS SE  
#-----------------------------------------------------------------------#

plot_scatter_map_vs_se <- function(estimates_data, true_beta, beta_value, se_value, beta_label, save_name, index) {
  
  y_range  <- range(estimates_data[[se_value]])
  offset_y <- (y_range[2] - y_range[1]) * 0.09
  x_range  <- range(estimates_data[[beta_value]])
  offset_x <- (x_range[2] - x_range[1]) * 0.15
  
  l_xaxis  <- bquote(beta[.(index)] * " MAP estimates")
  l_yaxis  <- bquote(beta[.(index)] * " SE")
  l_title  <- bquote(beta[.(index)] * " MAP Estimates vs. Standard Error")
  col_plot = "#72874EFF"
  col_line = "#023743FF"
  
  scatter_plot <- ggplot(estimates_data, aes_string(beta_value, se_value)) + 
    geom_point(col = col_plot) + 
    theme_minimal() + 
    labs(x = l_xaxis, y = l_yaxis, title = l_title) + 
    geom_vline(xintercept = true_beta, lty = "dashed", colour = col_line) + 
    annotate("text",
             x = true_beta + 0.7 * offset_x,
             y = offset_y + max(estimates_data[[se_value]]),
             label = paste0("True \u03B2", index, " = ", format(true_beta, digits = 3, nsmall = 2)),
             colour = col_line)
}

#-----------------------------------------------------------------------#
#             HISTOGRAM OF SES
#-----------------------------------------------------------------------#

plot_histogram_se <- function(estimates_data, se_value, beta_label, save_name, index) {
  
  l_xaxis  <- bquote(beta[.(index)] * " SE")
  l_title  <- bquote("Histogram of Standard Errors for " * beta[.(index)])
  col_line = "#023743FF"
  col_fill = "#72874EFF"
  
  hist_se_plot <- ggplot(estimates_data, aes_string(se_value)) + 
    geom_histogram(fill = col_fill, colour = col_line, alpha = 0.7) + 
    theme_minimal() + 
    labs(x = l_xaxis, y = "Frequency", title = l_title)
}

#-----------------------------------------------------------------------#
#             SAVE PLOTS AS A GRID
#-----------------------------------------------------------------------#

save_grid_plot <- function(grid_plot, filename, width = 12, height = 8) {
  ggsave(
    filename = file.path(plots_dir, filename),
    plot     = grid_plot,
    device   = "pdf",
    width    = width,
    height   = height
  )
}
