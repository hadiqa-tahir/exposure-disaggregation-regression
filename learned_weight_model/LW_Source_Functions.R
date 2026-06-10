# =============================================================================
# LW_Source_Functions.R
# Author:  Hadiqa Tahir
# Date:    June 2026
# =============================================================================
# Helper functions for LW_Model_Simulation_Study.R
#   - TMB model wrapper (cumulative_sum, TMB_model)
#   - Simulation helper (simulate_disease_data)
#   - Performance measures (Bias, EmpiricalSE, MSE, Coverage + MCEs)
#   - Plot functions (zipper, zipper original scale, contraction, beta estimates,
#                     histogram, scatter, SE histogram)
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

# Wrapper function to run the learned-weight TMB model
TMB_model <- function(cases, ndvi_and_cov, distances, startindex, n_pixels_vec,
                      priormean_log_alpha, priorsd_log_alpha) {
  
  datalist <- list(
    # Data
    cases        = cases,
    x            = as.matrix(ndvi_and_cov),
    distances    = distances,
    startindex   = startindex,
    n_pixels_vec = n_pixels_vec,
    
    # Priors: b0, b1 ~ N(0, 10^2);  log_alpha ~ N(priormean_log_alpha, priorsd_log_alpha^2)
    priormean_intercept = 0,
    priorsd_intercept   = 10,
    
    priormean_slope = 0,
    priorsd_slope   = 10,
    
    priormean_log_alpha = priormean_log_alpha,
    priorsd_log_alpha   = priorsd_log_alpha
  )
  print('Data Loaded.')
  print(names(datalist))
  
  parameters <- list(intercept = 0, slope = 0, log_alpha = 0.5)
  print('Parameters Loaded.')
  
  obj <- MakeADFun(datalist, parameters, DLL = "Learned_weights_model")
  print('Objective Function Loaded.')
  
  opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr,
               method = 'BFGS',
               control = list(maxit = 200, trace = 5))
  print('Optimised.')
  
  hess <- stats::optimHess(par = opt$par, fn = obj$fn, gr = obj$gr, maxit = 300)
  print('Calculated Hessian Matrix.')
  
  unc <- TMB::sdreport(obj, getJointPrecision = TRUE, hessian.fixed = hess)
  print('Calculated Uncertainty.')
  
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

Bias <- function(b.hat, b.true, n) {
  (1/n) * sum(b.hat - b.true)
}

MCE_bias <- function(n, b.i) {
  VB <- sum((b.i - mean(b.i))^2)
  sqrt(VB / (n * (n - 1)))
}

EmpiricalSE <- function(n, b.i) {
  VB <- sum((b.i - mean(b.i))^2)
  sqrt((1 / (n - 1)) * VB)
}

MCE_empirical <- function(n, b.i) {
  ESE <- EmpiricalSE(n, b.i)
  ESE / sqrt(2 * (n - 1))
}

MSE <- function(n, b.i, true.b) {
  (1/n) * sum((b.i - true.b)^2)
}

MCE_mse <- function(n, b.i, true.b) {
  VB <- sum((((b.i - true.b)^2) - MSE(n, b.i, true.b))^2)
  sqrt(VB / (n * (n - 1)))
}

ci_lower <- function(b, se.b, alpha) {
  b - qnorm(1 - alpha) * se.b
}

ci_upper <- function(b, se.b, alpha) {
  b + qnorm(1 - alpha) * se.b
}

coverage <- function(theta_low_hat, theta_upp_hat, theta_true, n) {
  (1/n) * sum(theta_low_hat <= theta_true & theta_true <= theta_upp_hat)
}

MCE_coverage <- function(coverage, n) {
  sqrt((coverage * (1 - coverage)) / n)
}

#-----------------------------------------------------------------------#
#               INTERNAL LABEL HELPERS                              ----
#-----------------------------------------------------------------------#

# bquote label for plot titles/axes — works for beta[i] and log(alpha)
.param_bquote <- function(index) {
  if (is.numeric(index)) bquote(beta[.(index)]) else bquote(log(alpha))
}

# Plain string label for annotate() — ggplot2 4.0.0 does not accept bquote() here
.param_annotate <- function(index, value) {
  if (is.numeric(index)) {
    paste0("\u03B2", index, " = ", format(value, digits = 3, nsmall = 2))
  } else {
    paste0("log(\u03B1) = ", format(value, digits = 3, nsmall = 2))
  }
}

#-----------------------------------------------------------------------#
#               ZIPPER PLOT                                         ----
#-----------------------------------------------------------------------#

zipper_plot <- function(parameter_name, true_value, data, n, index) {
  
  est <- data[, c(parameter_name, paste0("ci_lower_", parameter_name),
                  paste0("ci_upper_", parameter_name), paste0("se_", parameter_name))]
  est$Simulation <- 1:n
  
  est$within_CI <- ifelse(est[[paste0("ci_lower_", parameter_name)]] <= true_value &
                            est[[paste0("ci_upper_", parameter_name)]] >= true_value, TRUE, FALSE)
  
  est$zi   <- (est[[parameter_name]] - true_value) / est[[paste0("se_", parameter_name)]]
  est      <- est[order(abs(est$zi)), ]
  est$rank <- 1:nrow(est)
  
  title_text <- bquote("Ranked Credible Intervals for " * .(.param_bquote(index)))
  
  ggplot(est, aes(x = .data[[parameter_name]], y = rank)) +
    geom_pointrange(aes(xmin = .data[[paste0("ci_lower_", parameter_name)]],
                        xmax = .data[[paste0("ci_upper_", parameter_name)]], color = within_CI),
                    size = 0.5, fatten = 1) +
    geom_vline(xintercept = true_value, linetype = "dashed", color = "black") +
    annotate("text", x = true_value, y = max(est$rank),
             label  = .param_annotate(index, true_value),
             colour = "black", hjust = 1.1, vjust = 0) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       name   = "Within CI",
                       breaks = c(TRUE, FALSE),
                       labels = c("True Value Within CI", "True Value Outside CI")) +
    scale_y_continuous(limits = c(1, max(est$rank))) +
    theme_minimal() +
    labs(title = title_text, x = "MAP Estimate",
         y = expression("Rank Based on " * abs(z[i])))
}

# Zipper plot on the original (exponentiated) scale — used for log_alpha -> alpha
zipper_plot_original_scale <- function(parameter_name, true_value, data, n) {
  
  data[[parameter_name]]                          <- exp(data[[parameter_name]])
  data[[paste0("ci_lower_", parameter_name)]]     <- exp(data[[paste0("ci_lower_", parameter_name)]])
  data[[paste0("ci_upper_", parameter_name)]]     <- exp(data[[paste0("ci_upper_", parameter_name)]])
  true_value <- exp(true_value)
  
  est <- data[, c(parameter_name, paste0("ci_lower_", parameter_name),
                  paste0("ci_upper_", parameter_name), paste0("se_", parameter_name))]
  est$Simulation <- 1:n
  
  est$within_CI <- ifelse(est[[paste0("ci_lower_", parameter_name)]] <= true_value &
                            est[[paste0("ci_upper_", parameter_name)]] >= true_value, TRUE, FALSE)
  
  est$zi   <- (est[[parameter_name]] - true_value) / est[[paste0("se_", parameter_name)]]
  est      <- est[order(abs(est$zi)), ]
  est$rank <- 1:nrow(est)
  
  ggplot(est, aes(x = .data[[parameter_name]], y = rank)) +
    geom_pointrange(aes(xmin = .data[[paste0("ci_lower_", parameter_name)]],
                        xmax = .data[[paste0("ci_upper_", parameter_name)]], color = within_CI),
                    size = 0.5, fatten = 1) +
    geom_vline(xintercept = true_value, linetype = "dashed", color = "black") +
    annotate("text", x = true_value, y = max(est$rank),
             label  = paste0("\u03B1 = ", format(true_value, digits = 3, nsmall = 2)),
             colour = "black", hjust = 1.1, vjust = 0) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       name   = "Within CI",
                       breaks = c(TRUE, FALSE),
                       labels = c("True Value Within CI", "True Value Outside CI")) +
    scale_y_continuous(limits = c(1, max(est$rank))) +
    theme_minimal() +
    labs(title = bquote("Ranked Credible Intervals for " * alpha),
         x = expression(alpha ~ "(original scale)"),
         y = expression("Rank Based on " * abs(z[i])))
}

#-----------------------------------------------------------------------#
#               CONTRACTION VS Z-SCORE PLOT                         ----
#-----------------------------------------------------------------------#

plot_contraction_vs_zscore <- function(estimates_data, contraction_value, z_score_value, beta_label, index, save_name) {
  
  l_title <- bquote("Contraction vs Z-score for " * .(.param_bquote(index)))
  col_plot <- "#72874EFF"
  
  ggplot(estimates_data, aes_string(x = contraction_value, y = z_score_value)) +
    geom_point(alpha = 0.6, col = col_plot) +
    labs(x = "Posterior Contraction", y = "Posterior Z-score", title = l_title) +
    theme_minimal() +
    xlim(0, 1)
}

#-----------------------------------------------------------------------#
#          BETA ESTIMATE + 95% CREDIBLE INTERVAL PLOT               ----
#-----------------------------------------------------------------------#

plot_beta_estimates <- function(estimates_data, beta_value, true_beta, ci_lower, ci_upper,
                                simulation_column, beta_label, save_name, index) {
  
  y_range <- range(estimates_data[[ci_lower]], estimates_data[[ci_upper]])
  offset  <- (y_range[2] - y_range[1]) * 0.04
  
  l_title <- bquote(.(.param_bquote(index)) ~ "MAP Estimate and 95% Credible Intervals")
  l_yaxis <- bquote(.(.param_bquote(index)) ~ "MAP Estimate")
  col_plot <- "#72874EFF"
  col_line <- "black"
  
  ggplot() +
    geom_errorbar(data = estimates_data,
                  aes_string(x = simulation_column, ymin = ci_lower, ymax = ci_upper),
                  color = col_plot, alpha = 0.7) +
    geom_point(data = estimates_data,
               aes_string(x = simulation_column, y = beta_value),
               colour = col_plot) +
    geom_hline(yintercept = true_beta, lty = "dashed", colour = col_line) +
    annotate("text", x = 0, y = true_beta + offset,
             label  = .param_annotate(index, true_beta),
             colour = col_line) +
    theme_minimal() +
    labs(title = l_title, y = l_yaxis)
}

#-----------------------------------------------------------------------#
#               HISTOGRAM OF MAP ESTIMATES                          ----
#-----------------------------------------------------------------------#

plot_histogram <- function(estimates_data, beta_value, true_beta, beta_label, save_name, index, n_sims) {
  
  x_range <- range(estimates_data[[beta_value]])
  offset  <- (x_range[2] - x_range[1]) * 0.12
  
  l_xaxis <- .param_bquote(index)
  l_title <- bquote("Histogram of " * .(.param_bquote(index)) * " MAP estimates over " * .(n_sims) * " simulations")
  col_line <- "#023743FF"
  col_fill <- "#72874EFF"
  col_plot <- "#453947FF"
  
  ggplot(data = estimates_data, aes_string(beta_value)) +
    geom_histogram(binwidth = 0.1, fill = col_fill, color = col_plot, alpha = 0.7) +
    theme_minimal() +
    labs(x = l_xaxis, y = "Frequency", title = l_title) +
    geom_vline(xintercept = true_beta, linetype = "dashed", color = col_line) +
    annotate("text", x = true_beta + offset, y = Inf,
             label  = .param_annotate(index, true_beta),
             colour = col_line, vjust = 2)
}

#-----------------------------------------------------------------------#
#               SCATTER PLOT OF MAP VS SE                           ----
#-----------------------------------------------------------------------#

plot_scatter_map_vs_se <- function(estimates_data, true_beta, beta_value, se_value,
                                   beta_label, save_name, index) {
  
  y_range  <- range(estimates_data[[se_value]])
  offset_y <- (y_range[2] - y_range[1]) * 0.09
  x_range  <- range(estimates_data[[beta_value]])
  offset_x <- (x_range[2] - x_range[1]) * 0.15
  
  l_xaxis <- bquote(.(.param_bquote(index)) * " MAP estimates")
  l_yaxis <- bquote(.(.param_bquote(index)) * " SE")
  l_title <- bquote(.(.param_bquote(index)) * " MAP Estimates vs. Standard Error")
  col_plot <- "#72874EFF"
  col_line <- "#023743FF"
  
  ggplot(estimates_data, aes_string(beta_value, se_value)) +
    geom_point(col = col_plot) +
    theme_minimal() +
    labs(x = l_xaxis, y = l_yaxis, title = l_title) +
    geom_vline(xintercept = true_beta, lty = "dashed", colour = col_line) +
    annotate("text",
             x      = true_beta + 0.7 * offset_x,
             y      = offset_y + max(estimates_data[[se_value]]),
             label  = .param_annotate(index, true_beta),
             colour = col_line)
}

#-----------------------------------------------------------------------#
#               HISTOGRAM OF SES                                    ----
#-----------------------------------------------------------------------#

plot_histogram_se <- function(estimates_data, se_value, beta_label, save_name, index) {
  
  l_xaxis <- bquote(.(.param_bquote(index)) * " SE")
  l_title <- bquote("Histogram of Standard Errors for " * .(.param_bquote(index)))
  col_line <- "#023743FF"
  col_fill <- "#72874EFF"
  
  ggplot(estimates_data, aes_string(se_value)) +
    geom_histogram(fill = col_fill, colour = col_line, alpha = 0.7) +
    theme_minimal() +
    labs(x = l_xaxis, y = "Frequency", title = l_title)
}

#-----------------------------------------------------------------------#
#               SAVE PLOTS AS A GRID                                ----
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