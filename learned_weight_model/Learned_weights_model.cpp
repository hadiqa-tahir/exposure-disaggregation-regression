// Learned-Weight (LW) model.
// Extends the Fixed-Weight model by estimating the distance-decay exponent
// alpha from the data rather than fixing it at 1.
// Weights: w_i = (1/d_i^alpha) / sum(1/d_i^alpha)
// Linear predictor: n_j = sum_i [ w_i * (intercept + slope * NDVI_i) ]
// P(ARI_j) = logit^-1(n_j)

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  
  using namespace density;
  using namespace Eigen;
  
  // ------------------------------------------------------------------------ //
  // Data
  // ------------------------------------------------------------------------ //
  
  // Environmental/covariate data matrix — one row per pixel, in long format.
  // If a pixel falls within multiple individuals' buffers it appears multiple times.
  DATA_MATRIX(x);
  
  // Start index and length for each individual's pixel segment.
  // Used to slice the long-format pixel vectors back into per-individual chunks.
  DATA_IVECTOR(startindex);
  DATA_IVECTOR(n_pixels_vec);
  
  // Binary outcome (ARI = 1 / 0), one entry per individual.
  DATA_VECTOR(cases);
  
  // Euclidean distances from each pixel centroid to its individual's cluster centroid (metres).
  // Passed raw — weights are computed inside the model using the learned alpha.
  DATA_VECTOR(distances);
  
  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //
  
  PARAMETER(intercept);          // b0: overall intercept
  PARAMETER_VECTOR(slope);       // b1: NDVI slope (vector to generalise to multiple covariates)
  PARAMETER(log_alpha);          // log(alpha): log of the distance-decay exponent.
  // Estimated on the log scale to ensure alpha > 0.
  
  // Priors: intercept and slope ~ N(0, 10^2)
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  // Prior: log_alpha ~ N(priormean_log_alpha, priorsd_log_alpha^2)
  // Default: N(-1, 1^2) — centres alpha near exp(-1) ~ 0.37.
  DATA_SCALAR(priormean_log_alpha);
  DATA_SCALAR(priorsd_log_alpha);
  
  // Total number of individuals and pixels
  int n    = cases.size();
  int pixn = x.rows();
  
  // Back-transform alpha from log scale
  Type alpha = exp(log_alpha);
  REPORT(alpha);
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors
  // ------------------------------------------------------------------------ //
  
  // Initialise negative log-likelihood
  Type nll = 0.0;
  
  // Prior contributions for intercept, slope(s), and log_alpha
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for(int s = 0; s < slope.size(); s++){
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  nll -= dnorm(log_alpha, priormean_log_alpha, priorsd_log_alpha, true);
  
  // Jacobian adjustment for the log transform: d(alpha)/d(log_alpha) = alpha = exp(log_alpha).
  // Adding -log_alpha corrects the density from the log scale back to the original scale.
  nll -= log_alpha;
  
  Type nll_priors = nll;
  REPORT(nll_priors);
  
  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //
  
  // Compute unnormalised weights: u_i = 1 / d_i^alpha
  vector<Type> u(pixn);
  u = 1 / pow(distances, alpha);
  
  // Small constant to prevent division by zero in edge cases
  Type eps = 1e-50;
  
  // Normalise within each individual's buffer: w_i = u_i / sum_i(u_i)
  // so that weights sum to 1 for each individual.
  vector<Type> weights(pixn);
  for(int s = 0; s < n; s++) {
    vector<Type> u_segment = u.segment(startindex(s), n_pixels_vec(s));
    weights.segment(startindex(s), n_pixels_vec(s)) = u_segment / (u_segment.sum() + eps);
  }
  
  // Pixel-level linear predictor: w_i * (intercept + slope * NDVI_i)
  vector<Type> pixel_linear_pred(pixn);
  pixel_linear_pred = weights * (intercept + x * slope);
  
  REPORT(weights);
  REPORT(pixel_linear_pred);
  
  // Aggregate to individual level: n_j = sum_i [ w_i * (intercept + slope * NDVI_i) ]
  // then evaluate Bernoulli likelihood via dbinom_robust (takes logit-scale probability).
  vector<Type> pixel_linear_s;
  vector<Type> p(n);
  vector<Type> reportnll(n);
  
  for (int s = 0; s < n; s++) {
    // Extract this individual's pixel linear predictors
    pixel_linear_s = pixel_linear_pred.segment(startindex(s), n_pixels_vec(s)).array();
    
    // Sum to get the individual-level linear predictor (on logit scale)
    p[s] = sum(pixel_linear_s);
    
    // dbinom_robust takes logit(p) directly — no manual invlogit needed
    nll -= dbinom_robust(cases[s], Type(1), p[s], true);
  }
  
  Type nll_final = nll;
  REPORT(nll_final);
  REPORT(reportnll);
  REPORT(p);
  
  return nll;
}
