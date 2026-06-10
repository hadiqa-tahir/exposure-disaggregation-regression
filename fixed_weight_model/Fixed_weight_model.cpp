#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  
  using namespace density;
  using namespace Eigen;
  
  // ======================================================================== //
  // DATA (passed in from R)
  // ======================================================================== //
  
  // Pixel-level exposure matrix (rows = pixels, cols = covariates).
  // Each row is one NDVI pixel belonging to one individual.
  // Long format: individual j occupies rows startindex(j) to
  // startindex(j) + n_pixels_vec(j) - 1.
  DATA_MATRIX(x);
  
  // Index of the first pixel row belonging to individual j (0-based).
  DATA_IVECTOR(startindex);
  
  // Number of pixels belonging to individual j.
  DATA_IVECTOR(n_pixels_vec);
  
  // Binary outcome for each individual: 1 = case, 0 = no case.
  // Length = number of individuals (n).
  DATA_VECTOR(cases);
  
  // Weight for each pixel row, pre-computed in R. Flexible to change. 
  // Example case: 
  //         Inverse distance weighting
  //         For individual j: w_i = (1/d_i) / sum_i(1/d_i), so sum over j's pixels = 1.
  DATA_VECTOR(weights);
  
  // ======================================================================== //
  // PARAMETERS (estimated by TMB)
  // ======================================================================== //
  
  // Logistic regression intercept (beta_0).
  PARAMETER(intercept);
  
  // Logistic regression slope(s) (beta_1, beta_2,...). 
  PARAMETER_VECTOR(slope);
  
  // ======================================================================== //
  // PRIORS (Passed from R)
  // ======================================================================== //
  
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  // ======================================================================== //
  // STUDY DIMENSIONS
  // ======================================================================== //
  
  int n    = cases.size();  // number of individuals
  int pixn = x.rows();      // total number of pixels across all individuals
  
  // ======================================================================== //
  // NEGATIVE LOG-LIKELIHOOD: PRIORS
  // ======================================================================== //
  
  Type nll = 0.0;
  
  // Normal priors on intercept and slope(s).
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  
  Type nll_priors = nll;
  REPORT(nll_priors);
  
  // ======================================================================== //
  // NEGATIVE LOG-LIKELIHOOD: DATA
  // ======================================================================== //
  
  // Step 1: compute the weighted linear predictor for every pixel.
  // pixel_linear_pred_i = w_i * (beta_0 + beta_1 * NDVI_i)
  vector<Type> pixel_linear_pred(pixn);
  pixel_linear_pred = weights * (intercept + x * slope);
  REPORT(pixel_linear_pred);
  
  // Step 2: for each individual j, sum pixel contributions to get the individual-level linear predictor:
  // eta_j = sum_{i in j} w_i * (beta_0 + beta_1 * NDVI_i)
  // P(ARI_j = 1) = invlogit(eta_j) = 1 / (1 + exp(-eta_j))
  vector<Type> pixel_linear_s;  // temporary pixel slice for individual j
  vector<Type> p(n);            // individual-level linear predictors (on logit scale)
  vector<Type> reportnll(n);
  
  for (int s = 0; s < n; s++) {
    
    // Extract pixel rows belonging to individual s.
    pixel_linear_s = pixel_linear_pred.segment(startindex(s), n_pixels_vec(s)).array();
    
    // Sum to get eta_j (logit-scale probability for individual s).
    p[s] = sum(pixel_linear_s);
    
    // Bernoulli log-likelihood. dbinom_robust takes the logit probability directly
    // No need to call invlogit manually.
    nll -= dbinom_robust(cases[s], Type(1), p[s], true);
  }
  
  Type nll_final = nll;
  REPORT(nll_final);
  REPORT(reportnll);
  REPORT(p);          // logit-scale linear predictors for each individual
  
  return nll;
}
