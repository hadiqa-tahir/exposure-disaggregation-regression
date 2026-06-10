# =============================================================================#
# generate_resampled_data.R
# =============================================================================#
#   This script bootstraps (samples with replacement) from the 1,114 synthetic 
#   clusters to produce data of 50,000 and 100,000 individuals.
#
#   Real DHS data cannot be shared publicly. To reproduce the paper's results,
#   apply for the Albania 2017-18 DHS data at https://dhsprogram.com
#
# =============================================================================#

library(arrow)
library(dplyr)
library(parallel)
library(here)

# =============================================================================#
# 1: Load "Observed" Data                                                 ----
# =============================================================================#

df <- arrow::read_parquet(here::here("data", "df_1114_synthetic.parquet"))
cat("Loaded original dataset:", nrow(df), "rows,", length(unique(df$HOUSE)), "individuals\n")

# =============================================================================#
# 2: Generate Dataset                                              ----
# =============================================================================#

set.seed(4000)
N_sample_size <- 50000            # or 100,000 - user input
sampled_houses <- data.frame(
  HOUSE_orig = sample(unique(df$HOUSE), N_sample_size, replace = TRUE),
  HOUSE = 1:N_sample_size
)

df_resampled <- sampled_houses %>%
  left_join(df %>% rename(HOUSE_orig = HOUSE), by = "HOUSE_orig") %>%
  select(-HOUSE_orig)
cat("dataset dimensions:", nrow(df_resampled), "rows,", length(unique(df_resampled$HOUSE)), "individuals\n")

# Save it
path <- paste0("df_", N_sample_size / 1000, "K_synthetic.parquet")
arrow::write_parquet(df_resampled, here::here("data", path))
cat("Saved:", path, "\n")
