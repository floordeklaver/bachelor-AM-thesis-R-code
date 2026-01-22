# This file contains the R code for the pilot-plus-t-test simulation. It includes the five different ways the SD was calculated.
# Written by M.A.Y. de Klaver 
# SD variation functions provided by Dr. Rianne de Heide

system.time({
  
  ## ---- Step 0: Settings ----
  set.seed(42) # For reproducibility
  
  np <- 12 # pilot sample size
  sigma <- 1 # true population SD
  mu0 <- 0 # null mean
  
  # Number of simulated pilots per true effect size
  n_rep <- 10000
  
  # Vector of true standardized effect sizes (including 0)
  delta_values <- seq(0.2, 1.0, by = 0.1)
  
  ## ---- Step 1: Prepare empty results table ----
  pilot_results <- data.frame(
    delta_true = delta_values,
    
    # raw SD results
    delta_hat_mean_raw = NA_real_,
    delta_hat_sd_raw = NA_real_,
    n_z_based_mean_raw = NA_real_,
    n_z_based_sd_raw = NA_real_,
    n_power_t_mean_raw = NA_real_,
    n_power_t_sd_raw = NA_real_,
    s_raw_mean = NA_real_,
    s_raw_sd = NA_real_,
    
    # SD = 1
    delta_hat_mean_s1 = NA_real_,
    delta_hat_sd_s1 = NA_real_,
    n_z_based_mean_s1 = NA_real_,
    n_z_based_sd_s1 = NA_real_,
    n_power_t_mean_s1 = NA_real_,
    n_power_t_sd_s1 = NA_real_,
    s_s1_mean = NA_real_,
    s_s1_sd = NA_real_,
    
    # bias-corrected SD
    delta_hat_mean_corr = NA_real_,
    delta_hat_sd_corr = NA_real_,
    n_z_based_mean_corr = NA_real_,
    n_z_based_sd_corr = NA_real_,
    n_power_t_mean_corr = NA_real_,
    n_power_t_sd_corr = NA_real_,
    s_corr_mean = NA_real_,
    s_corr_sd = NA_real_,
    
    # stabilized SD
    delta_hat_mean_stab = NA_real_,
    delta_hat_sd_stab = NA_real_,
    n_z_based_mean_stab = NA_real_,
    n_z_based_sd_stab = NA_real_,
    n_power_t_mean_stab = NA_real_,
    n_power_t_sd_stab = NA_real_,
    s_stab_mean = NA_real_,
    s_stab_sd = NA_real_,
    
    # Stein shrinkage
    delta_hat_mean_shrink = NA_real_,
    delta_hat_sd_shrink = NA_real_,
    n_z_based_mean_shrink = NA_real_,
    n_z_based_sd_shrink = NA_real_,
    n_power_t_mean_shrink = NA_real_,
    n_power_t_sd_shrink = NA_real_,
    s_shrink_mean = NA_real_,
    s_shrink_sd = NA_real_
  )
  
  ## ---- Step 2: Sample size formulas ----
  alpha <- 0.05
  power_target <- 0.80
  
  ## 1. z-based formula (normal approximation)
  n_z_based <- function(delta, alpha = 0.05, power = 0.80) {
    d <- abs(delta)
    if (is.na(d) || d == 0) return(NA_real_)
    
    z_alpha <- qnorm(1 - alpha / 2)
    z_beta <- qnorm(power)
    
    ((z_alpha + z_beta) / d)^2
  }
  
  ## 2. R's power.t.test function
  n_power_t <- function(delta, alpha = 0.05, power = 0.80) {
    d <- abs(delta)
    if (is.na(d) || d == 0) return(NA_real_)
    
    power.t.test(delta = d,
                 sd = 1, # delta is standardised
                 sig.level = alpha,
                 power = power,
                 type = "one.sample",
                 alternative = "two.sided")$n
  }
  
  ## ---- SD correction helper functions ----
  c4 <- function(n) sqrt(2/(n-1)) * gamma(n/2) / gamma((n-1)/2)
  sd_unbiased <- function(x) c4(length(x)) * sd(x) # A correction for the bias: multiply by c(n)
  
  sd_stab <- function(s, n, tau2 = 1/n) sqrt(s^2 + tau2) # Shrinkage with a “stabilized SD”
  
  shrink_sigma <- function(s, n, sigma0 = 1) {
    w <- max(0, (n - 3) / n)
    w * s + (1 - w) * sigma0
  } # Stein Shrinkage
  
  ## ---- Step 3: Loop over all true effect sizes and repetitions ----
  for (i in seq_along(delta_values)) {
    
    delta_true <- delta_values[i]
    
    # Storage for this delta over n_rep simulated pilots
    dhat_raw <- numeric(n_rep)
    dhat_s1 <- numeric(n_rep)
    dhat_corr <- numeric(n_rep)
    dhat_stab <- numeric(n_rep)
    dhat_shrink <- numeric(n_rep)
    
    n_z_raw <- numeric(n_rep)
    n_z_s1 <- numeric(n_rep)
    n_z_corr <- numeric(n_rep)
    n_z_stab <- numeric(n_rep)
    n_z_shrink <- numeric(n_rep)
    
    n_pt_raw <- numeric(n_rep)
    n_pt_s1 <- numeric(n_rep)
    n_pt_corr <- numeric(n_rep)
    n_pt_stab <- numeric(n_rep)
    n_pt_shrink <- numeric(n_rep)
    
    s_raw_vec <- numeric(n_rep)
    s_s1_vec <- numeric(n_rep)
    s_corr_vec <- numeric(n_rep)
    s_stab_vec <- numeric(n_rep)
    s_shrink_vec <- numeric(n_rep)
    
    for (r in seq_len(n_rep)) {
      
      # Simulate the pilot data for this repetition
      x_pilot <- rnorm(np, mean = mu0 + delta_true, sd = sigma)
      
      # Sample mean and sample SD from the pilot
      m <- mean(x_pilot)
      s_raw <- sd(x_pilot)
      
      # All SD alternatives
      s1 <- 1
      s_corr <- sd_unbiased(x_pilot)
      s_stab <- sd_stab(s_raw, np)
      s_sh <- shrink_sigma(s_raw, np, sigma0 = 1)
      
      # Store SDs
      s_raw_vec[r] <- s_raw
      s_s1_vec[r] <- s1
      s_corr_vec[r] <- s_corr
      s_stab_vec[r] <- s_stab
      s_shrink_vec[r] <- s_sh
      
      # Estimated standardized effects
      dhat_raw[r] <- (m - mu0) / s_raw
      dhat_s1[r] <- (m - mu0) / s1
      dhat_corr[r] <- (m - mu0) / s_corr
      dhat_stab[r] <- (m - mu0) / s_stab
      dhat_shrink[r] <- (m - mu0) / s_sh
      
      # Required sample sizes
      n_z_raw[r] <- n_z_based(dhat_raw[r])
      n_z_s1[r] <- n_z_based(dhat_s1[r])
      n_z_corr[r] <- n_z_based(dhat_corr[r])
      n_z_stab[r] <- n_z_based(dhat_stab[r])
      n_z_shrink[r] <- n_z_based(dhat_shrink[r])
      
      n_pt_raw[r] <- n_power_t(dhat_raw[r])
      n_pt_s1[r] <- n_power_t(dhat_s1[r])
      n_pt_corr[r] <- n_power_t(dhat_corr[r])
      n_pt_stab[r] <- n_power_t(dhat_stab[r])
      n_pt_shrink[r] <- n_power_t(dhat_shrink[r])
    }
    
    ## ---- Step 4: Summaries ----
    
    # raw
    pilot_results$delta_hat_mean_raw[i] <- mean(dhat_raw, na.rm = TRUE)
    pilot_results$delta_hat_sd_raw[i] <- sd(dhat_raw, na.rm = TRUE)
    pilot_results$s_raw_mean[i] <- mean(s_raw_vec, na.rm = TRUE)
    pilot_results$s_raw_sd[i] <- sd(s_raw_vec, na.rm = TRUE)
    pilot_results$n_z_based_mean_raw[i] <- mean(n_z_raw, na.rm = TRUE)
    pilot_results$n_z_based_sd_raw[i] <- sd(n_z_raw, na.rm = TRUE)
    pilot_results$n_power_t_mean_raw[i] <- mean(n_pt_raw, na.rm = TRUE)
    pilot_results$n_power_t_sd_raw[i] <- sd(n_pt_raw, na.rm = TRUE)
    
    # SD = 1
    pilot_results$delta_hat_mean_s1[i] <- mean(dhat_s1, na.rm = TRUE)
    pilot_results$delta_hat_sd_s1[i] <- sd(dhat_s1, na.rm = TRUE)
    pilot_results$s_s1_mean[i] <- mean(s_s1_vec, na.rm = TRUE)
    pilot_results$s_s1_sd[i] <- sd(s_s1_vec, na.rm = TRUE)
    pilot_results$n_z_based_mean_s1[i] <- mean(n_z_s1, na.rm = TRUE)
    pilot_results$n_z_based_sd_s1[i] <- sd(n_z_s1, na.rm = TRUE)
    pilot_results$n_power_t_mean_s1[i] <- mean(n_pt_s1, na.rm = TRUE)
    pilot_results$n_power_t_sd_s1[i] <- sd(n_pt_s1, na.rm = TRUE)
    
    # bias corrected
    pilot_results$delta_hat_mean_corr[i] <- mean(dhat_corr, na.rm = TRUE)
    pilot_results$delta_hat_sd_corr[i] <- sd(dhat_corr, na.rm = TRUE)
    pilot_results$s_corr_mean[i] <- mean(s_corr_vec, na.rm = TRUE)
    pilot_results$s_corr_sd[i] <- sd(s_corr_vec, na.rm = TRUE)
    pilot_results$n_z_based_mean_corr[i] <- mean(n_z_corr, na.rm = TRUE)
    pilot_results$n_z_based_sd_corr[i] <- sd(n_z_corr, na.rm = TRUE)
    pilot_results$n_power_t_mean_corr[i] <- mean(n_pt_corr, na.rm = TRUE)
    pilot_results$n_power_t_sd_corr[i] <- sd(n_pt_corr, na.rm = TRUE)
    
    # stabilized
    pilot_results$delta_hat_mean_stab[i] <- mean(dhat_stab, na.rm = TRUE)
    pilot_results$delta_hat_sd_stab[i] <- sd(dhat_stab, na.rm = TRUE)
    pilot_results$s_stab_mean[i] <- mean(s_stab_vec, na.rm = TRUE)
    pilot_results$s_stab_sd[i] <- sd(s_stab_vec, na.rm = TRUE)
    pilot_results$n_z_based_mean_stab[i] <- mean(n_z_stab, na.rm = TRUE)
    pilot_results$n_z_based_sd_stab[i] <- sd(n_z_stab, na.rm = TRUE)
    pilot_results$n_power_t_mean_stab[i] <- mean(n_pt_stab, na.rm = TRUE)
    pilot_results$n_power_t_sd_stab[i] <- sd(n_pt_stab, na.rm = TRUE)
    
    # shrinkage
    pilot_results$delta_hat_mean_shrink[i] <- mean(dhat_shrink, na.rm = TRUE)
    pilot_results$delta_hat_sd_shrink[i] <- sd(dhat_shrink, na.rm = TRUE)
    pilot_results$s_shrink_mean[i] <- mean(s_shrink_vec, na.rm = TRUE)
    pilot_results$s_shrink_sd[i] <- sd(s_shrink_vec, na.rm = TRUE)
    pilot_results$n_z_based_mean_shrink[i] <- mean(n_z_shrink, na.rm = TRUE)
    pilot_results$n_z_based_sd_shrink[i] <- sd(n_z_shrink, na.rm = TRUE)
    pilot_results$n_power_t_mean_shrink[i] <- mean(n_pt_shrink, na.rm = TRUE)
    pilot_results$n_power_t_sd_shrink[i] <- sd(n_pt_shrink, na.rm = TRUE)
  }
  
  ## ---- Step 5: Add total sample sizes (pilot + planned) ----
  pilot_results$n_z_based_total_mean_raw <- ceiling(pilot_results$n_z_based_mean_raw + np)
  pilot_results$n_z_based_total_mean_s1 <- ceiling(pilot_results$n_z_based_mean_s1 + np)
  pilot_results$n_z_based_total_mean_corr <- ceiling(pilot_results$n_z_based_mean_corr + np)
  pilot_results$n_z_based_total_mean_stab <- ceiling(pilot_results$n_z_based_mean_stab + np)
  pilot_results$n_z_based_total_mean_shrink <- ceiling(pilot_results$n_z_based_mean_shrink + np)
  
  pilot_results$n_power_t_total_mean_raw <- ceiling(pilot_results$n_power_t_mean_raw + np)
  pilot_results$n_power_t_total_mean_s1 <- ceiling(pilot_results$n_power_t_mean_s1 + np)
  pilot_results$n_power_t_total_mean_corr <- ceiling(pilot_results$n_power_t_mean_corr + np)
  pilot_results$n_power_t_total_mean_stab <- ceiling(pilot_results$n_power_t_mean_stab + np)
  pilot_results$n_power_t_total_mean_shrink <- ceiling(pilot_results$n_power_t_mean_shrink + np)
  
  
## ---- Tables per method ----

# 1) raw SD (original behaviour)
table_raw <- data.frame(
  delta_true = pilot_results$delta_true,
  delta_hat_mean = pilot_results$delta_hat_mean_raw,
  delta_hat_sd = pilot_results$delta_hat_sd_raw,
  s_mean = pilot_results$s_raw_mean,
  s_sd = pilot_results$s_raw_sd,
  n_z_based_mean = pilot_results$n_z_based_mean_raw,
  n_z_based_sd = pilot_results$n_z_based_sd_raw,
  n_power_t_mean = pilot_results$n_power_t_mean_raw,
  n_power_t_sd = pilot_results$n_power_t_sd_raw,
  n_z_based_total_mean = pilot_results$n_z_based_total_mean_raw,
  n_power_t_total_mean = pilot_results$n_power_t_total_mean_raw
)

# 2) SD = 1
table_s1 <- data.frame(
  delta_true = pilot_results$delta_true,
  delta_hat_mean = pilot_results$delta_hat_mean_s1,
  delta_hat_sd = pilot_results$delta_hat_sd_s1,
  s_mean = pilot_results$s_s1_mean,
  s_sd = pilot_results$s_s1_sd,
  n_z_based_mean = pilot_results$n_z_based_mean_s1,
  n_z_based_sd = pilot_results$n_z_based_sd_s1,
  n_power_t_mean = pilot_results$n_power_t_mean_s1,
  n_power_t_sd = pilot_results$n_power_t_sd_s1,
  n_z_based_total_mean = pilot_results$n_z_based_total_mean_s1,
  n_power_t_total_mean = pilot_results$n_power_t_total_mean_s1
)

# 3) bias-corrected SD
table_corr <- data.frame(
  delta_true = pilot_results$delta_true,
  delta_hat_mean = pilot_results$delta_hat_mean_corr,
  delta_hat_sd = pilot_results$delta_hat_sd_corr,
  s_mean = pilot_results$s_corr_mean,
  s_sd = pilot_results$s_corr_sd,
  n_z_based_mean = pilot_results$n_z_based_mean_corr,
  n_z_based_sd = pilot_results$n_z_based_sd_corr,
  n_power_t_mean = pilot_results$n_power_t_mean_corr,
  n_power_t_sd = pilot_results$n_power_t_sd_corr,
  n_z_based_total_mean = pilot_results$n_z_based_total_mean_corr,
  n_power_t_total_mean = pilot_results$n_power_t_total_mean_corr
)

# 4) stabilized SD
table_stab <- data.frame(
  delta_true = pilot_results$delta_true,
  delta_hat_mean = pilot_results$delta_hat_mean_stab,
  delta_hat_sd = pilot_results$delta_hat_sd_stab,
  s_mean = pilot_results$s_stab_mean,
  s_sd = pilot_results$s_stab_sd,
  n_z_based_mean = pilot_results$n_z_based_mean_stab,
  n_z_based_sd = pilot_results$n_z_based_sd_stab,
  n_power_t_mean = pilot_results$n_power_t_mean_stab,
  n_power_t_sd = pilot_results$n_power_t_sd_stab,
  n_z_based_total_mean = pilot_results$n_z_based_total_mean_stab,
  n_power_t_total_mean = pilot_results$n_power_t_total_mean_stab
)

# 5) Stein shrinkage
table_shrink <- data.frame(
  delta_true = pilot_results$delta_true,
  delta_hat_mean = pilot_results$delta_hat_mean_shrink,
  delta_hat_sd = pilot_results$delta_hat_sd_shrink,
  s_mean = pilot_results$s_shrink_mean,
  s_sd = pilot_results$s_shrink_sd,
  n_z_based_mean = pilot_results$n_z_based_mean_shrink,
  n_z_based_sd = pilot_results$n_z_based_sd_shrink,
  n_power_t_mean = pilot_results$n_power_t_mean_shrink,
  n_power_t_sd = pilot_results$n_power_t_sd_shrink,
  n_z_based_total_mean = pilot_results$n_z_based_total_mean_shrink,
  n_power_t_total_mean = pilot_results$n_power_t_total_mean_shrink
)

}) # End system.time

# Print the different tables
table_raw
table_s1
table_corr
table_stab
table_shrink