# This file contains the R code for the e-value-t-test simulations.
# Written by M.A.Y. de Klaver

system.time({
  ## ---- Step 0: Import ----
  
  # Load e-value functions from Ettest.R
  source("C:/Users/User/OneDrive - University of Twente/Documenten/UT/AM/2025-2026/B-ass _ e-values/RStudio/Ettest.R")
  
  # Load nmax values
  library(readxl)
  e_nmax_table <- read_excel("C:/Users/User/OneDrive - University of Twente/Documenten/UT/AM/2025-2026/B-ass _ e-values/e_nmax_table_kopie.xlsx")
  colnames(e_nmax_table) <- c("delta", "nmax_e1", "nmax_e2")
  
  
  ## ---- Step 1: Settings ----
  set.seed(42) # For reproducibility
  
  sigma <- 1 # true population SD
  mu0 <- 0 # null mean
  alpha <- 0.05
  power_target <- 0.80
  
  # Number of simulated pilots per true effect size
  n_rep <- 10000
  
  # Vector of true standardized effect sizes included in the table
  delta_values <- sort(e_nmax_table$delta)
  
  ## ---- Step 2: Prepare empty results table ----
  e_results <- data.frame(
    delta_true = delta_values,
    power_e1 = NA_real_, # power from e-value version 1 (estimated theta)
    mean_n_e1 = NA_real_, # mean sample size needed to cross 1/alpha for e1
    sd_n_e1 = NA_real_, # SD of sample size needed to cross 1/alpha for e1
    
    power_e2 = NA_real_, # power from e-value version 2 (oracle theta)
    mean_n_e2 = NA_real_, # mean sample size needed to cross 1/alpha for e2
    sd_n_e2 = NA_real_ # SD of sample size needed to cross 1/alpha for e2
  )
  
  
  ## ---- Step 3: Helper functions ----
  
  ## 1. retrive nmax from the table for a given delta and version
  get_nmax <- function(delta, version = c("e1", "e2")) {
    # "e1": nmax_e1  (estimated theta)
    # "e2": nmax_e2  (oracle theta)
    version <- match.arg(version)
    
    row <- e_nmax_table[e_nmax_table$delta == delta, ] # Find the row in the table that matches this delta
    
    if (nrow(row) == 0) {
      stop("No nmax found in e_nmax_table for delta = ", delta)
    }
    
    # Return the appropriate nmax
    if (version == "e1") return(row$nmax_e1)
    if (version == "e2") return(row$nmax_e2)
  }
  
  ## 2. simulate one sequential e-value study for both e1 and e2
  simulate_one_e_trial_both <- function(delta, alpha = 0.05, mu0 = 0, sigma = 1) {
    
    # Design max sample sizes (from the table)
    nmax_e1 <- get_nmax(delta, "e1")
    nmax_e2 <- get_nmax(delta, "e2")
    nmax_total <- max(nmax_e1, nmax_e2)
    
    # Generate one full data path (same data used for both tests)
    mu_true <- mu0 + delta * sigma
    x <- rnorm(nmax_total, mean = mu_true, sd = sigma)
    
    thresh <- 1 / alpha # e-value decision threshold
    
    # Defaults: no rejection within design nmax
    reject_e1 <- FALSE
    reject_e2 <- FALSE
    n_stop_e1 <- nmax_e1
    n_stop_e2 <- nmax_e2
    
    # Check e-values sequentially
    for (k in 2:nmax_total) {
      xk <- x[1:k]
      
      if (!reject_e1 && k <= nmax_e1) {
        e1 <- e.t.test1(xk)
        if (e1 > thresh) {
          reject_e1 <- TRUE
          n_stop_e1 <- k
        }
      }
      
      if (!reject_e2 && k <= nmax_e2) {
        e2 <- e.t.test2(xk, thetaR = delta)
        if (e2 > thresh) {
          reject_e2 <- TRUE
          n_stop_e2 <- k
        }
      }
      
      # Stop early if both methods are done (rejected or reached their nmax)
      if ((reject_e1 || k > nmax_e1) && (reject_e2 || k > nmax_e2)) break
    }
    
    list(
      reject_e1 = reject_e1,
      n_stop_e1 = n_stop_e1,
      reject_e2 = reject_e2,
      n_stop_e2 = n_stop_e2
    )
  }
  
  
  ## ---- Step 4: Loop over all deltas and repetitions ----
  for (i in seq_along(delta_values)) {
    
    delta_true <- delta_values[i]
    
    # Storage for this delta over n_rep simulated studies
    n_stops_e1 <- numeric(n_rep)
    n_stops_e2 <- numeric(n_rep)
    rejects_e1 <- logical(n_rep)
    rejects_e2 <- logical(n_rep)
    
    for (r in seq_len(n_rep)) {
      
      out <- simulate_one_e_trial_both(delta = delta_true,
                                       alpha = alpha,
                                       mu0 = mu0,
                                       sigma = sigma)
      
      # Store outcomes for this repetition
      n_stops_e1[r] <- out$n_stop_e1
      n_stops_e2[r] <- out$n_stop_e2
      rejects_e1[r] <- out$reject_e1
      rejects_e2[r] <- out$reject_e2
    }
    
    # Summaries for this true delta
    e_results$power_e1[i] <- mean(rejects_e1)
    e_results$mean_n_e1[i] <- mean(n_stops_e1)
    e_results$sd_n_e1[i] <- sd(n_stops_e1)
    
    e_results$power_e2[i] <- mean(rejects_e2)
    e_results$mean_n_e2[i] <- mean(n_stops_e2)
    e_results$sd_n_e2[i] <- sd(n_stops_e2)
    
  }
}) # End system.time

e_results # print results


## ---- Plots for e-results ----

# y-limits that work for both methods (ignoring NA)
y_max <- max(
  e_results$mean_n_e1 + e_results$sd_n_e1,
  e_results$mean_n_e2 + e_results$sd_n_e2,
  na.rm = TRUE
)
y_min <- min(
  e_results$mean_n_e1 - e_results$sd_n_e1,
  e_results$mean_n_e2 - e_results$sd_n_e2,
  na.rm = TRUE
)
y_min <- max(0, y_min)

y_lim_global <- range(
  c(
    e_results$mean_n_e1 - e_results$sd_n_e1,
    e_results$mean_n_e1 + e_results$sd_n_e1,
    e_results$mean_n_e2 - e_results$sd_n_e2,
    e_results$mean_n_e2 + e_results$sd_n_e2,
    e_nmax_table$nmax_e1,
    e_nmax_table$nmax_e2
  ),
  na.rm = TRUE
)

y_lim_global[1] <- max(0, y_lim_global[1])

# limits for the zoom-ins
idx_zoom <- e_results$delta_true >= 0.5 & e_results$delta_true <= 1.0

y_lim_zoom <- range(
  c(
    e_results$mean_n_e1[idx_zoom] - e_results$sd_n_e1[idx_zoom],
    e_results$mean_n_e1[idx_zoom] + e_results$sd_n_e1[idx_zoom],
    e_results$mean_n_e2[idx_zoom] - e_results$sd_n_e2[idx_zoom],
    e_results$mean_n_e2[idx_zoom] + e_results$sd_n_e2[idx_zoom]
  ),
  na.rm = TRUE
)

y_lim_zoom[1] <- max(0, y_lim_zoom[1])


  ## ---- Plot 1: Statistical plots (points + error bars) ----
  
  ## 1) e1: full range (0.2–1.0)
  idx_e1_full <- e_results$delta_true >= 0.2 & e_results$delta_true <= 1.0
  
  plot(
    e_results$delta_true[idx_e1_full],
    e_results$mean_n_e1[idx_e1_full],
    type = "p",
    pch = 16,
    xlab = expression(delta~"(true)"),
    ylab = "Sample size at stop (mean ± SD)",
    ylim = y_lim_global,
    main = "e1: Stopping sample size vs true delta"
  )
  
  arrows(
    x0 = e_results$delta_true[idx_e1_full],
    y0 = e_results$mean_n_e1[idx_e1_full] - e_results$sd_n_e1[idx_e1_full],
    x1 = e_results$delta_true[idx_e1_full],
    y1 = e_results$mean_n_e1[idx_e1_full] + e_results$sd_n_e1[idx_e1_full],
    angle = 90, code = 3, length = 0.05
  )
  
  
  ## 2) e1: zoom range (0.5–1.0)
  idx_e1_zoom <- e_results$delta_true >= 0.5 & e_results$delta_true <= 1.0
  
  plot(
    e_results$delta_true[idx_zoom],
    e_results$mean_n_e1[idx_zoom],
    type = "p",
    pch = 16,
    ylim = y_lim_zoom
  )
  
  arrows(
    x0 = e_results$delta_true[idx_zoom],
    y0 = e_results$mean_n_e1[idx_zoom] - e_results$sd_n_e1[idx_zoom],
    x1 = e_results$delta_true[idx_zoom],
    y1 = e_results$mean_n_e1[idx_zoom] + e_results$sd_n_e1[idx_zoom],
    angle = 90, code = 3, length = 0.05
  )
  
  ## 3) e2: full range (0.2–1.0)
  idx_e2_full <- e_results$delta_true >= 0.2 & e_results$delta_true <= 1.0
  
  plot(
    e_results$delta_true[idx_e2_full],
    e_results$mean_n_e2[idx_e2_full],
    type = "p",
    pch = 16,
    xlab = expression(delta~"(true)"),
    ylab = "Sample size at stop (mean ± SD)",
    ylim = y_lim_global,
    main = "e2: stopping sample size vs true delta"
  )
  
  arrows(
    x0 = e_results$delta_true[idx_e2_full],
    y0 = e_results$mean_n_e2[idx_e2_full] - e_results$sd_n_e2[idx_e2_full],
    x1 = e_results$delta_true[idx_e2_full],
    y1 = e_results$mean_n_e2[idx_e2_full] + e_results$sd_n_e2[idx_e2_full],
    angle = 90, code = 3, length = 0.05
  )
  
  
  ## 4) e2: zoom range (0.5–1.0)
  idx_e2_zoom <- e_results$delta_true >= 0.5 & e_results$delta_true <= 1.0
  plot(
    e_results$delta_true[idx_zoom],
    e_results$mean_n_e2[idx_zoom],
    type = "p",
    pch = 16,
    ylim = y_lim_zoom
  )
  
  arrows(
    x0 = e_results$delta_true[idx_zoom],
    y0 = e_results$mean_n_e2[idx_zoom] - e_results$sd_n_e2[idx_zoom],
    x1 = e_results$delta_true[idx_zoom],
    y1 = e_results$mean_n_e2[idx_zoom] + e_results$sd_n_e2[idx_zoom],
    angle = 90, code = 3, length = 0.05
  )
  
  
  
  ## ---- Plot 2: Simplified plots (lines through points + nmax) ----
  
  # colours (grayscale-safe)
  col_e1 <- "#1B4F72"      # dark blue
  col_e2 <- "#7D2E2E"      # dark red
  col_e1_nmax <- "#5DADE2" # light blue
  col_e2_nmax <- "#C0392B" # light red
  
  plot(
    e_results$delta_true,
    e_results$mean_n_e1,
    type = "b",
    pch = 16,
    col = col_e1,
    xlab = expression(delta~"(true)"),
    ylab = "Sample size",
    ylim = y_lim_global,
    main = "E-value stopping sample size vs true delta"
  )
  
  # mean stopping n (e2)
  lines(e_results$delta_true, e_results$mean_n_e2, type = "b", pch = 17, col = col_e2)
  
  # nmax lines (clearly distinct, both non-solid)
  lines(e_nmax_table$delta, e_nmax_table$nmax_e1, lty = 2, lwd = 2, col = col_e1_nmax)  # dashed
  lines(e_nmax_table$delta, e_nmax_table$nmax_e2, lty = 3, lwd = 2, col = col_e2_nmax)  # dot
  
  legend(
    "topright",
    legend = c(
      "e1 mean stop",
      "e2 mean stop",
      "e1 nmax",
      "e2 nmax"
    ),
    pch = c(16, 17, NA, NA),
    lty = c(1, 1, 2, 3),
    lwd = c(1, 1, 2, 2),
    col = c(col_e1, col_e2, col_e1_nmax, col_e2_nmax),
    bty = "n"
  )
  
  ## ---- Extra plot for poster: e1 statistics + nmax ----
  
  idx_e1_full <- e_results$delta_true >= 0.2 & e_results$delta_true <= 1.0
  
  plot(
    e_results$delta_true[idx_e1_full],
    e_results$mean_n_e1[idx_e1_full],
    type = "p",
    pch = 16,
    col = col_e1,
    xlab = expression(delta~"(true)"),
    ylab = "Sample size at stop",
    ylim = y_lim_global,
    main = "Stopping sample size (mean ± SD) with nmax"
  )
  
  # error bars (mean ± SD)
  arrows(
    x0 = e_results$delta_true[idx_e1_full],
    y0 = e_results$mean_n_e1[idx_e1_full] - e_results$sd_n_e1[idx_e1_full],
    x1 = e_results$delta_true[idx_e1_full],
    y1 = e_results$mean_n_e1[idx_e1_full] + e_results$sd_n_e1[idx_e1_full],
    angle = 90, code = 3, length = 0.05,
    col = col_e1
  )
  
  # nmax line for e1
  lines(
    e_nmax_table$delta,
    e_nmax_table$nmax_e1,
    lty = 2,
    lwd = 2,
    col = col_e1_nmax
  )
  
  legend(
    "topright",
    legend = c("mean stop", "mean ± SD", "nmax"),
    pch = c(16, NA, NA),
    lty = c(NA, 1, 2),
    lwd = c(NA, 1, 2),
    col = c(col_e1, col_e1, col_e1_nmax),
    bty = "n"
  )
