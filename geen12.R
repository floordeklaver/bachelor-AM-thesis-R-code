# This file contains the R code used to investigate the instability observed in the pilot-plus-t-test simulation results.
# Written by Dr. Rianne de Heide

set.seed(42)

# Parameters
n        <- 12
sigma    <- 1
deltas   <- seq(0.2, 1.0, by = 0.1)
n_sims   <- 100000

# Function to simulate estimated effect size
simulate_d <- function(delta, n, sigma, n_sims) {
  mu <- delta * sigma
  
  d_hat <- replicate(n_sims, {
    x <- rnorm(n, mean = mu, sd = sigma)
    mean(x) / sd(x)
  })
  
  data.frame(
    delta_true = delta,
    d_hat      = d_hat
  )
}

# Run simulations
sim_results <- do.call(
  rbind,
  lapply(deltas, simulate_d, n = n, sigma = sigma, n_sims = n_sims)
)

# Summarize performance
summary_results <- aggregate(
  d_hat ~ delta_true,
  data = sim_results,
  FUN = function(x) c(
    mean = mean(x),
    bias = mean(x) - unique(sim_results$delta_true[sim_results$d_hat %in% x]),
    sd   = sd(x),
    rmse = sqrt(mean((x - unique(sim_results$delta_true[sim_results$d_hat %in% x]))^2))
  )
)

# Clean summary table
summary_table <- data.frame(
  delta_true = summary_results$delta_true,
  mean_d_hat = summary_results$d_hat[, "mean"],
  sd_d_hat   = summary_results$d_hat[, "sd"],
  rmse_d_hat = summary_results$d_hat[, "rmse"]
)

print(summary_table)


library(ggplot2)

ggplot(sim_results, aes(x = d_hat)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(aes(xintercept = delta_true),
             linetype = "dashed", color = "red") +
  facet_wrap(~ delta_true, scales = "free") +
  labs(
    title = "Sampling Distribution of estimated effect size (n = 12)",
    x = "Estimated effect size (d)",
    y = "Density"
  )

library(ggplot2)

# Select deltas for poster
poster_deltas <- c(0.2, 0.5, 1.0)

poster_data <- subset(sim_results, delta_true %in% poster_deltas)

ggplot(poster_data, aes(x = d_hat)) +
  geom_density(
    fill = "#4C72B0",
    alpha = 0.5,
    linewidth = 0.9
  ) +
  geom_vline(
    aes(xintercept = delta_true, color = "True effect", linetype = "True effect"),
    linewidth = 1.0
  ) +
  facet_wrap(~ delta_true, nrow = 1) +
  scale_color_manual(values = c("True effect" = "#C44E52"), name = NULL) +
  scale_linetype_manual(values = c("True effect" = "dashed"), name = NULL) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.0)),
    linetype = "none"   # keep a single legend (no duplicates)
  ) +
  labs(
    title = "Sampling distribution of estimated effect size",
    x = "Estimated effect size",
    y = "Density"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    strip.text = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    legend.position = c(0.98, 0.98),      # top-right inside
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.width = unit(1.6, "lines"),
    legend.text = element_text(size = 14)
  )


ggsave(
  "sampling_distribution_d_hat_poster_row.pdf",
  width = 12,
  height = 4.5,
  dpi = 300
)
