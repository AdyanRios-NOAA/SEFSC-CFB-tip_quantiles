# ==============================================================================
# LOBSTER CARAPACE LENGTH: RAPID INDICATOR & TRIAGE SCRIPT
# ==============================================================================
# Prerequisites: 
# - Ensure your data frame is named `lobster_data`
# - Must contain columns `length_mm` and `year`
# ==============================================================================

# Install missing packages if needed:
# install.packages(c("dplyr", "ggplot2", "changepoint", "purrr", "here"))

# 1. Load Required Packages
library(dplyr)
library(ggplot2)
library(changepoint)
library(purrr)

# Read in TIP data
lobster_data <- read.csv(
  here::here("data", "pr_carspinylob_tip_8124_20251209_c.csv")
) |>
  dplyr::filter(year >= 1985)

# 2. Calculate Annual Medians and Bootstrapped Precision
annual_indicators <- lobster_data |>
  filter(!is.na(length_mm), !is.na(year)) |>
  group_by(year) |>
  summarise(
    n = n(),
    median_length = median(length_mm),
    # Non-parametric bootstrap of the median's sampling distribution (1,000 iterations)
    se_median = sd(replicate(1000, median(sample(length_mm, replace = TRUE)))),
    .groups = "drop"
  ) |>
  mutate(
    # Year-specific standard error envelope (+/- 1 SE)
    lower_se = median_length - se_median,
    upper_se = median_length + se_median
  ) |>
  arrange(year)

# 3. Changepoint Analysis (Detecting Regimes)
# PELT algorithm detects structural mean shifts in the time series
cpt_results <- cpt.mean(
  annual_indicators$median_length, 
  method = "PELT", 
  penalty = "SIC"
)
break_indices <- cpts(cpt_results)
break_years <- annual_indicators$year[break_indices]
cat("Detected structural breakpoints occurring after years:", break_years, "\n\n")

# 4. Establish Reference Baseline & Robust Monitoring Bounds (2007 - Present)
baseline_years <- annual_indicators |>
  filter(year >= 2007) 

ref_median <- median(baseline_years$median_length)

# Calculate robust scale, falling back to NIQR or SD if MAD collapses to 0 due to integer ties
raw_mad <- median(abs(baseline_years$median_length - ref_median))

if (raw_mad > 0) {
  ref_dispersion <- 1.4826 * raw_mad
  metric_used <- "NMAD"
} else {
  niqr <- IQR(baseline_years$median_length) / 1.349
  if (niqr > 0) {
    ref_dispersion <- niqr
    metric_used <- "Normalized IQR"
  } else {
    ref_dispersion <- sd(baseline_years$median_length)
    metric_used <- "Inter-annual SD"
  }
}

thresholds <- list(
  upper_warn   = ref_median + (1.0 * ref_dispersion),
  lower_warn   = ref_median - (1.0 * ref_dispersion),
  upper_action = ref_median + (2.0 * ref_dispersion),
  lower_action = ref_median - (2.0 * ref_dispersion)
)

cat(paste0("--- Baseline Parameters (2007-Present) ---\n",
           "Metric Used:           ", metric_used, "\n",
           "Reference Median (M0): ", round(ref_median, 2), " mm\n",
           "Inter-annual Spread:   ", round(ref_dispersion, 2), " mm\n",
           "+/- 1 Spread Warning:  [", round(thresholds$lower_warn, 2), ", ", round(thresholds$upper_warn, 2), "] mm\n",
           "+/- 2 Spread Action:   [", round(thresholds$lower_action, 2), ", ", round(thresholds$upper_action, 2), "] mm\n\n"))

# 5. Rapid Monitoring Classification
annual_monitored <- annual_indicators |>
  mutate(
    status = case_when(
      median_length > thresholds$upper_action | median_length < thresholds$lower_action ~ "Action (Outside +/- 2 Spread)",
      median_length > thresholds$upper_warn   | median_length < thresholds$lower_warn   ~ "Warning (Outside +/- 1 Spread)",
      TRUE ~ "Target / Stable"
    )
  )

# 6. Automated Triage Function for New Data
triage_lobster_year <- function(new_lengths, year_label = "Current Year") {
  # Force input to an atomic numeric vector and remove NAs
  clean_lengths <- as.numeric(unlist(new_lengths))
  clean_lengths <- clean_lengths[!is.na(clean_lengths)]
  
  n_sample <- length(clean_lengths)
  if (n_sample == 0) stop("No valid length measurements found for this year.")
  
  med <- median(clean_lengths)
  boot_se <- sd(replicate(1000, median(sample(clean_lengths, replace = TRUE))))
  
  status_msg <- if (med < thresholds$lower_action | med > thresholds$upper_action) {
    "RED: Action Required (Size truncated/abnormal - prioritize full assessment)"
  } else if (med < thresholds$lower_warn | med > thresholds$upper_warn) {
    "YELLOW: Warning (Moderate departure - check CPUE/recruitment)"
  } else {
    "GREEN: Stable (Within historical high-biomass baseline)"
  }
  
  cat(paste0(
    "--- Rapid Assessment Triage: ", year_label, " ---\n",
    "Sample Size (n):    ", n_sample, "\n",
    "Median Length:      ", round(med, 2), " mm (SE: +/- ", round(boot_se, 2), " mm)\n",
    "95% CI (Median):    [", round(med - 1.96 * boot_se, 2), ", ", round(med + 1.96 * boot_se, 2), "] mm\n",
    "Triage Status:      ", status_msg, "\n\n"
  ))
}

# Run the 2024 triage check
cat("Running 2024 Triage Check...\n")
lobster_data |>
  filter(year == 2024) |>
  pull(length_mm) |>
  triage_lobster_year(year_label = "2024")

# 7. Generate the Monitoring Control Chart
plot_indicator <- ggplot(annual_monitored, aes(x = year, y = median_length)) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = thresholds$lower_action, ymax = thresholds$upper_action, 
           fill = "goldenrod1", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = thresholds$lower_warn, ymax = thresholds$upper_warn, 
           fill = "forestgreen", alpha = 0.15) +
  geom_hline(yintercept = ref_median, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_line(color = "steelblue4", linewidth = 0.7) +
  geom_errorbar(aes(ymin = lower_se, ymax = upper_se), width = 0.25, color = "grey20", linewidth = 0.5) +
  geom_point(aes(color = status, size = n)) +
  scale_color_manual(values = c(
    "Target / Stable" = "forestgreen",
    "Warning (Outside +/- 1 Spread)" = "orange3",
    "Action (Outside +/- 2 Spread)" = "firebrick3"
  )) +
  scale_size_continuous(range = c(2, 5), name = "Port Sample (n)") +
  labs(
    title = "Lobster Carapace Length Monitoring Indicator",
    subtitle = paste0("Baseline: 2007-Present | Target Median = ", round(ref_median, 1), " mm"),
    x = "Year",
    y = "Median Carapace Length (mm)",
    color = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(plot_indicator)



# 7. Generate the Monitoring Control Chart
plot_indicator_0724 <- annual_monitored |>
  filter(year >= 2007) |>
  ggplot(aes(x = year, y = median_length)) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = thresholds$lower_action, ymax = thresholds$upper_action, 
           fill = "goldenrod1", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = thresholds$lower_warn, ymax = thresholds$upper_warn, 
           fill = "forestgreen", alpha = 0.15) +
  geom_hline(yintercept = ref_median, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_line(color = "steelblue4", linewidth = 0.7) +
  geom_errorbar(aes(ymin = lower_se, ymax = upper_se), width = 0.25, color = "grey20", linewidth = 0.5) +
  geom_point(aes(color = status, size = n)) +
  scale_color_manual(values = c(
    "Target / Stable" = "forestgreen",
    "Warning (Outside +/- 1 Spread)" = "orange3",
    "Action (Outside +/- 2 Spread)" = "firebrick3"
  )) +
  scale_size_continuous(range = c(2, 5), name = "Port Sample (n)") +
  labs(
    title = "Lobster Carapace Length Monitoring Indicator",
    subtitle = paste0("Baseline: 2007-Present | Target Median = ", round(ref_median, 1), " mm"),
    x = "Year",
    y = "Median Carapace Length (mm)",
    color = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(plot_indicator_0724)
