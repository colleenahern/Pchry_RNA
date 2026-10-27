# Colleen Ahern
# Making Pchry CO2 figure
# 07/21/2026

library(tidyverse)
library(ggpubr)
library(cowplot)

# ── Path ──────────────────────────────────────────────────────────────────────
FILE <- "/Users/colleenahern/Documents/Pchry/Biogas/Pchry_biogas_org_07222026.csv"

# ── Read raw wide-format data ───────────────────────────────────────────────────
raw <- read_csv(FILE, col_names = TRUE, show_col_types = FALSE)

# =============================================================================
# The file has 5 condition "blocks" side by side, each with its own replicate
# columns + summary columns. Column names are now unique per block, so we can
# reference them directly by name rather than guessing positions.
# Note: Cellulose and GF_CA share "Pchry_cellulose_time"; HDPE and the blank
# Pchry control share "Pchry_HDPE_time" (confirmed from the header row).
# rep_prefix/n_reps are used to pull individual replicate values for the
# endpoint significance bar plot (t-tests need per-replicate data, not just
# the pre-computed avg/std columns).
# =============================================================================

block_map <- list(
  Cellulose              = list(time = "Pchry_cellulose_time", avg = "Pchry_cellulose_AA_avg", std = "Pchry_cellulose_AA_std",
                                gc_avg = "Pchry_cellulose_gC_avg", gc_std = "Pchry_cellulose_gC_std",
                                rep_prefix = "Pchry_cellulose_", reps = c(2, 6)),
  "Lowly acetylated CA (GF CA)"  = list(time = "Pchry_cellulose_time", avg = "Pchry_GFCA_AA_avg", std = "Pchry_GFCA_AA_std",
                                gc_avg = "Pchry_GFCA_gC_avg", gc_std = "Pchry_GFCA_gC_std",
                                rep_prefix = "Pchry_GF_CA_", reps = c(6, 7, 8)),
  "Highly acetylated CA (SA CA)" = list(time = "Pchry_SACA_time", avg = "Pchry_SACA_AA_avg", std = "Pchry_SACA_AA_std",
                                gc_avg = "Pchry_SACA_gC_avg", gc_std = "Pchry_SACA_gC_std",
                                rep_prefix = "Pchry_SA_CA_", reps = c(7, 8)),
  HDPE                   = list(time = "Pchry_HDPE_time", avg = "Pchry_HDPE_AA_avg", std = "Pchry_HDPE_AA_std",
                                gc_avg = "Pchry_HDPE_gC_avg", gc_std = "Pchry_HDPE_gC_std",
                                rep_prefix = "Pchry_HDPE_", reps = c(2, 3)),
  "No substrate"         = list(time = "Pchry_HDPE_time", avg = "Pchry_AA_avg", std = "Pchry_AA_std",
                                gc_avg = "Pchry_gC_avg", gc_std = "Pchry_gC_std",
                                rep_prefix = "Pchry_", reps = c(6, 7, 8))
)

# ── Helper: pull out one condition's block as a tidy data frame ───────────────
extract_block <- function(raw, cols, condition_name, metric = c("AA", "gC")) {
  metric <- match.arg(metric)
  avg_col <- if (metric == "AA") cols$avg else cols$gc_avg
  std_col <- if (metric == "AA") cols$std else cols$gc_std
  
  tibble(
    condition = condition_name,
    time      = raw[[cols$time]],
    avg       = raw[[avg_col]],
    std       = raw[[std_col]]
  ) %>%
    filter(!is.na(time), !is.na(avg))  # drop empty trailing rows
}

# ── Helper: pull individual replicate values at the FINAL timepoint for one condition ──
extract_endpoint <- function(raw, cols, condition_name, metric = c("AA", "gC")) {
  metric  <- match.arg(metric)
  avg_col <- if (metric == "AA") cols$avg else cols$gc_avg
  
  # Last row where this condition has a non-NA average = its final measured timepoint
  last_row <- max(which(!is.na(raw[[avg_col]])))
  
  suffix   <- if (metric == "AA") "" else "_gC"
  rep_cols <- paste0(cols$rep_prefix, cols$reps, suffix)
  rep_cols <- rep_cols[rep_cols %in% colnames(raw)]  # guard against any naming mismatch
  
  tibble(
    condition = condition_name,
    replicate = rep_cols,
    value     = as.numeric(raw[last_row, rep_cols])
  ) %>%
    filter(!is.na(value))
}

# =============================================================================
# Build long-format data frame
# ── Change metric = "AA" to metric = "gC" below to plot grams-carbon instead ──
# =============================================================================

# ── Shared color mapping (named, so colors stay tied to condition identity) ──
condition_colors <- c(
  "Cellulose"              = "#1B9E77",
  "Lowly acetylated CA (GF CA)"    = "#D95F02",
  "Highly acetylated CA (SA CA)"   = "#E7298A",
  "HDPE"                   = "#7570B3",
  "No substrate"           = "grey40"
)
condition_levels <- names(condition_colors)

METRIC <- "gC"   # "AA" = raw accumulated biogas; "gC" = grams-carbon equivalent

growth_df <- imap_dfr(block_map, ~ extract_block(raw, .x, .y, metric = METRIC))

# ── Set condition order for consistent legend/color ordering ──────────────────
growth_df$condition <- factor(growth_df$condition, levels = condition_levels)

# =============================================================================
# Growth curve plot
# =============================================================================

y_label <- if (METRIC == "AA") expression("Cumulative CO"[2]*" (mL)") else expression("Cumulative CO"[2]*" (gC)")

growth_plot <- ggplot(growth_df, aes(x = time, y = avg, color = condition)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = avg - std, ymax = avg + std), width = 50, alpha = 0.6) +
  scale_color_manual(values = condition_colors) +
  labs(
    x     = "Time (hours)",
    y     = y_label,
    color = "Substrate"
  ) +
  theme_classic(base_size = 12) +
  theme(
    text         = element_text(size = 12),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 12)
  )
growth_plot

growth_plot
ggsave("/Users/colleenahern/Documents/Pchry/Figures/Pchry_growth_curve_10262026.pdf",
       plot = growth_plot, width = 8, height = 6)


# =============================================================================
# Student's t-test comparison 
# =============================================================================

comparisons_list <- list(
  c("Cellulose", "GF CA"),
  c("Cellulose", "SA CA"),
  c("Cellulose", "HDPE"),
  c("Cellulose", "No substrate"),
  c("GF CA", "SA CA"),
  c("GF CA", "HDPE"),
  c("GF CA", "No substrate")
)

endpoint_df <- imap_dfr(block_map, ~ extract_endpoint(raw, .x, .y, metric = METRIC))

endpoint_df$condition <- recode(endpoint_df$condition,
                                "Lowly acetylated CA (GF CA)" = "GF CA",
                                "Highly acetylated CA (SA CA)" = "SA CA"
)
condition_levels <- recode(condition_levels,
                           "Lowly acetylated CA (GF CA)" = "GF CA",
                           "Highly acetylated CA (SA CA)" = "SA CA"
)
endpoint_df$condition <- factor(endpoint_df$condition, levels = condition_levels)

# rename condition_colors keys to match the short names
names(condition_colors)[names(condition_colors) == "Lowly acetylated CA (GF CA)"] <- "GF CA"
names(condition_colors)[names(condition_colors) == "Highly acetylated CA (SA CA)"] <- "SA CA"

bar_y_label <- if (METRIC == "AA") expression("Total CO"[2]*" (mL)") else expression("Total CO"[2]*" (gC)")

endpoint_bar <- ggbarplot(endpoint_df, x = "condition", y = "value",
                          add = "mean_sd", fill = "condition",
                          xlab = "Substrate") +
  ylab(bar_y_label) +
  scale_fill_manual(values = condition_colors) +
  scale_x_discrete(labels = function(x) ifelse(x == "No substrate", gsub(" ", "\n", x), x)) +
  stat_compare_means(method = "t.test",
                     method.args = list(var.equal = TRUE),
                     label = "p.signif",
                     comparisons = comparisons_list) +
  theme_classic(base_size = 12) +
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )
endpoint_bar

combined_plot <- plot_grid(growth_plot, endpoint_bar,
                           ncol = 1, nrow = 2,
                           rel_heights = c(.75, 1),
                           labels = c("A)", "B)"))
combined_plot

# ── Save ────────────────────────────────────────────────────────────────────
ggsave("/Users/colleenahern/Documents/Pchry/Figures/Pchry_growth_curve_08152026.pdf",
       plot = combined_plot, width = 8, height = 10)
ggsave("/Users/colleenahern/Documents/Pchry/Figures/Pchry_growth_curve_08152026.tiff",
       plot = combined_plot, width = 8, height = 12, dpi = 300)


# =============================================================================
# Welch's t-test results (for comparison against the Student's test used in the plot)
# =============================================================================
run_ttest <- function(pair, var_equal) {
  d1 <- endpoint_df$value[endpoint_df$condition == pair[1]]
  d2 <- endpoint_df$value[endpoint_df$condition == pair[2]]
  tt <- t.test(d1, d2, var.equal = var_equal)
  tibble(group1 = pair[1], group2 = pair[2], n1 = length(d1), n2 = length(d2), p = tt$p.value)
}

welch_pvals <- map_dfr(comparisons_list, run_ttest, var_equal = FALSE) %>%
  mutate(p.signif = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  ))

print(welch_pvals, n = Inf)
