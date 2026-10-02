# Colleen Ahern
# CAZyme family-level (class) summary bar chart
# Shows which CAZy classes (GH, CE, GT, AA, CBM, PL) dominate each category
# (GF CA-specific / Cellulose-specific / Shared), built from all_caz

library(tidyverse)

# =============================================================================
# 1. Extract CAZy class from each gene's CAZ_anno
# =============================================================================
# A gene can have multiple CAZyme family hits (e.g., "GH31_1; CE1") — unnest
# so each family hit is counted separately, then extract just the class
# prefix (letters only, e.g. "GH31_1" -> "GH", "AA3_2" -> "AA", "CBM48" -> "CBM")

caz_long <- all_caz %>%
  dplyr::select(GeneID, category, CAZ_anno) %>%
  separate_rows(CAZ_anno, sep = "; ") %>%          # one row per CAZyme family hit
  mutate(caz_class = str_extract(CAZ_anno, "^[A-Za-z]+")) %>%
  filter(!is.na(caz_class))

# ── Count genes (family hits) per class, per category ─────────────────────────
class_counts <- caz_long %>%
  dplyr::count(category, caz_class, name = "n_genes") %>%
  complete(category, caz_class, fill = list(n_genes = 0))

# =============================================================================
# 2. Plot — grouped bar chart
# =============================================================================

category_colors <- c(
  "GF CA-specific"             = "saddlebrown",
  "Shared (GF CA + Cellulose)" = "turquoise1",
  "Cellulose-specific"        = "springgreen"
)

caz_class_plot <- ggplot(class_counts, aes(x = caz_class, y = n_genes, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = category_colors, name = "Category") +
  labs(
    x = "CAZyme class",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 12) +
  theme(
    text                  = element_text(size = 12),
    plot.title            = element_text(size = 12, hjust = 0.5, face = "bold"),
    axis.text             = element_text(size = 12),
    axis.title            = element_text(size = 12),
    legend.position       = c(0.005, 0.995),  # top-left corner
    legend.justification  = c(0, 1),        # anchor the legend's own top-left corner to that point
    legend.background = element_rect(fill = "white", colour = NA),  # keep white fill, remove auto border
    legend.margin         = margin(5, 8, 5, 8),
    legend.title          = element_text(size = 12, face = "bold"),
    legend.text           = element_text(size = 12),
    legend.key.size       = unit(0.4, "cm"),
    panel.border          = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line             = element_blank(),
    aspect.ratio          = 1,
    plot.margin           = margin(5, 5, 5, 5)
  ) 
caz_class_plot

# =============================================================================
# 4. Save
# =============================================================================

ggsave(plot = caz_class_plot,
       filename = file.path(FIG_DIR, paste0("Pchry_CAZy_class_grouped_bar_", DATE, ".pdf")),
       width = 5, height = 5, units = "in", dpi = 300)
ggsave(plot = caz_class_plot,
       filename = file.path(FIG_DIR, paste0("Pchry_CAZy_class_grouped_bar_", DATE, ".tiff")),
       width = 8, height = 6, units = "in", dpi = 200, bg = "white")
