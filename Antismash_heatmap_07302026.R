# Heatmap for antiSMASH biosynthetic gene clusters (cleaned)

library(tidyverse)
library(pheatmap)
library(grid)

BASE    <- "/Users/colleenahern/Documents/Pchry"
FIG_DIR <- file.path(BASE, "Figures")
DATE    <- "07302026"

# =============================================================================
# 1. Import DESeq2 results
# =============================================================================

vary_43caz <- read_csv(file.path(BASE, "DESeq2_analysis/Tables/vary43_jgi_allconds_anno_cazanno_08272025.csv"),
                       show_col_types = FALSE)

# =============================================================================
# 2. Mask non-significant LFCs to 0 (finds padj column by name, not position)
# =============================================================================

mask_nonsig_lfc <- function(df, alpha = 0.05) {
  lfc_cols <- grep("_log2FoldChange$", colnames(df), value = TRUE)
  for (lfc_col in lfc_cols) {
    padj_col <- gsub("_log2FoldChange$", "_padj", lfc_col)
    stopifnot(padj_col %in% colnames(df))
    not_sig <- df[[padj_col]] >= alpha | is.na(df[[padj_col]])
    df[[lfc_col]][not_sig] <- 0
  }
  df
}

vary_43caz_masked <- mask_nonsig_lfc(vary_43caz)

# =============================================================================
# 3. Merge with antiSMASH biosynthetic gene cluster annotations
# =============================================================================

cbgs <- read.delim(file.path(BASE, "Antismash_results/Core_antiSMASH_biosynthetic_genes_12052025.txt"))

cbgs$antiSMASH_BGC_Type <- gsub("terpene", "Terpene", cbgs$antiSMASH_BGC_Type)
cbgs$antiSMASH_BGC_Type <- gsub("Terpene ", "Terpene", cbgs$antiSMASH_BGC_Type)
cbgs$antiSMASH_BGC_Type <- gsub("betalactone", "Betalactone", cbgs$antiSMASH_BGC_Type)
unique(cbgs$antiSMASH_BGC_Type)

# NOTE: merge by ProteinID — confirmed as the shared key between cbgs and vary_43caz_masked
vary_43 <- merge(cbgs, vary_43caz_masked, by = "ProteinID")

# =============================================================================
# 4. Row annotation (BGC type)
# =============================================================================

# ── Clean up row labels: remove "proteinID = " prefix if present ──────────────
vary_43$Name <- gsub("proteinID = ", "", vary_43$Name)

# ── Capitalize just the first letter of each Name (leave the rest as-is) ──────
vary_43$Name <- paste0(toupper(substr(vary_43$Name, 1, 1)), substr(vary_43$Name, 2, nchar(vary_43$Name)))

annotation_rowpc <- data.frame(BGC_Type = vary_43$antiSMASH_BGC_Type)
rownames(annotation_rowpc) <- vary_43$Name

ann_colorspc <- list(
  BGC_Type = c(
    `NRPS-like`           = "#FE00FA",
    `Terpene`             = "#16FF32",
    NRPS                  = "#FBE426",
    T3PKS                 = "#FEAF16",
    Betalactone           = "#1CFFCE",
    `Terpene-precursor`   = "#B00068",
    `NRPS-like, Terpene`  = "#AA0DFE",
    `NRPS-like, T1PKS`    = "#DEA0FD"
  )
)

# =============================================================================
# 5. Column labels and data matrix
# =============================================================================

lfc_cols <- c("SACA_v_blank_log2FoldChange", "GFCA_v_blank_log2FoldChange",
              "cellulose_v_blank_log2FoldChange", "GFCA_v_cellulose_log2FoldChange")

newnames <- c("SA CA vs.\n No substrate", "GF CA vs.\n No substrate",
              "Cellulose vs.\n No substrate", "GF CA vs.\n Cellulose")

df_figpc <- data.matrix(vary_43[, lfc_cols])
rownames(df_figpc) <- vary_43$Name

# =============================================================================
# 6. Heatmap
# =============================================================================

draw_colnames_horizontal <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"),
           vjust = 1.1, hjust = 0.5, rot = 0, gp = gpar(...))
}
assignInNamespace("draw_colnames", draw_colnames_horizontal, ns = asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(299)

pc_regulation_heatmap <- pheatmap(
  df_figpc,
  col                  = my_palette,
  cellwidth            = 80,
  cellheight           = 15,
  cluster_rows         = TRUE,
  cluster_cols         = FALSE,
  fontsize             = 10,
  show_rownames        = TRUE,
  labels_col           = newnames,
  annotation_row       = annotation_rowpc,
  annotation_colors    = ann_colorspc,
  annotation_names_row = FALSE,
  border_color         = "black",
  breaks               = seq(-4, 4, length.out = 300)
)

pc_regulation_heatmap$gtable$heights[length(pc_regulation_heatmap$gtable$heights)] <- unit(1.7, "cm")
grid.newpage()
grid.draw(pc_regulation_heatmap$gtable)

# =============================================================================
# 7. Save
# =============================================================================

ggsave(plot = pc_regulation_heatmap,
       filename = file.path(FIG_DIR, paste0("Pchry_antismash_heatmap_", DATE, ".pdf")),
       width = 12.5, height = 9, units = "in", dpi = 300)
ggsave(plot = pc_regulation_heatmap,
       filename = file.path(FIG_DIR, paste0("Pchry_antismash_heatmap_", DATE, ".tiff")),
       width = 12, height = 9, units = "in", dpi = 200, bg = "white")

message("Done.")
