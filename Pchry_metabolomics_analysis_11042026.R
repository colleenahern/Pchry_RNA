# Colleen Ahern
# Making Pchry metabolomics heatmaps
# 11/04/2026

# make final heatmap: metabolites upregulated for both Pchry + substrate vs Pchry + no substrate control and Pchry + substrate vs abiotic substrate control
# same as 08032026 code but performed log2 transformation and no scaling

# Colleen Ahern
# Pchry metabolomics: bio-vs-bio, bio-vs-abio, and their intersection
# Fully self-contained — run top to bottom

library(tidyverse)
library(pheatmap)
library(grid)
library(dplyr)

FIG_DIR  <- "/Users/colleenahern/Documents/Pchry/Figures"
ANNO_PATH <- "/Users/colleenahern/Documents/Pchry/Metabolomics_analysis/annotations_filtered_grouped_08032026_fin.csv"

my_palette <- colorRampPalette(c("blue", "black", "red"))(299)

# =============================================================================
# 1. Reusable pipeline: read files, annotate, build masked LFC matrix
# =============================================================================

read_and_prefix <- function(filepath) {
  fname  <- basename(filepath)
  prefix <- str_extract(fname, "(?<=^Pchry_)[^_]+")
  df <- read_csv(filepath, show_col_types = FALSE)
  other_cols <- colnames(df)[-1]
  df %>% rename_with(~ paste0(prefix, "_", .x), .cols = all_of(other_cols))
}

build_metabolomics_data <- function(dir, conditions, display_names, anno_path) {
  csv_files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  all_dfs <- purrr::map(csv_files, read_and_prefix)
  
  combined <- purrr::reduce(all_dfs, full_join, by = "...1") %>% 
    dplyr::rename(FeatureID = ...1)
  
  combined$Name_clean <- combined$FeatureID %>% 
    gsub("^F_\\d+_", "", .) %>% 
    gsub("_M_H_H2O$|_M_NH4$|_M_K$|_M_H$", "", .)
  
  combined$rowID <- str_extract(combined$FeatureID, "(?<=^F_)\\d+") %>% 
    as.numeric()
  
  anno <- read_csv(anno_path, show_col_types = FALSE)
  anno$rowID <- as.numeric(anno$rowID)
  
  combined <- combined %>% 
    dplyr::left_join(anno %>% dplyr::select(rowID, Annotation), by = "rowID")
  
  combined <- combined %>% dplyr::filter(!is.na(Annotation))
  
  combined$Annotation_fin <- ifelse(
    is.na(combined$Annotation), 
    paste0("F_", combined$rowID, ": ", combined$Name_clean), 
    paste0("F_", combined$rowID, ": ", combined$Annotation)
  )
  combined$Annotation_fin <- gsub("\\s*\\|.*$", "", combined$Annotation_fin)
  combined$Annotation_fin <- gsub("HYDROXYCINNAMIC ACID", "Hydroxycinnamic acid", combined$Annotation_fin)
  combined$Annotation_fin <- gsub("HYDROXYPHENYLLACTIC ACID", "Hydroxyphenyllactic acid", combined$Annotation_fin)
  combined$Annotation_fin <- gsub("Massbank:PB000509 ", "", combined$Annotation_fin)
  
  combined <- combined %>% 
    dplyr::relocate(rowID, Annotation, Annotation_fin, .after = FeatureID)
  
  lfc_matrix <- combined %>% 
    column_to_rownames("Annotation_fin") %>% 
    dplyr::select(all_of(paste0(conditions, "_logFC")))
  
  padj_matrix <- combined %>% 
    column_to_rownames("Annotation_fin") %>% 
    dplyr::select(all_of(paste0(conditions, "_adj.P.Val")))
  
  masked_lfc <- lfc_matrix
  for (i in seq_along(conditions)) {
    not_sig <- is.na(padj_matrix[[i]]) | padj_matrix[[i]] >= 0.05
    masked_lfc[not_sig, i] <- 0
  }
  colnames(masked_lfc) <- display_names
  
  sig_positive <- (masked_lfc > 0) & !is.na(masked_lfc)
  
  list(combined = combined, masked_lfc = masked_lfc, sig_positive = sig_positive)
}

# =============================================================================
# 2. Bio-vs-bio comparison: Pchry + substrate vs. Pchry + no substrate
# =============================================================================

biobio_conditions <- c("Cellulose", "GFCA", "SACA", "RCG", "PDB")
biobio_display <- c("Cellulose vs.\nNo substrate", "GF CA vs.\nNo substrate",
                    "SA CA vs.\nNo substrate", "RCG vs.\nNo substrate", "PDB vs.\nNo substrate")

biobio_data <- build_metabolomics_data(
  dir            = "//Users/colleenahern/Documents/Pchry/Metabolomics_analysis/Metaboanalyst_results/11042026_analysis_biovbio_filt",
  conditions     = biobio_conditions,
  display_names  = biobio_display,
  anno_path      = ANNO_PATH
)

biobio_sig_filtered <- biobio_data$masked_lfc[rowSums(biobio_data$sig_positive, na.rm = TRUE) > 0, ]
message("Bio-vs-bio: ", nrow(biobio_sig_filtered), " significantly upregulated metabolites")

biobio_hm <- pheatmap(
  as.matrix(biobio_sig_filtered),
  col           = my_palette,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  cellheight    = 20,
  cellwidth     = 80,
  fontsize_row  = 10,
  border_color  = "black",
  breaks        = seq(-10, 10, length.out = 300),
  fontsize      = 12,
  fontsize_col  = 12,
  angle_col     = 0
)

biobio_hm

pdf(file.path(FIG_DIR, "Metabolomics_heatmap_filt_bio_v_bio_11042026.pdf"),
    width = 14, height = 4)

grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.45, just = "center"))
grid.draw(biobio_hm$gtable)
popViewport()

dev.off()

# =============================================================================
# 3. Bio-vs-abio comparison: Pchry + substrate vs. abiotic substrate alone
# =============================================================================

bioabio_conditions <- c("Cellulose", "GFCA", "SACA", "RCG", "PDB", "Media")
bioabio_display <- c("Cellulose", "GF CA", "SA CA", "RCG", "PDB", "Media")

bioabio_data <- build_metabolomics_data(
  dir            = "/Users/colleenahern/Documents/Pchry/Metabolomics_analysis/Metaboanalyst_results/11042026_analysis_biovabio_filt",
  conditions     = bioabio_conditions,
  display_names  = bioabio_display,
  anno_path      = ANNO_PATH
)

bioabio_sig_filtered <- bioabio_data$masked_lfc[rowSums(bioabio_data$sig_positive, na.rm = TRUE) > 0, ]
message("Bio-vs-abio: ", nrow(bioabio_sig_filtered), " significantly upregulated metabolites")

bioabio_hm <- pheatmap(
  as.matrix(bioabio_sig_filtered),
  col           = my_palette,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  cellheight    = 12,
  cellwidth     = 50,
  fontsize_row  = 10,
  border_color  = "black",
  breaks        = seq(-10, 10, length.out = 300),
  fontsize      = 12,
  fontsize_col  = 12,
  angle_col     = "0"
)

bioabio_hm

pdf(file.path(FIG_DIR, "Metabolomics_heatmap_filt_bio_v_abio_11042026.pdf"),
    width = 16, height = 8)

grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.475, just = "center"))
grid.draw(bioabio_hm$gtable)
popViewport()

dev.off()

# =============================================================================
# 4. Intersection: metabolites significantly upregulated in BOTH comparisons,
#    for the same substrate condition — i.e., genuinely produced by
#    P. chrysosporium metabolizing that substrate, not from general fungal
#    metabolism (excluded via bio-vs-bio) or substrate leaching (excluded
#    via bio-vs-abio).
# =============================================================================

shared_conditions <- c("Cellulose", "GF CA", "PDB", "RCG", "SA CA")  # Media has no bio-vs-bio counterpart

# Map each shared condition to its ACTUAL column name in each dataset
# (biobio columns are named "X vs.\nNo substrate"; bioabio columns are plain "X")
biobio_col_lookup   <- c("Cellulose" = "Cellulose vs.\nNo substrate",
                         "GF CA"     = "GF CA vs.\nNo substrate",
                         "SA CA"     = "SA CA vs.\nNo substrate",
                         "RCG"       = "RCG vs.\nNo substrate",
                         "PDB"       = "PDB vs.\nNo substrate")

biobio_rowid  <- biobio_data$combined$rowID[match(rownames(biobio_data$masked_lfc), biobio_data$combined$Annotation_fin)]
bioabio_rowid <- bioabio_data$combined$rowID[match(rownames(bioabio_data$masked_lfc), bioabio_data$combined$Annotation_fin)]

intersection_flags <- map_dfc(shared_conditions, function(cond) {
  biobio_col <- biobio_col_lookup[[cond]]
  
  sig_biobio  <- tibble(rowID = biobio_rowid,  sig = biobio_data$sig_positive[, biobio_col])
  sig_bioabio <- tibble(rowID = bioabio_rowid, sig = bioabio_data$sig_positive[, cond])
  
  matched <- sig_biobio %>%
    left_join(sig_bioabio, by = "rowID", suffix = c("_biobio", "_bioabio")) %>%
    mutate(both_sig = sig_biobio & sig_bioabio) %>%
    pull(both_sig)
  
  set_names(list(matched), cond)
})

rows_to_keep <- rowSums(intersection_flags, na.rm = TRUE) > 0
message(sum(rows_to_keep), " metabolites significantly upregulated in BOTH comparisons for at least one condition")

final_matrix <- biobio_data$masked_lfc[rows_to_keep, biobio_col_lookup[shared_conditions]]
colnames(final_matrix) <- shared_conditions  # rename to plain condition names for the final figure

# ── Reorder columns for display: Cellulose, GF CA, SA CA, RCG, PDB ────────────
final_matrix <- final_matrix[, c("Cellulose", "GF CA", "SA CA", "RCG", "PDB")]
rownames(final_matrix) <- gsub(
  "4-acetyloxy-8-\\(3-oxo-2-pent-2-enylcyclopenten-1-yl\\)octanoic acid",
  "4-acetyloxy-8-(3-oxo-2-pent-2-\nenylcyclopenten-1-yl)octanoic acid",
  rownames(final_matrix)
)

final_hm <- pheatmap(
  as.matrix(final_matrix),
  col           = my_palette,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  cellheight    = 40,
  cellwidth     = 60,
  fontsize_row  = 10,
  border_color  = "black",
  breaks        = seq(-4, 4, length.out = 300),
  fontsize      = 12,
  fontsize_col  = 12,
  angle_col     = "0"
)

final_hm

pdf(file.path(FIG_DIR, "Metabolomics_heatmap_intersection_11042026.pdf"),
    width = 8, height = 2)

grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, just = "center"))
grid.draw(final_hm$gtable)
popViewport()

dev.off()

# Make supplementary tables showing mz and retention time values for all significantly upregulated metabolites for the bio v bio and bio v abio comparisons
library(tidyverse)

FEAT_PATH <- "/Users/colleenahern/Downloads/Supernatants-ec58bb01006c49b28da858e504f6f135-aligned_features.csv"
OUT_DIR   <- "/Users/colleenahern/Documents/Pchry/Metabolomics_analysis/Metaboanalyst_results"

BIOBIO_DIR  <- "/Users/colleenahern/Documents/Pchry/Metabolomics_analysis/Metaboanalyst_results/11042026_analysis_biovbio_filt"
BIOABIO_DIR <- "/Users/colleenahern/Documents/Pchry/Metabolomics_analysis/Metaboanalyst_results/11042026_analysis_biovabio_filt"

# ---- Feature metadata: rowID -> m/z, RT ----
feat_meta <- read_csv(FEAT_PATH, show_col_types = FALSE) %>%
  dplyr::select(rowID  = `row ID`,
                mz     = `row m/z`,
                RT_min = `row retention time`) %>%
  mutate(rowID = as.numeric(rowID))

# ---- Read one comparison file ----
read_comparison <- function(filepath) {
  comparison <- basename(filepath) %>%
    str_remove("_covariate_result_\\d+\\.csv$") %>%
    str_remove("\\.csv$")
  
  read_csv(filepath, show_col_types = FALSE) %>%
    dplyr::rename(Feature = 1) %>%
    mutate(Comparison = comparison, .before = 1)
}

# ---- Stack, attach m/z and RT, keep only significant upregulated ----
build_sig_table <- function(dir) {
  list.files(dir, pattern = "\\.csv$", full.names = TRUE) %>%
    purrr::map_dfr(read_comparison) %>%
    dplyr::filter(adj.P.Val < 0.05, logFC > 0) %>%
    mutate(rowID = as.numeric(str_extract(Feature, "(?<=^F_)\\d+"))) %>%
    dplyr::left_join(feat_meta, by = "rowID") %>%
    dplyr::arrange(Comparison, desc(logFC)) %>%
    dplyr::select(Comparison, Feature, mz, RT_min, logFC, P.Value, adj.P.Val)
}

biobio_sig  <- build_sig_table(BIOBIO_DIR)
bioabio_sig <- build_sig_table(BIOABIO_DIR)

# ---- Checks ----
nrow(biobio_sig)
nrow(bioabio_sig)
sum(is.na(biobio_sig$mz))    # should be 0
sum(is.na(bioabio_sig$mz))   # should be 0
table(biobio_sig$Comparison)
table(bioabio_sig$Comparison)
head(biobio_sig, 10)

# ---- Write out ----
write_csv(biobio_sig,  file.path(OUT_DIR, "biovbio_significant_11042026.csv"))
write_csv(bioabio_sig, file.path(OUT_DIR, "biovabio_significant_11042026.csv"))
message("Done.")