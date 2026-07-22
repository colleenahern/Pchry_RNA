# Colleen Ahern
# 07/21/2026
# Complete DESeq2 analysis on Pchry count matrix from featureCounts
# Exclude SACA5 sample (based on MDS)
# Tutorial reference: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# ── Libraries ─────────────────────────────────────────────────────────────────

library(DESeq2)
library(tidyverse)  

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE <- "/Users/colleenahern/Documents/Pchry/DESeq2_analysis"
DATE <- "07212026"

# =============================================================================
# HELPERS
# =============================================================================

# ── Run one pairwise contrast, tidy it, and count significant genes ──────────
run_contrast <- function(dds, group_var, level_a, level_b, alpha = 0.05) {
  label <- paste0(level_a, "_v_", level_b)
  
  res <- results(dds, contrast = c(group_var, level_a, level_b))
  res_df <- as.data.frame(as.matrix(res))
  colnames(res_df) <- paste(label, colnames(res_df), sep = "_")
  
  n_sig <- sum(res_df[[paste0(label, "_padj")]] < alpha, na.rm = TRUE)
  message(label, ": ", n_sig, " genes significant at padj < ", alpha)
  
  list(df = res_df, label = label, n_sig = n_sig)
}

# =============================================================================
# 1. Import count data and sample info
# =============================================================================

pcseqdata   <- read.delim(file.path(BASE, "featcount/pchry_countmat_jgi_07312025.txt"),
                          stringsAsFactors = FALSE)
pcsampleinfo <- read.delim(file.path(BASE, "pcsampleinfo.txt"),
                           stringsAsFactors = TRUE)

dim(pcseqdata)
pcsampleinfo

# ── Store GeneID as rownames, drop non-count columns (chr, start, end, strand, length) ──
pccountdata <- pcseqdata
rownames(pccountdata) <- pccountdata$Geneid
pccountdata <- pccountdata[, -c(1:6)]

# ── Match column names to sample info ─────────────────────────────────────────
pcsampleinfo$libraryName_jgi <- gsub("[/-]", ".", pcsampleinfo$libraryName_jgi)

stopifnot(all(colnames(pccountdata) %in% pcsampleinfo$libraryName_jgi))

pccountdata <- pccountdata[, as.character(pcsampleinfo$libraryName_jgi)]
stopifnot(all(colnames(pccountdata) == pcsampleinfo$libraryName_jgi))

colnames(pccountdata) <- pcsampleinfo$sampleName
rownames(pcsampleinfo) <- pcsampleinfo$sampleName

pccts     <- pccountdata
pccoldata <- pcsampleinfo
stopifnot(all(rownames(pccoldata) == colnames(pccts)))

# =============================================================================
# 2. Remove outlier sample(s)
# =============================================================================

samples_to_remove <- "SACA5"

pccts_subset     <- pccts[, !(colnames(pccts) %in% samples_to_remove)]
pccoldata_subset <- pccoldata[!(rownames(pccoldata) %in% samples_to_remove), ]

stopifnot(all(rownames(pccoldata_subset) == colnames(pccts_subset)))
stopifnot(!(samples_to_remove %in% colnames(pccts_subset)))  # confirm removal actually happened
message("Removed sample(s): ", paste(samples_to_remove, collapse = ", "),
        " | Samples remaining: ", ncol(pccts_subset))

# =============================================================================
# 3. Build DESeqDataSet
# =============================================================================

pcdds <- DESeqDataSetFromMatrix(countData = pccts_subset,
                                colData   = pccoldata_subset,
                                design    = ~ groupName)

# ── Pre-filter genes with zero counts across all samples ──────────────────────
n_before <- nrow(pcdds)
pcdds <- pcdds[rowSums(counts(pcdds)) > 0, ]
message("Pre-filter: ", n_before, " -> ", nrow(pcdds), " genes")

# ── Set factor levels (blank = reference) ──────────────────────────────────────
pcdds$groupName <- factor(pcdds$groupName, levels = c("blank", "cellulose", "GFCA", "SACA"))

# ── Run DESeq2 ──────────────────────────────────────────────────────────────────
pcdds <- DESeq(pcdds)

# =============================================================================
# 4. Pairwise contrasts (6 total)
# =============================================================================

comparisons <- list(
  c("SACA",      "blank"),
  c("GFCA",      "blank"),
  c("cellulose", "blank"),
  c("cellulose", "SACA"),
  c("GFCA",      "SACA"),
  c("GFCA",      "cellulose")
)

results_list <- lapply(comparisons, function(pair) {
  run_contrast(pcdds, "groupName", pair[1], pair[2])
})
names(results_list) <- sapply(results_list, `[[`, "label")

# ── Summary of significant gene counts per comparison ──────────────────────────
sig_summary <- tibble(
  comparison  = names(results_list),
  n_sig_genes = sapply(results_list, `[[`, "n_sig")
)
print(sig_summary)

# ── Confirm gene order matches across all result data frames ──────────────────
gene_order_ok <- all(sapply(results_list, function(x) all(rownames(x$df) == rownames(results_list[[1]]$df))))
stopifnot(gene_order_ok)

# =============================================================================
# 5. Combine into one wide table
# =============================================================================

Pc_DESeq <- bind_cols(lapply(results_list, `[[`, "df"))
Pc_DESeq <- Pc_DESeq %>%
  rownames_to_column("GeneID") %>%
  relocate(GeneID)  # GeneID first, regardless of how many comparison columns exist

write_csv(Pc_DESeq, file.path(BASE, paste0("Tables/Pc_DESeq_jgi_", DATE, ".csv")))

message("Done. Output written to: ", file.path(BASE, paste0("Tables/Pc_DESeq_jgi_", DATE, ".csv")))

