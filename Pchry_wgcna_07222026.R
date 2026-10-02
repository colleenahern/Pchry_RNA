# Colleen Ahern
# WGCNA analysis on Pchry RNA-seq data + GFCA/Cellulose CAZyme specificity heatmap
# Cleaned and combined: 07/22/2026

library(tidyverse)
library(WGCNA)
library(DESeq2)
library(gridExtra)
library(CorLevelPlot)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(cowplot)
library(ggpubr)
library(ape)

BASE <- "/Users/colleenahern/Documents/Pchry"
FIG_DIR <- file.path(BASE, "Figures")
TAB_DIR <- file.path(BASE, "DESeq2_analysis/Tables")
DATE    <- "07222026"

# =============================================================================
# 1. Read in and prepare data
# =============================================================================

counts <- read.delim(file.path(BASE, "DESeq2_analysis/featcount/pchry_countmat_jgi_07312025.txt"), stringsAsFactors = FALSE)
data <- counts[, c(1, 7:26)] %>% column_to_rownames(var = "Geneid")

pcsampleinfo <- read.delim(file.path(BASE, "DESeq2_analysis/pcsampleinfo.txt"), stringsAsFactors = TRUE)
pcsampleinfo$libraryName_jgi <- gsub("[/-]", ".", pcsampleinfo$libraryName_jgi)
row.names(pcsampleinfo) <- pcsampleinfo$libraryName_jgi

# =============================================================================
# 2. QC — outlier detection
# =============================================================================

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

data <- data[gsg$goodGenes == TRUE, ]
message("Genes after goodSamplesGenes filter: ", nrow(data))

htree <- hclust(dist(t(data)), method = "average")
plot(htree)  # SACA5 clusters as an outlier

pca <- prcomp(t(data))
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0("PC1: ", pca.var.percent[1], " %"),
       y = paste0("PC2: ", pca.var.percent[2], " %"))  # SACA5 confirmed as outlier

samples.to.be.excluded <- "HISAT2_jgi.RNA15.Pchry.SACA5_trimmed_aligned.bam"
data.subset <- data[, !(colnames(data) %in% samples.to.be.excluded)]

# =============================================================================
# 3. Normalization
# =============================================================================

colData <- pcsampleinfo %>% filter(!row.names(.) %in% samples.to.be.excluded)

stopifnot(all(rownames(colData) %in% colnames(data.subset)))
stopifnot(all(colnames(data.subset) %in% rownames(colData)))

data.subset <- data.subset[, as.character(colData$libraryName_jgi)]
stopifnot(all(rownames(colData) == colnames(data.subset)))

colnames(data.subset) <- colData$sampleName
rownames(colData)     <- colData$sampleName
stopifnot(all(rownames(colData) == colnames(data.subset)))

dds <- DESeqDataSetFromMatrix(countData = data.subset, colData = colData, design = ~ 1)

# Remove genes with < 15 counts in more than 75% of samples (WGCNA RNA-seq FAQ recommendation)
dds75 <- dds[rowSums(counts(dds) >= 15) >= round(ncol(dds) * 0.75), ]
message("Genes after WGCNA count filter: ", nrow(dds75))

dds_norm <- vst(dds75)
norm.counts <- t(assay(dds_norm))

# =============================================================================
# 4. Network construction
# =============================================================================

power <- c(1:10, seq(12, 50, by = 2))

sft <- pickSoftThreshold(norm.counts, powerVector = power, networkType = "signed", verbose = 5)
sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

norm.counts[] <- sapply(norm.counts, as.numeric)
soft_power <- 24

temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# =============================================================================
# 5. Module eigengenes
# =============================================================================

module_eigengenes <- bwnet$MEs
head(module_eigengenes)
table(bwnet$colors)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)
# grey = genes that don't fall into any other module

# =============================================================================
# 5b. Module membership (kME)
# =============================================================================

module_membership <- signedKME(norm.counts, module_eigengenes)

own_module_kME <- tibble(
  GeneID = rownames(module_membership),
  module = bwnet$colors,
  kME    = sapply(seq_len(nrow(module_membership)), function(i) {
    col <- paste0("kME", bwnet$colors[i])
    if (col %in% colnames(module_membership)) module_membership[i, col] else NA
  })
)

own_module_kME %>%
  group_by(module) %>%
  summarize(n = n(), median_kME = median(kME, na.rm = TRUE),
            pct_below_0.7 = round(mean(kME < 0.7, na.rm = TRUE) * 100, 1)) %>%
  print(n = Inf)

# =============================================================================
# 6. Relate modules to traits
# =============================================================================

traits <- colData %>%
  mutate(
    ca_biodeg        = ifelse(grepl("GFCA", groupName), 1, 0),
    cellulose_biodeg = ifelse(grepl("cellulose", groupName), 1, 0)
  ) %>%
  dplyr::select(ca_biodeg, cellulose_biodeg)

nSamples <- nrow(norm.counts)
nGenes   <- ncol(norm.counts)

module.trait.corr       <- cor(module_eigengenes, traits, use = "p")
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, traits, by = "row.names") %>%
  column_to_rownames(var = "Row.names") %>%
  dplyr::rename(Cellulose = cellulose_biodeg, `GF CA` = ca_biodeg)

trait_cols  <- c("Cellulose", "GF CA")
module_cols <- setdiff(colnames(heatmap.data), trait_cols)

wgcnaplot <- CorLevelPlot(heatmap.data,
                          x = trait_cols,
                          y = module_cols,
                          fontCorval = 1, fontLabX = 1, fontLabY = 1,
                          cexCorval = 1, cexLabX = 1, cexLabY = 1,
                          col = c("blue1", "skyblue", "white", "pink", "red"))
wgcnaplot

# pdf(file.path(FIG_DIR, paste0("WGCNA_GFCAcellulose_corrplot_", DATE, ".pdf")), width = 8, height = 8)
# print(wgcnaplot)
# dev.off()

# tiff(file.path(FIG_DIR, paste0("WGCNA_GFCAcellulose_corrplot_", DATE, ".tiff")),
#      width = 8, height = 8, units = "in", res = 300)
# print(wgcnaplot)
# dev.off()

# =============================================================================
# 7. Module size barplot
# =============================================================================

module.gene.mapping <- as.data.frame(bwnet$colors)

barp <- as.data.frame(table(bwnet$colors))
barp$Var1 <- paste0("ME", barp$Var1)

bpwgcna <- ggplot(barp, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  ylim(0, 5000) +
  scale_fill_manual(values = c(
    "MEblack" = "black", "MEblue" = "blue", "MEbrown" = "brown", "MEgreen" = "green",
    "MEgrey" = "grey", "MEmagenta" = "magenta", "MEpink" = "deeppink", "MEpurple" = "purple",
    "MEred" = "red", "MEturquoise" = "turquoise", "MEyellow" = "yellow"
  )) +
  labs(y = "Module size (N)", x = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24),
        legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(bpwgcna)
bpwgcna

# ggsave(plot = bpwgcna, filename = file.path(FIG_DIR, paste0("Pchry_wgcna_barplot_", DATE, ".pdf")),
#        width = 12, height = 10, units = "in", dpi = 300)
# ggsave(plot = bpwgcna, filename = file.path(FIG_DIR, paste0("Pchry_wgcna_barplot_", DATE, ".tiff")),
#        width = 12, height = 10, units = "in", dpi = 200, bg = "white")

# =============================================================================
# 8. Extract genes per module of interest, annotate, merge with DESeq2 results
# =============================================================================

gff <- read.gff(file.path(BASE, "DESeq2_analysis/Phchr4_2_GeneCatalog_20210419.gff3"), GFF3 = TRUE)
gff <- gff[gff$type == "gene", ]

gene_anno <- gff %>%
  dplyr::select(seqid, everything()) %>%
  mutate(attributes = str_split(attributes, ";")) %>%
  unnest_wider(where(is.list), names_sep = "") %>%
  dplyr::select(attributes1, attributes4, attributes5, attributes6) %>%
  dplyr::rename(GeneID = attributes1, Protein_name = attributes4,
                ProteinID = attributes5, TranscriptID = attributes6) %>%
  mutate(
    GeneID       = gsub("ID=", "", GeneID),
    Protein_name = gsub("product_name=", "", Protein_name),
    ProteinID    = gsub("proteinId=", "", ProteinID),
    TranscriptID = gsub("transcriptId=", "", TranscriptID)
  )

caz_anno <- read.csv(file.path(BASE, "DESeq2_analysis/Pchry_dbcan/cazlist_08072025.csv"))[, c("ProteinID", "CAZ_anno")]
geneanno <- merge(gene_anno, caz_anno, by = "ProteinID", all.x = TRUE)

vary_43    <- read_csv(file.path(TAB_DIR, "vary43_jgi_allconds_anno_cazanno_07222026.csv"), show_col_types = FALSE)
vary_43sub <- vary_43[, grepl("GeneID|log2FoldChange|_padj", colnames(vary_43))]

process_module <- function(color, save_date = DATE, kME_threshold = 0.7) {
  genes <- module.gene.mapping %>%
    filter(bwnet$colors == color) %>%
    rownames() %>%
    as.data.frame()
  colnames(genes) <- "GeneID"
  
  annotated <- merge(genes, geneanno, by = "GeneID")
  deseq_merged <- merge(annotated, vary_43sub, by = "GeneID") %>%
    left_join(own_module_kME %>% dplyr::select(GeneID, kME), by = "GeneID") %>%
    mutate(across(where(is.numeric), ~round(., 3))) %>%
    arrange(desc(kME))
  
  write_csv(deseq_merged, file.path(TAB_DIR, paste0("WGCNA_me", color, "_deseq_", save_date, ".csv")))
  
  high_conf <- deseq_merged %>% filter(kME >= kME_threshold)
  write_csv(high_conf, file.path(TAB_DIR, paste0("WGCNA_me", color, "_deseq_highconf_", save_date, ".csv")))
  
  n_weak <- sum(deseq_merged$kME < kME_threshold, na.rm = TRUE)
  message("ME", color, ": ", nrow(deseq_merged), " genes total, ",
          nrow(high_conf), " high-confidence (kME >= ", kME_threshold, "), ",
          n_weak, " below threshold")
  
  deseq_merged
}

mepurple_deseq <- process_module("purple")  # + correlated with GFCA
mebrown_deseq  <- process_module("brown")   # + correlated with GFCA
meblack_deseq  <- process_module("black")   # + correlated with GFCA
megreen_deseq  <- process_module("green")   # + correlated with cellulose
meblue_deseq   <- process_module("blue")    # + correlated with both GFCA and cellulose

# =============================================================================
# 9. GFCA-specific / Cellulose-specific / Shared CAZyme gene classification
# =============================================================================

# ── Mask LFC to 0 for non-significant genes (finds padj col by name, not position) ──
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

mepurple_deseq_b <- mask_nonsig_lfc(mepurple_deseq)
mebrown_deseq_b  <- mask_nonsig_lfc(mebrown_deseq)
meblack_deseq_b  <- mask_nonsig_lfc(meblack_deseq)
megreen_deseq_b  <- mask_nonsig_lfc(megreen_deseq)
meblue_deseq_b   <- mask_nonsig_lfc(meblue_deseq)

# NOTE: purplecaz currently has no direction/significance filter, unlike the
# other three below — add one (matching the pattern used for browncaz) before
# including purple in the combined heatmap, or it isn't really "GFCA-specific"
# by the same definition the other categories use.
purplecaz <- mepurple_deseq_b %>%
  filter(!is.na(CAZ_anno), kME >= 0.7) %>%
  mutate(Name = paste0(Protein_name, " (proteinID = ", ProteinID, ")"))

browncaz <- mebrown_deseq_b %>%
  filter(!is.na(CAZ_anno), kME >= 0.7, GFCA_v_blank_log2FoldChange > 0, GFCA_v_blank_padj < 0.05) %>%
  mutate(Name = paste0(Protein_name, " (proteinID = ", ProteinID, ")"))

blackcaz <- meblack_deseq_b %>%
  filter(!is.na(CAZ_anno), kME >= 0.7, GFCA_v_blank_log2FoldChange > 0, GFCA_v_blank_padj < 0.05) %>%
  mutate(Name = paste0(Protein_name, " (proteinID = ", ProteinID, ")"))

greencaz <- megreen_deseq_b %>%
  filter(!is.na(CAZ_anno), kME >= 0.7, cellulose_v_blank_log2FoldChange > 0, cellulose_v_blank_padj < 0.05) %>%
  mutate(Name = paste0(Protein_name, " (proteinID = ", ProteinID, ")"))

bluecaz <- meblue_deseq_b %>%
  filter(!is.na(CAZ_anno), kME >= 0.7,
         GFCA_v_blank_log2FoldChange > 0, GFCA_v_blank_padj < 0.05,
         cellulose_v_blank_log2FoldChange > 0, cellulose_v_blank_padj < 0.05) %>%
  mutate(Name = paste0(Protein_name, " (proteinID = ", ProteinID, ")"))

# =============================================================================
# 10. Combined heatmap (one figure, dendrogram, category shown via row annotation)
# =============================================================================

all_caz <- bind_rows(
  browncaz %>% mutate(category = "GF CA-specific"),
  greencaz %>% mutate(category = "Cellulose-specific"),
  bluecaz  %>% mutate(category = "Shared (GF CA + Cellulose)")
  # , purplecaz %>% mutate(category = "GF CA-specific")  # add a direction/significance filter to purplecaz first (see note above)
)

table(all_caz$category)

# Row label: "Protein name (CAZ family(ies); ProteinID)"
# e.g. one CAZ family:      "Alpha-glucosidase (GH31_1; 6288750)"
#      multiple CAZ families: "... (GH31_1, GH2; 6288750)"
all_caz$Name <- paste0(all_caz$Protein_name, " (", gsub("; ", ", ", all_caz$CAZ_anno), "; ", all_caz$ProteinID, ")")

# ── Fix URL-encoded characters and known JGI annotation typos ─────────────────
all_caz$Name <- gsub("%27", "'", all_caz$Name)
all_caz$Name <- gsub("%2C", ",", all_caz$Name)
all_caz$Name <- gsub("interation domain-domain", "interaction domain", all_caz$Name)

# ── Capitalize just the first letter of each Name (leave the rest as-is) ──────
all_caz$Name <- paste0(toupper(substr(all_caz$Name, 1, 1)), substr(all_caz$Name, 2, nchar(all_caz$Name)))

all_caz <- as.data.frame(all_caz)
rownames(all_caz) <- all_caz$Name

lfc_cols <- c("SACA_v_blank_log2FoldChange", "GFCA_v_blank_log2FoldChange",
              "cellulose_v_blank_log2FoldChange", "GFCA_v_cellulose_log2FoldChange")

newnames <- c("SA CA vs.\nNo substrate", "GF CA vs.\nNo substrate",
              "Cellulose vs.\nNo substrate", "GF CA vs.\nCellulose")

annotation_row <- all_caz %>%
  dplyr::select(category) %>%
  dplyr::rename(Category = category)

ann_colors <- list(
  Category = c(
    "GF CA-specific"             = "saddlebrown",
    "Shared (GF CA + Cellulose)" = "turquoise1",
    "Cellulose-specific"        = "green"
  )
)

# ── Custom pheatmap column-label renderer (horizontal, slightly raised) ───────
draw_colnames_horizontal <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"),
           vjust = 1.1, hjust = 0.5, rot = 0, gp = gpar(...))
}
assignInNamespace("draw_colnames", draw_colnames_horizontal, ns = asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(299)

shared_max <- max(abs(as.matrix(all_caz[, lfc_cols])), na.rm = TRUE)
shared_breaks <- seq(-shared_max, shared_max, length.out = 300)

specificity_hm <- pheatmap(
  data.matrix(all_caz[, lfc_cols]),
  col                  = my_palette,
  cellwidth            = 80,
  cellheight           = 12,
  cluster_rows         = TRUE,   # dendrogram groups genes by expression similarity
  cluster_cols         = FALSE,
  fontsize             = 12,
  fontsize_col         = 12,
  fontsize_row         = 10,
  show_rownames        = TRUE,
  labels_col           = newnames,
  border_color         = "black",
  breaks               = seq(-5, 5, length.out = 300),
  annotation_row       = annotation_row,
  annotation_colors    = ann_colors,
  annotation_names_row = FALSE
)

specificity_hm$gtable$heights[length(specificity_hm$gtable$heights)] <- unit(1.2, "cm")
grid.newpage()
grid.draw(specificity_hm$gtable)

ggsave(plot = specificity_hm, filename = file.path(FIG_DIR, paste0("Pchry_wgcna_GFCAcellulose_heatmap_", DATE, ".pdf")),
       width = 16.5, height = 10, units = "in", dpi = 300)
ggsave(plot = specificity_hm, filename = file.path(FIG_DIR, paste0("Pchry_wgcna_GFCAcellulose_heatmap_", DATE, ".tiff")),
       width = 16.5, height = 10, units = "in", dpi = 300, bg = "white")

message("Done.")


# Combine MDS plot (A) and GFCA/Cellulose specificity heatmap (B) into one figure in powerpoint

# Total significantly differentially expressed genes per comparison (padj < 0.05)

n_sig_gfca_v_blank <- sum(vary_43caz$GFCA_v_blank_padj < 0.05, na.rm = TRUE)
n_sig_cellulose_v_blank <- sum(vary_43caz$cellulose_v_blank_padj < 0.05, na.rm = TRUE)
n_sig_gfca_v_cellulose <- sum(vary_43caz$GFCA_v_cellulose_padj < 0.05, na.rm = TRUE)

message("GF CA vs. blank: ", n_sig_gfca_v_blank, " significant DEGs")
message("Cellulose vs. blank: ", n_sig_cellulose_v_blank, " significant DEGs")
message("GF CA vs. Cellulose: ", n_sig_gfca_v_cellulose, " significant DEGs")

step1 <- all_caz %>%
  dplyr::select(category, CAZ_anno) %>%
  filter(!is.na(CAZ_anno))

step2 <- step1 %>%
  separate_rows(CAZ_anno, sep = ";\\s*")

str(step2$CAZ_anno)
str(step2$category)

step3 <- step2 %>%
  mutate(CAZ_class = str_extract(CAZ_anno, "^[A-Z]+"))

str(step3$CAZ_class)
class(step3)

caz_breakdown <- step3 %>%
  dplyr::count(category, CAZ_class) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)

caz_breakdown
