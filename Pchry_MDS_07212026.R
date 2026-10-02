# Colleen Ahern
# 07/21/2026 
# MDS plot from Pchry featureCounts count matrix

# ── Libraries ─────────────────────────────────────────────────────────────────

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(NMF)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE <- "/Users/colleenahern/Documents/Pchry/DESeq2_analysis"

# =============================================================================
# 1. Import count data and sample info
# =============================================================================

pcseqdata    <- read.delim(file.path(BASE, "featcount/pchry_countmat_jgi_07312025.txt"),
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

# =============================================================================
# 2. Build DGEList
# =============================================================================

pcy <- DGEList(pccountdata)

group <- factor(pcsampleinfo$groupName)
pcy$samples$group <- group
pcy$samples

# =============================================================================
# 3. Filter low-expression genes using filterByExpr
# =============================================================================

# Filter genes
keep <- filterByExpr(pcy, group = group)
summary(keep)
pcy <- pcy[keep, , keep.lib.sizes = FALSE]

# =============================================================================
# 4. Library size and normalization diagnostics
# =============================================================================

pcy$samples$lib.size

barplot(pcy$samples$lib.size / 1e06, names = colnames(pcy), las = 2, ann = FALSE, cex.names = 0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# ── Pre-normalization logCPM distribution ──────────────────────────────────────
logcounts <- cpm(pcy, log = TRUE)
boxplot(logcounts, xlab = "", ylab = "Log2 counts per million", las = 2)
abline(h = median(logcounts), col = "blue")
title("Boxplots of logCPMs (unnormalised)")

# ── Apply TMM normalization ────────────────────────────────────────────────────
pcy <- calcNormFactors(pcy)
pcy$samples

# ── Post-normalization logCPM distribution ─────────────────────────────────────
logcounts <- cpm(pcy, log = TRUE)
boxplot(logcounts, xlab = "", ylab = "Log2 counts per million", las = 2)
abline(h = median(logcounts), col = "blue")
title("Boxplots of logCPMs (normalised)")

# =============================================================================
# 5. MDS plot
# =============================================================================
# NOTE: confirm the substrate-to-color/label mapping below matches your actual
# groupName factor levels (check with `levels(group)`) before trusting the legend.

levels(group)

col.cell <- c("grey40", "#1B9E77", "#D95F02", "#E7298A")[group]
data.frame(group = pcsampleinfo$groupName, col.cell)

# ── Reusable plotting function so TIFF and PDF versions stay identical ────────
draw_mds_plot <- function() {
  mdsplot <- plotMDS(pcy, col = col.cell, pch = 16)
  legend("top", fill = c("grey40", "#1B9E77", "#D95F02", "#E7298A"),
         legend = c("No substrate", "Cellulose", "GF CA", "SA CA"),
         title = expression(bold("Substrate")))
  invisible(mdsplot)
}

# ── TIFF ────────────────────────────────────────────────────────────────────
tiff("/Users/colleenahern/Documents/Pchry/Figures/Pchry_MDS_08232026.tiff",
     units = "in", width = 7, height = 7, res = 300)
draw_mds_plot()
dev.off()

# ── PDF ─────────────────────────────────────────────────────────────────────
pdf("/Users/colleenahern/Documents/Pchry/Figures/Pchry_MDS_08232026.pdf",
    width = 7, height = 7)
draw_mds_plot()
dev.off()





# Build mds_plot_gg from the pcy
# object already created above, matching your existing colors/labels.

library(ggplot2)

# ── Get MDS coordinates without plotting (plot = FALSE just returns the data) ──
mds_data <- plotMDS(pcy, plot = FALSE)

mds_df <- data.frame(
  x = mds_data$x,
  y = mds_data$y,
  group = pcsampleinfo$groupName
)

# ── Match your existing group levels (check this matches `levels(group)` from your script) ──
mds_df$group <- factor(mds_df$group, levels = levels(group))

mds_plot_gg <- ggplot(mds_df, aes(x = x, y = y, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("grey40", "#1B9E77", "#D95F02", "#E7298A"),
    labels = c("No substrate", "Cellulose", "GF CA", "SA CA")
  ) +
  labs(
    x = paste0("Leading logFC dim 1 (", round(mds_data$var.explained[1] * 100), "%)"),
    y = paste0("Leading logFC dim 2 (", round(mds_data$var.explained[2] * 100), "%)"),
    color = "Substrate"
  ) +
  theme_classic(base_size = 12) +
  theme(
    text                  = element_text(size = 12),
    axis.title            = element_text(size = 12),
    axis.text             = element_text(size = 12),
    legend.title          = element_text(size = 12, face = "bold"),
    legend.text           = element_text(size = 12),
    legend.key.size       = unit(0.4, "cm"),
    legend.position       = c(0.5, 0.88),  # inset into the blank middle-top area — adjust if it overlaps points
    legend.justification  = c(0.5, 0.5),
    # legend.background     = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    legend.margin         = margin(5, 8, 5, 8),
    panel.border          = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line             = element_blank(),
    aspect.ratio          = 1,
    plot.margin           = margin(5, 5, 5, 5)
  )

mds_plot_gg

ggsave(plot = mds_plot_gg,
       filename = file.path(FIG_DIR, paste0("Pchry_MDSgg_", DATE, ".pdf")),
       width = 5, height = 5, units = "in", dpi = 300)
å