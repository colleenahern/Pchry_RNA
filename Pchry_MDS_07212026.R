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
# 4. Filter low-expression genes using filterByExpr
# =============================================================================

# Filter genes
keep <- filterByExpr(pcy, group = group)
summary(keep)
pcy <- pcy[keep, , keep.lib.sizes = FALSE]

# =============================================================================
# 5. Library size and normalization diagnostics
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
# 6. MDS plot
# =============================================================================
# NOTE: confirm the substrate-to-color/label mapping below matches your actual
# groupName factor levels (check with `levels(group)`) before trusting the legend.

levels(group)

col.cell <- c("blue", "green", "purple", "orange")[group]
data.frame(group = pcsampleinfo$groupName, col.cell)

# ── Reusable plotting function so TIFF and PDF versions stay identical ────────
draw_mds_plot <- function() {
  mdsplot <- plotMDS(pcy, col = col.cell, pch = 16)
  legend("top", fill = c("blue", "green", "purple", "orange"),
         legend = c("No substrate", "Cellulose", "Lowly-acetylated CA (GFCA)", "Highly-acetylated CA (SACA)"),
         title = expression(bold("Substrate")))
  invisible(mdsplot)
}

# ── TIFF ────────────────────────────────────────────────────────────────────
tiff(filename = file.path(BASE, "Figures/Pchry_MDS_07212026.tiff"),
     units = "in", width = 7, height = 7, res = 300)
draw_mds_plot()
dev.off()

# ── PDF ─────────────────────────────────────────────────────────────────────
pdf(file = file.path(BASE, "Figures/Pchry_MDS_07212026.pdf"),
    width = 7, height = 7)
draw_mds_plot()
dev.off()

message("Done. MDS plot saved to:\n  ",
        file.path(BASE, "Figures/Pchry_MDS_07212026.tiff"), "\n  ",
        file.path(BASE, "Figures/Pchry_MDS_07212026.pdf"))
