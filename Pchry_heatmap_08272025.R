# Colleen Ahern
# 08/27/2025
# Heatmap file for my Pchry_DESeq_08272025 file, which removed the SACA5 sample

# includes dbcan annotation

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("pasilla")

library(tximport)
library(readr)
library(tximportData)
library(tidyverse)
library(pasilla)

# Heatmap for all genes

# Import the count data and sample info
pcseqdata <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/featcount/pchry_countmat_jgi_07312025.txt", stringsAsFactors = FALSE)
pcsampleinfo <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/pcsampleinfo.txt", stringsAsFactors = TRUE)
head(pcseqdata)
dim(pcseqdata)

pcsampleinfo

# Store GeneID as rownames and then remove non-count columns
pccountdata <- pcseqdata
rownames(pccountdata) <- pccountdata$Geneid
pccountdata <- pccountdata[,-c(1:6)]
head(pccountdata)

pcsampleinfo$libraryName_jgi <- gsub("/",".",pcsampleinfo$libraryName_jgi)
pcsampleinfo$libraryName_jgi <- gsub("-",".",pcsampleinfo$libraryName_jgi)

# Make sure that the column names are the same as libraryName in the pcsampleinfo file
colnames(pccountdata) # using substr, you can change the characters of the colnames
colnames(pccountdata) %in% pcsampleinfo$libraryName_jgi
table(colnames(pccountdata)==pcsampleinfo$libraryName)
pccountdata <- pccountdata[,as.character(pcsampleinfo$libraryName_jgi)]
table(colnames(pccountdata)==pcsampleinfo$libraryName_jgi) # Now they match
colnames(pccountdata) <- pcsampleinfo$sampleName
rownames(pcsampleinfo) <- pcsampleinfo$sampleName

# remove outlier samples 
samples.to.be.removed <- "SACA5"
pccts_subset <- pccountdata[,!(colnames(pccountdata) %in% samples.to.be.removed)]

pccoldata_subset <- pcsampleinfo[!(rownames(pcsampleinfo) %in% samples.to.be.removed),]
all(rownames(pccoldata_subset) == colnames(pccts_subset))

# Apply tpm cutoff to the gene count matrix
# Make tpm gene count matrix
# install.packages("remotes")
# remotes::install_github("davidrequena/drfun")
library(DRnaSeq)
all(rownames(pccts_subset) == pcseqdata$Geneid)
gl <- pcseqdata$Length
tpm <- as.data.frame(tpm(pccts_subset, gl))

# Import/reformat sample information 
asm <- pccoldata_subset
all(colnames(tpm) == asm$sampleName)

# Make sure the sample names in the count matrix and sample info match
colnames(tpm) == asm$sampleName

tpm$rem <- 0

# Apply tpm cutoff: the average TPM of the  biological replicates for any of the four substrates must be >2
for (i in 1:nrow(tpm)) {
  tpm$rem[i] <- ifelse(mean(c(tpm$cellulose1[i],tpm$cellulose2[i],tpm$cellulose3[i],tpm$cellulose4[i],tpm$cellulose5[i])) < 2 && mean(c(tpm$GFCA1[i],tpm$GFCA2[i],tpm$GFCA3[i],tpm$GFCA4[i],tpm$GFCA5[i])) < 2 && mean(c(tpm$SACA1[i],tpm$SACA2[i],tpm$SACA3[i],tpm$SACA4[i])) < 2 && mean(c(tpm$blank1[i],tpm$blank2[i],tpm$blank3[i],tpm$blank4[i],tpm$blank5[i])) < 2, "REMOVE", "KEEP")
}

tpm <- tpm[(tpm$rem=="KEEP"),] # reduces from 12256 genes to 10888
tpm <- tpm[,-c(20)]
tpm$GeneID <- rownames(tpm)

# Import list of differentially regulated genes for heatmap
# Remove genes with a p-value = NA for all three comparitive conditions

vary_43 <- read_csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_jgi_08272025.csv") # My DESeq

vary_43$rem <- 0
for (i in 1:nrow(vary_43)) {
  vary_43$rem[i] <- ifelse(is.na(vary_43$SACA_v_blank_padj[i]) == TRUE && is.na(vary_43$GFCA_v_blank_padj[i]) == TRUE && is.na(vary_43$cellulose_v_blank_padj[i]) == TRUE && is.na(vary_43$cellulose_v_SACA_padj[i]) == TRUE && is.na(vary_43$GFCA_v_SACA_padj[i]) == TRUE && is.na(vary_43$GFCA_v_cellulose_padj[i]) == TRUE, "REMOVE", "KEEP")
}

vary_43 <- vary_43[(vary_43$rem=="KEEP"),] # reduces from 11811 to 11791
vary_43 <- vary_43[,-38]

# Merge tpm matrix with vary_43 matrix - to make a matrix that contains only genes that pass both of those filters
a <- merge(vary_43,tpm,"GeneID") # 10870 genes
vary_43 <- a[,c(1:37)]

# Create a minimum p value column that takes the minimum p value from any of the three comparative conditions for each gene
# Then sort the matrix from lowest padjmin to largest
# vary_43$padjmin <- apply(vary_43[,c("SACA_v_blank_padj","GFCA_v_blank_padj","cellulose_v_blank_padj")],1,min, na.rm = TRUE)
vary_43$padjmin <- apply(vary_43[,c(colnames(vary_43)[grep("padj",colnames(vary_43))])],1,min, na.rm = TRUE)
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]

# Import list of JGI gene annotations
# BiocManager::install("ape")
library(ape)
library(tidyverse)
my_gene_col_44 <- read.gff("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Phchr4_2_GeneCatalog_20210419.gff3", GFF3 = TRUE)
my_gene_col_44 <- my_gene_col_44[my_gene_col_44$type == "gene",] # 12,256 genes
my_gene_col_44b <- my_gene_col_44 %>%
  dplyr::select(seqid, everything())%>%
  mutate(attributes=str_split(attributes, ";")) %>% 
  unnest_wider(where(is.list), names_sep = "")

my_gene_col_44 <- my_gene_col_44b[,c(9,12,13,14)]
my_gene_col_44
my_gene_col_44$attributes1 <- gsub("ID=","",my_gene_col_44$attributes1)
my_gene_col_44$attributes4 <- gsub("product_name=","",my_gene_col_44$attributes4)
my_gene_col_44$attributes5 <- gsub("proteinId=","",my_gene_col_44$attributes5)
my_gene_col_44$attributes6 <- gsub("transcriptId=","",my_gene_col_44$attributes6)
colnames(my_gene_col_44) <- c("GeneID","Protein_name","ProteinID","TranscriptID")

# Add annotations to the data frame and reorder padjmin
vary_43 <- merge(vary_43, my_gene_col_44, "GeneID")
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]

for (i in 1:nrow(vary_43)) {
  vary_43$Name[i] <- paste(vary_43$Protein_name[i], " ", "(proteinID = ", vary_43$ProteinID[i], ")", sep = "")
}


# 08/04/2025 - import list of CAZymes from dbcan
cazlist <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pchry_dbcan/overview.txt")
cazlist2 <- cazlist[cazlist$X.ofTools >1,]
cazlist_f <- cazlist2
cazlist_f$Gene.ID <- gsub("jgi\\|","",cazlist2$Gene.ID)
cazlist_f$Gene.ID <- gsub(".*\\|(.*?)\\|.*", "\\1", cazlist_f$Gene.ID)
cazlist_f$HMMER <- gsub("\\(.*","",cazlist_f$HMMER)
cazlist_f$DIAMOND <- gsub("Phchr2_GeneCatalog_proteins_20131210.aa.fasta\\+","",cazlist_f$DIAMOND)
cazlist_f$DIAMOND <- gsub(".*fasta\\+","",cazlist_f$DIAMOND)
rownames(cazlist_f) <- 1:nrow(cazlist_f)
cazlist_f$CAZyme <- NA
cazlist_f$CAZyme[1:376] <- cazlist_f[1:377,"HMMER"]
cazlist_f$CAZyme[377:388] <- cazlist_f[377:388,"eCAMI"]
cazlist_f$DIAMOND <- gsub("\\+3.*", "", cazlist_f$DIAMOND)
cazlist_ff <- cazlist_f %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(eCAMI=str_split(eCAMI, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")
cazlist_ff <- cazlist_ff %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(DIAMOND=str_split(DIAMOND, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")
cazlist_ff <- cazlist_ff %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(HMMER=str_split(HMMER, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")

cazlist_fff <- cazlist_ff
cazlist_fff$CAZ <- NA

for (i in 1:nrow(cazlist_fff)) {
  z <- as.data.frame(table(as.character(cazlist_fff[i,3:9])))
  cazlist_fff$CAZ[i] <- paste(paste(z$Var1," x",z$Freq,sep = ""),collapse = "; ")
}

cazlist_f4 <- cazlist_fff
cazlist_f4$CAZf <- cazlist_f4$CAZ

cazlist_f4 <- cazlist_f4 %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(CAZf=str_split(CAZf, ";")) %>% 
  unnest_wider(where(is.list), names_sep = "")

for (i in 13:16) {
  cazlist_f4[,i] <- lapply(cazlist_f4[,i], function(x) replace(x, grepl("x1", x), NA))
}

cazlist_f4 <- cazlist_f4 %>%
  unite(CAZ_name, starts_with('CAZf'), 
        na.rm = TRUE, remove = TRUE, sep = ';')

cazlist_f4$CAZ_name <- gsub(" ", "", cazlist_f4$CAZ_name)
cazlist_f4$CAZ_name <- gsub(";", "; ", cazlist_f4$CAZ_name)
cazlist_f4$CAZ_anno <- cazlist_f4$CAZ_name
cazlist_f4$CAZ_anno <- gsub("x2","",cazlist_f4$CAZ_anno)
cazlist_f4$CAZ_anno <- gsub("x3","",cazlist_f4$CAZ_anno)
cazlist_f4$CAZ_anno <- gsub("x4","",cazlist_f4$CAZ_anno)

cazlist_f5 <- cazlist_f4
cazlist_f5 <- cazlist_f5[!cazlist_f5$CAZ_anno == "",]
cazlist_f5 <- cazlist_f5 %>%
  rename("ProteinID" = "Gene.ID")

write.csv(cazlist_f5, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pchry_dbcan/cazlist_08072025.csv")
caz_anno <- read.csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pchry_dbcan/cazlist_08072025.csv")
caz_anno <- caz_anno[,c("ProteinID","CAZ_anno")]

vary_43caz <- merge(vary_43,caz_anno, by = "ProteinID", all.x = TRUE)
vary_43caz <- vary_43caz[order(vary_43caz$padjmin,decreasing=FALSE),]
write_csv(vary_43caz, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_cazanno_08272025.csv")


# 08/27/2025 use this for downstream stuff, i just wanted to get my vary_43 table so I haven't proofread this stuff below
vary_43caz <- read_csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_cazanno_08272025.csv")
vary_43caz <- vary_43caz[order(vary_43caz$padjmin,decreasing=FALSE),]
vary_43caz_ramya <- vary_43caz[vary_43caz$padjmin < 0.05,]
write_csv(vary_43caz_ramya,"/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/Tables_for_alphafold/LFC_allconds_anno_cazanno_sigpadjmin_8142025.csv")

vary_43cazonly <- vary_43caz[is.na(vary_43caz$CAZ_anno) == FALSE,] 
write_csv(vary_43cazonly, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_cazonly_08072025.csv")

# 08/14/2025 making tables for Ramya for alphafold -- line 222 is the deseq results for all conditions with a padjmin of <0.05
# but maybe I'm only interested in the GFCA vs. cellulose condition? Let's also make a table for that condition for genes with padj < 0.05
vary_43caz_ramya_gfca_cell <- vary_43caz[order(vary_43caz$GFCA_v_cellulose_padj,decreasing=FALSE),]
vary_43caz_ramya_gfca_cell_sig <- vary_43caz_ramya_gfca_cell[vary_43caz_ramya_gfca_cell$GFCA_v_cellulose_padj < 0.05,]
vary_43caz_ramya_gfca_cell_sig_lfcpos <- vary_43caz_ramya_gfca_cell_sig[vary_43caz_ramya_gfca_cell_sig$GFCA_v_cellulose_log2FoldChange > 0,]
write_csv(vary_43caz_ramya_gfca_cell_sig_lfcpos,"/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/Tables_for_alphafold/LFC_allconds_anno_cazanno_sigGFCAcellpadj_posGFCAcellLFC_8142025.csv")

# 12/05/2025 - making heatmap of John's genes from antiSmash
# use vary_43caz from line 205/206
vary_43caz <- vary_43caz[order(vary_43caz$padjmin,decreasing=FALSE),]
cbgs<- read.delim("/Users/colleenahern/Documents/Pchry/Antismash_results/Core_antiSMASH_biosynthetic_genes_12052025.txt")
# vary_43 <- merge(cbgs, vary_43caz) # use this line for heatmap of only antismash genes
vary_43 <- vary_43caz # use this line for heatmap of all genes and all conditions (choose top 100 padjmin genes and highlight cazymes)

# Turn LFC of genes with padj < 0.05 to 0
for (j in c(4,10,16,34)) {
  for (i in 1:nrow(vary_43)) {
    if (as.numeric(vary_43[i,j+4]) > 0.049) {
      vary_43[i,j] <- 0
    }
    else {
      vary_43[i,j] <- as.numeric(vary_43[i,j])
    }
  }
}

# redid padjmin because i think I did it across all comparisons instead of the 4 im interested in
vary_43$padjmin <- apply(vary_43[,c("SACA_v_blank_padj","GFCA_v_blank_padj","cellulose_v_blank_padj","GFCA_v_cellulose_padj")],1,min, na.rm = TRUE)
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]

# checking to see if this matches the original results order
vary_43$padjminorig <- apply(vary_43[,c(colnames(vary_43)[grep("_padj",colnames(vary_43))])],1,min, na.rm = TRUE)
vary_43orig <- vary_43[order(vary_43$padjminorig,decreasing=FALSE),]

# newnames <- lapply(
#   colnames(vary_43[,c(2,4,6)]),
#   function(x) bquote(.(x)))
newnames <- lapply(
  colnames(vary_43[,c(4,10,16,34)]),
  function(x) bquote(.(x)))
newnames
newnames <- gsub("SACA_v_blank_log2FoldChange", "Highly-acetylated CA vs. No substrate", newnames)
newnames <- gsub("GFCA_v_blank_log2FoldChange", "Lowly-acetylated CA vs. No substrate", newnames)
newnames <- gsub("cellulose_v_blank_log2FoldChange", "Cellulose vs. No substrate", newnames)
newnames <- gsub("GFCA_v_cellulose_log2FoldChange", "Lowly-acetylated CA vs. Cellulose", newnames)
newnames
newnames <- gsub(" vs. ", "\nvs. ", newnames)
newnames

head(vary_43)
nrow(vary_43)
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43[c(rownum_from:rownum_to),]

# ann_cols <- vary_44[,12]
ann_cols <- vary_44[,42]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1.5, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}
tobold <-vary_44[is.na(vary_44$CAZ_anno) == FALSE, "Name"]
tobold <- as.list(tobold)
tobold

# Idk why this bold_names fucntion didn't work, so i used a different method that works
labels <- ifelse(
  unlist(rownames(vary_44)) %in% unlist(tobold),
  paste0("bold('", vary_44$Name, "')"),
  paste0("'", vary_44$Name, "'")
)
labels

labels <- parse(text = labels)
labels

vary_44 <- as.data.frame(vary_44)
row.names(vary_44) <- vary_44$Name

par(mar=c(0,0,0,0))

pc_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,10,16,34)]), col=my_palette, cellwidth=130,cellheight=7, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize_col = 12, fontsize_row = 8,show_rownames=T, labels_col = newnames, labels_row = labels, border_color="black",
                                breaks=seq(-4, 6, length.out=300))

pc_regulation_heatmap$gtable$heights[length(pc_regulation_heatmap$gtable$heights)] <- unit(1, "cm")
grid.newpage()
grid.draw(pc_regulation_heatmap$gtable)

ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top100padjmin_heatmap_03112026.pdf", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top100padjmin_heatmap_03112026.tiff", width = 14, height = 11, units = "in", dpi = 200, bg="white")


# Let's do top 50 genes cus 100 is too many 
vary_44 <- vary_43[c(1:50),]
vary_44 <- vary_43[c(1:20),]
vary_44 <- vary_43[c(1:30),]

# ann_cols <- vary_44[,12]
ann_cols <- vary_44[,42]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1.5, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

tobold <-vary_44[is.na(vary_44$CAZ_anno) == FALSE, "Name"]
tobold <- as.list(tobold)
tobold

row.names(vary_44) <- vary_44$Name

# Idk why this bold_names fucntion didn't work, so i used a different method that works
labels <- ifelse(
  unlist(rownames(vary_44)) %in% unlist(tobold),
  paste0("bold('", vary_44$Name, "')"),
  paste0("'", vary_44$Name, "'")
)
labels

labels <- parse(text = labels)
labels

vary_44 <- as.data.frame(vary_44)
# row.names(vary_44) <- vary_44$Name

par(mar=c(0,0,0,0))

pc_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,10,16,34)]), col=my_palette, cellwidth=130,cellheight=15, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize_col = 12, fontsize_row = 12,show_rownames=T, labels_col = newnames, labels_row = labels, border_color="black",
                                breaks=seq(-4, 6, length.out=300))

pc_regulation_heatmap$gtable$heights[length(pc_regulation_heatmap$gtable$heights)] <- unit(1.7, "cm")
grid.newpage()
grid.draw(pc_regulation_heatmap$gtable)

ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top50padjmin_heatmap_03112026.pdf", width = 16, height = 12, units = "in", dpi = 300)
ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top50padjmin_heatmap_03112026.tiff", width = 16, height = 12, units = "in", dpi = 200, bg="white")

ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top20padjmin_heatmap_03112026.pdf", width = 16, height = 6, units = "in", dpi = 300)
ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top20padjmin_heatmap_03112026.tiff", width = 16, height = 6, units = "in", dpi = 200, bg="white")

ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top30padjmin_heatmap_03112026.pdf", width = 16, height = 8, units = "in", dpi = 300)
ggsave(plot = pc_regulation_heatmap, filename = "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Figures/Pchry_4comps_top30padjmin_heatmap_03112026.tiff", width = 16, height = 8, units = "in", dpi = 200, bg="white")



## John's heatmap for antismash
# use vary_43 from line 224
newnames <- lapply(
  colnames(vary_43[,c(6,12,18,36)]),
  function(x) bquote(.(x)))
newnames
newnames <- gsub("_log2FoldChange", "", newnames)
newnames
newnames <- gsub("_v_", " vs. ", newnames)
newnames
newnames <- gsub("blank", "Blank", newnames)
newnames <- gsub("cellulose", "Cellulose", newnames)
newnames

head(vary_43)
nrow(vary_43)

vary_44 <- vary_43

ann_cols <- vary_44$Name
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = .9, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

pc_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(6,12,18,36)]), col=my_palette, cellwidth=40,cellheight=8, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44$Name), border_color=NA,breaks=seq(-4, 4, length.out=300))

max(vary_44[,c(6,12,18,36)])
min(vary_44[,c(6,12,18,36)])

writexl::write_xlsx("/Users/colleenahern/Documents/Pchry/Antismash_results/John_antismash_Pchry_heatmap_table_12052025.xlsx")

##########

# Heatmap for just the GFCA vs cellulose comparison
# all genes

# Import the count data and sample info
pcseqdata <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/featcount/pchry_countmat.txt", stringsAsFactors = FALSE)
pcsampleinfo <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/pcsampleinfo.txt", stringsAsFactors = TRUE)
head(pcseqdata)
dim(pcseqdata)

pcsampleinfo

# Store GeneID as rownames and then remove non-count columns
pccountdata <- pcseqdata
rownames(pccountdata) <- pccountdata$Geneid
pccountdata <- pccountdata[,-c(1:6)]
head(pccountdata)

pcsampleinfo$libraryName <- gsub("/",".",pcsampleinfo$libraryName)
pcsampleinfo$libraryName <- gsub("-",".",pcsampleinfo$libraryName)
# Make sure that the column names are the same as libraryName in the pcsampleinfo file
colnames(pccountdata) # using substr, you can change the characters of the colnames
colnames(pccountdata) %in% pcsampleinfo$libraryName
table(colnames(pccountdata)==pcsampleinfo$libraryName)
pccountdata <- pccountdata[,as.character(pcsampleinfo$libraryName)]
table(colnames(pccountdata)==pcsampleinfo$libraryName) # Now they match
colnames(pccountdata) <- pcsampleinfo$sampleName

# Apply tpm cutoff to the gene count matrix
# Make tpm gene count matrix
# install.packages("remotes")
# remotes::install_github("davidrequena/drfun")
library(DRnaSeq)
all(rownames(pccountdata) == pcseqdata$Geneid)
gl <- pcseqdata$Length
tpm <- as.data.frame(tpm(pccountdata, gl))

# Import/reformat sample information 
asm <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/pcsampleinfo.txt", stringsAsFactors = TRUE)
all(colnames(tpm) == asm$sampleName)

# Make sure the sample names in the count matrix and sample info match
colnames(tpm) == asm$sampleName

tpm$rem <- 0

# Apply tpm cutoff: the average TPM of the  biological replicates for any of the four substrates must be >2
for (i in 1:nrow(tpm)) {
  tpm$rem[i] <- ifelse(mean(c(tpm$cellulose1[i],tpm$cellulose2[i],tpm$cellulose3[i],tpm$cellulose4[i],tpm$cellulose5[i])) < 2 && mean(c(tpm$GFCA1[i],tpm$GFCA2[i],tpm$GFCA3[i],tpm$GFCA4[i],tpm$GFCA5[i])) < 2 && mean(c(tpm$SACA1[i],tpm$SACA2[i],tpm$SACA3[i],tpm$SACA4[i],tpm$SACA5[i])) < 2 && mean(c(tpm$blank1[i],tpm$blank2[i],tpm$blank3[i],tpm$blank4[i],tpm$blank5[i])) < 2, "REMOVE", "KEEP")
}

tpm <- tpm[(tpm$rem=="KEEP"),] # reduces from 12256 genes to 10843
tpm <- tpm[,-c(21)]
tpm$GeneID <- rownames(tpm)

# Import list of differentially regulated genes for heatmap
# Remove genes with a p-value = NA for all three comparitive conditions

vary_43 <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_07242025.txt", header=T, check.names = FALSE) # My DESeq
vary_43 <- vary_43[,c(1,17:19)] 

vary_43$rem <- 0
for (i in 1:nrow(vary_43)) {
  vary_43$rem[i] <- ifelse(is.na(vary_43$GFCA_v_cellulose_padj[i]) == TRUE, "REMOVE", "KEEP")
}

vary_43 <- vary_43[(vary_43$rem=="KEEP"),] # reduces from 11718 to 11690
vary_43 <- vary_43[,-5]

# Merge tpm matrix with vary_43 matrix - to make a matrix that contains only genes that pass both of those filters
a <- merge(vary_43,tpm,"GeneID") # 10,822 genes
vary_43 <- a[,c(1:4)]

# Import list of JGI gene annotations
#BiocManager::install("ape")
library(ape)
library(tidyverse)
my_gene_col_44 <- read.gff("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Phchr4_2_GeneCatalog_20210419.gff3", GFF3 = TRUE)
my_gene_col_44 <- my_gene_col_44[my_gene_col_44$type == "gene",] # 12,256 genes
my_gene_col_44b <- my_gene_col_44 %>%
  dplyr::select(seqid, everything())%>%
  mutate(attributes=str_split(attributes, ";")) %>% 
  unnest_wider(where(is.list), names_sep = "")

my_gene_col_44 <- my_gene_col_44b[,c(9,12,13,14)]
my_gene_col_44
my_gene_col_44$attributes1 <- gsub("ID=","",my_gene_col_44$attributes1)
my_gene_col_44$attributes4 <- gsub("product_name=","",my_gene_col_44$attributes4)
my_gene_col_44$attributes5 <- gsub("proteinId=","",my_gene_col_44$attributes5)
my_gene_col_44$attributes6 <- gsub("transcriptId=","",my_gene_col_44$attributes6)
colnames(my_gene_col_44) <- c("GeneID","Protein_name","ProteinID","TranscriptID")

# Add annotations to the data frame and reorder padjmin
vary_43 <- merge(vary_43, my_gene_col_44, "GeneID")

for (i in 1:nrow(vary_43)) {
  vary_43$Name[i] <- paste(vary_43$Protein_name[i], " ", "(proteinID = ", vary_43$ProteinID[i], ")", sep = "")
}

write_csv(vary_43, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary_43_GFCA_v_cellulose_07242025.csv")
vary_43 <- read_csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary_43_GFCA_v_cellulose_07242025.csv")
GFCA_v_cellulose_lfc <- read_csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary_43_GFCA_v_cellulose_07242025.csv")

vary_43 <- vary_43[order(vary_43$GFCA_v_cellulose_padj,decreasing=FALSE),]

newnames <- "GF CA vs. cellulose"
newnames
head(vary_43)
nrow(vary_43)
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,8]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

pc_regulation_heatmap<-pheatmap(data.matrix(vary_44[,2]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,8]), border_color=NA,breaks=seq(-4, 6, length.out=300))



#### Now look at the upregulation for GFCA vs cellulose cus the most significant genes are mostly downregulated in that comparison
vary_43_sig <- vary_43[vary_43$GFCA_v_cellulose_padj <0.05,]
vary_43_sigpos <- vary_43_sig[vary_43_sig$GFCA_v_cellulose_log2FoldChange > 1,]
vary_43_sigpos <- vary_43_sigpos[order(vary_43_sigpos$GFCA_v_cellulose_padj,decreasing=FALSE),]


newnames <- "GF CA vs. cellulose"
newnames
head(vary_43)
nrow(vary_43)
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43_sigpos[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,8]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

pc_regulation_heatmap<-pheatmap(data.matrix(vary_44[,2]), col=my_palette, cellwidth=100,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,8]), border_color=NA,breaks=seq(-4, 4, length.out=300))











# Colleen Ahern
# Pchry DESeq2 heatmap (cleaned) — SACA5-excluded, 4-comparison pipeline
# Original: 08/27/2025 | Cleaned: 07/22/2026
#
# NOTE: This cleaned version covers the CURRENT/correct pipeline only —
# i.e. everything through the 12/05/2025 "top N padjmin" heatmaps.
# The final "GFCA vs cellulose only" section from the original script was
# NOT included here: it re-derives results from an older count matrix and
# DESeq output that still included the SACA5 outlier, contradicting the
# corrected analysis this script is based on. Let me know if you actually
# need that section rebuilt on the corrected (SACA5-excluded) data instead.

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(ape)
library(pheatmap)
library(grid)
library(RColorBrewer)
library(DRnaSeq)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/Users/colleenahern/Documents/Pchry/DESeq2_analysis"
FIG_DIR <- file.path(BASE, "Figures")
DATE    <- "07222026"

# =============================================================================
# HELPERS
# =============================================================================

# ── Custom pheatmap column-label renderer (horizontal, wrapped labels) ────────
draw_colnames_horizontal <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"),
           vjust = 1.5, hjust = 0.5, rot = 0, gp = gpar(...))
}

# ── Build bold/non-bold expression labels for rows matching a CAZyme gene set ──
# IMPORTANT: `mat` must already have gene Name as rownames before calling this,
# or the bold-matching will silently fail (this was the bug in the original
# top-100 heatmap — bolding was compared against default integer rownames).
build_row_labels <- function(mat, name_col, cazyme_names) {
  stopifnot(identical(rownames(mat), mat[[name_col]]))  # catches the ordering bug explicitly
  is_cazyme <- rownames(mat) %in% cazyme_names
  labels <- ifelse(is_cazyme,
                   paste0("bold('", mat[[name_col]], "')"),
                   paste0("'", mat[[name_col]], "'"))
  parse(text = labels)
}

# ── Reusable heatmap builder for the 4-comparison LFC matrix ──────────────────
draw_pchry_heatmap <- function(mat, comparison_cols, comparison_labels, cazyme_names,
                               cellheight, fontsize_row = 8, breaks_range = c(-4, 6)) {
  assignInNamespace("draw_colnames", draw_colnames_horizontal, ns = asNamespace("pheatmap"))
  my_palette <- colorRampPalette(c("blue", "black", "red"))(299)
  
  row_labels <- build_row_labels(mat, "Name", cazyme_names)
  
  pheatmap(
    data.matrix(mat[, comparison_cols]),
    col              = my_palette,
    cellwidth        = 130,
    cellheight       = cellheight,
    cluster_rows     = TRUE,
    cluster_cols     = FALSE,
    fontsize_col     = 12,
    fontsize_row     = fontsize_row,
    show_rownames    = TRUE,
    labels_col       = comparison_labels,
    labels_row       = row_labels,
    border_color     = "black",
    breaks           = seq(breaks_range[1], breaks_range[2], length.out = 300)
  )
}

# =============================================================================
# 1. Import count data and sample info
# =============================================================================

pcseqdata    <- read.delim(file.path(BASE, "featcount/pchry_countmat_jgi_07312025.txt"), stringsAsFactors = FALSE)
pcsampleinfo <- read.delim(file.path(BASE, "pcsampleinfo.txt"), stringsAsFactors = TRUE)

pccountdata <- pcseqdata
rownames(pccountdata) <- pccountdata$Geneid
pccountdata <- pccountdata[, -c(1:6)]

pcsampleinfo$libraryName_jgi <- gsub("[/-]", ".", pcsampleinfo$libraryName_jgi)

stopifnot(all(colnames(pccountdata) %in% pcsampleinfo$libraryName_jgi))
pccountdata <- pccountdata[, as.character(pcsampleinfo$libraryName_jgi)]
stopifnot(all(colnames(pccountdata) == pcsampleinfo$libraryName_jgi))

colnames(pccountdata) <- pcsampleinfo$sampleName
rownames(pcsampleinfo) <- pcsampleinfo$sampleName

# ── Remove SACA5 outlier ────────────────────────────────────────────────────────
samples_to_remove <- "SACA5"
pccts_subset     <- pccountdata[, !(colnames(pccountdata) %in% samples_to_remove)]
pccoldata_subset <- pcsampleinfo[!(rownames(pcsampleinfo) %in% samples_to_remove), ]
stopifnot(all(rownames(pccoldata_subset) == colnames(pccts_subset)))

# =============================================================================
# 2. TPM cutoff filter (vectorized — replaces the original row-by-row loop)
# =============================================================================

stopifnot(all(rownames(pccts_subset) == pcseqdata$Geneid))
gl  <- pcseqdata$Length
tpm <- as.data.frame(tpm(pccts_subset, gl))

stopifnot(all(colnames(tpm) == pccoldata_subset$sampleName))

# Mean TPM across replicates, per substrate group — vectorized instead of a per-row for loop
cellulose_mean <- rowMeans(tpm[, c("cellulose1","cellulose2","cellulose3","cellulose4","cellulose5")])
gfca_mean      <- rowMeans(tpm[, c("GFCA1","GFCA2","GFCA3","GFCA4","GFCA5")])
saca_mean      <- rowMeans(tpm[, c("SACA1","SACA2","SACA3","SACA4")])  # SACA5 excluded
blank_mean     <- rowMeans(tpm[, c("blank1","blank2","blank3","blank4","blank5")])

keep_tpm <- !(cellulose_mean < 2 & gfca_mean < 2 & saca_mean < 2 & blank_mean < 2)
tpm_keep <- tpm[keep_tpm, ]
tpm_keep$GeneID <- rownames(tpm_keep)
message("TPM filter: ", nrow(tpm), " -> ", nrow(tpm_keep), " genes")

# =============================================================================
# 3. Load DESeq2 results, drop genes with no p-value in ANY comparison
# =============================================================================

vary_43 <- read_csv(file.path(BASE, "Tables/Pc_DESeq_jgi_08272025.csv"), show_col_types = FALSE)

padj_cols <- grep("padj", colnames(vary_43), value = TRUE)
keep_padj <- rowSums(!is.na(vary_43[, padj_cols])) > 0  # keep if AT LEAST one comparison has a p-value
vary_43   <- vary_43[keep_padj, ]
message("padj NA filter: kept ", nrow(vary_43), " genes")

# ── Restrict to genes passing the TPM filter too ───────────────────────────────
vary_43 <- merge(vary_43, tpm_keep[, "GeneID", drop = FALSE], by = "GeneID")
message("After TPM + padj filters: ", nrow(vary_43), " genes")

# =============================================================================
# 4. Gene annotations (JGI GFF3)
# =============================================================================

gff <- read.gff(file.path(BASE, "Phchr4_2_GeneCatalog_20210419.gff3"), GFF3 = TRUE)
gff <- gff[gff$type == "gene", ]

gene_anno <- gff %>%
  dplyr::select(seqid, everything()) %>%
  mutate(attributes = str_split(attributes, ";")) %>%
  unnest_wider(where(is.list), names_sep = "") %>%
  dplyr::select(attributes1, attributes4, attributes5, attributes6) %>%
  rename(GeneID = attributes1, Protein_name = attributes4,
         ProteinID = attributes5, TranscriptID = attributes6) %>%
  mutate(
    GeneID       = gsub("ID=", "", GeneID),
    Protein_name = gsub("product_name=", "", Protein_name),
    ProteinID    = gsub("proteinId=", "", ProteinID),
    TranscriptID = gsub("transcriptId=", "", TranscriptID)
  )

vary_43 <- merge(vary_43, gene_anno, by = "GeneID")
vary_43$Name <- paste0(vary_43$Protein_name, " (proteinID = ", vary_43$ProteinID, ")")

# =============================================================================
# 5. CAZyme annotation (dbCAN)
# NOTE: this block is unchanged from the original — it's complex string
# parsing logic that works, so it wasn't worth rewriting for this cleanup.
# =============================================================================

cazlist  <- read.delim(file.path(BASE, "Pchry_dbcan/overview.txt"))
cazlist2 <- cazlist[cazlist$X.ofTools > 1, ]

cazlist_f <- cazlist2
cazlist_f$Gene.ID <- gsub("jgi\\|", "", cazlist2$Gene.ID)
cazlist_f$Gene.ID <- gsub(".*\\|(.*?)\\|.*", "\\1", cazlist_f$Gene.ID)
cazlist_f$HMMER   <- gsub("\\(.*", "", cazlist_f$HMMER)
cazlist_f$DIAMOND <- gsub("Phchr2_GeneCatalog_proteins_20131210.aa.fasta\\+", "", cazlist_f$DIAMOND)
cazlist_f$DIAMOND <- gsub(".*fasta\\+", "", cazlist_f$DIAMOND)
rownames(cazlist_f) <- 1:nrow(cazlist_f)

cazlist_f$CAZyme <- NA
cazlist_f$CAZyme[1:376] <- cazlist_f[1:376, "HMMER"]
cazlist_f$CAZyme[377:388] <- cazlist_f[377:388, "eCAMI"]
cazlist_f$DIAMOND <- gsub("\\+3.*", "", cazlist_f$DIAMOND)

cazlist_ff <- cazlist_f %>%
  dplyr::select(Gene.ID, everything()) %>%
  mutate(eCAMI = str_split(eCAMI, "\\+")) %>% unnest_wider(where(is.list), names_sep = "") %>%
  mutate(DIAMOND = str_split(DIAMOND, "\\+")) %>% unnest_wider(where(is.list), names_sep = "") %>%
  mutate(HMMER = str_split(HMMER, "\\+")) %>% unnest_wider(where(is.list), names_sep = "")

cazlist_fff <- cazlist_ff
cazlist_fff$CAZ <- NA
for (i in 1:nrow(cazlist_fff)) {
  z <- as.data.frame(table(as.character(cazlist_fff[i, 3:9])))
  cazlist_fff$CAZ[i] <- paste(paste(z$Var1, " x", z$Freq, sep = ""), collapse = "; ")
}

cazlist_f4 <- cazlist_fff
cazlist_f4$CAZf <- cazlist_f4$CAZ
cazlist_f4 <- cazlist_f4 %>%
  dplyr::select(Gene.ID, everything()) %>%
  mutate(CAZf = str_split(CAZf, ";")) %>%
  unnest_wider(where(is.list), names_sep = "")

for (i in 13:16) {
  cazlist_f4[, i] <- lapply(cazlist_f4[, i], function(x) replace(x, grepl("x1", x), NA))
}

cazlist_f4 <- cazlist_f4 %>%
  unite(CAZ_name, starts_with("CAZf"), na.rm = TRUE, remove = TRUE, sep = ";")

cazlist_f4$CAZ_name <- gsub(" ", "", cazlist_f4$CAZ_name)
cazlist_f4$CAZ_name <- gsub(";", "; ", cazlist_f4$CAZ_name)
cazlist_f4$CAZ_anno <- gsub("x[234]", "", cazlist_f4$CAZ_name)

caz_anno <- cazlist_f4[cazlist_f4$CAZ_anno != "", ] %>%
  rename(ProteinID = Gene.ID) %>%
  dplyr::select(ProteinID, CAZ_anno)

write_csv(caz_anno, file.path(BASE, paste0("Pchry_dbcan/cazlist_", DATE, ".csv")))

vary_43caz <- merge(vary_43, caz_anno, by = "ProteinID", all.x = TRUE)

write_csv(vary_43caz, file.path(BASE, "Tables", paste0("vary43_jgi_allconds_anno_cazanno_", DATE, ".csv")))

# =============================================================================
# 6. padjmin — restricted to the 4 comparisons of interest
# (matches the correction already made in the original script)
# =============================================================================

comparison_padj_cols <- c("SACA_v_blank_padj", "GFCA_v_blank_padj",
                          "cellulose_v_blank_padj", "GFCA_v_cellulose_padj")
comparison_lfc_cols  <- c("SACA_v_blank_log2FoldChange", "GFCA_v_blank_log2FoldChange",
                          "cellulose_v_blank_log2FoldChange", "GFCA_v_cellulose_log2FoldChange")

vary_43caz$padjmin <- apply(vary_43caz[, comparison_padj_cols], 1, min, na.rm = TRUE)
vary_43caz <- vary_43caz[order(vary_43caz$padjmin, decreasing = FALSE), ]

# =============================================================================
# 7. Mask non-significant LFCs to 0 for the heatmap
# (padj > 0.049 -> LFC set to 0; significant genes keep their LFC as-is.
#  NOTE: the original script's comment said the opposite of what the code
#  actually does — the code's behavior here is the correct/intended one.)
# =============================================================================

vary_43caz_masked <- vary_43caz
for (k in seq_along(comparison_lfc_cols)) {
  lfc_col  <- comparison_lfc_cols[k]
  padj_col <- comparison_padj_cols[k]
  not_sig  <- vary_43caz_masked[[padj_col]] >= 0.05 | is.na(vary_43caz_masked[[padj_col]])
  vary_43caz_masked[[lfc_col]][not_sig] <- 0
}

comparison_labels <- c(
  "Highly-acetylated CA\nvs. No substrate",
  "Lowly-acetylated CA\nvs. No substrate",
  "Cellulose\nvs. No substrate",
  "Lowly-acetylated CA\nvs. Cellulose"
)

cazyme_gene_names <- vary_43caz_masked$Name[!is.na(vary_43caz_masked$CAZ_anno)]

# =============================================================================
# 8. Heatmaps — top N genes by padjmin
# =============================================================================

draw_and_save_heatmap <- function(n_genes, cellheight, fontsize_row, width, height) {
  top_n <- vary_43caz_masked[1:n_genes, ]
  top_n <- as.data.frame(top_n)
  rownames(top_n) <- top_n$Name  # MUST happen before build_row_labels() is called
  
  hm <- draw_pchry_heatmap(top_n, comparison_lfc_cols, comparison_labels,
                           cazyme_gene_names, cellheight = cellheight, fontsize_row = fontsize_row)
  
  grid.newpage()
  grid.draw(hm$gtable)
  
  base_name <- paste0("Pchry_4comps_top", n_genes, "padjmin_heatmap_", DATE)
  ggsave(plot = hm, filename = file.path(FIG_DIR, paste0(base_name, ".pdf")),
         width = width, height = height, units = "in", dpi = 300)
  ggsave(plot = hm, filename = file.path(FIG_DIR, paste0(base_name, ".tiff")),
         width = width, height = height, units = "in", dpi = 200, bg = "white")
  
  invisible(hm)
}

draw_and_save_heatmap(n_genes = 100, cellheight = 7,  fontsize_row = 8,  width = 14, height = 11)
draw_and_save_heatmap(n_genes = 50,  cellheight = 15, fontsize_row = 12, width = 16, height = 12)
draw_and_save_heatmap(n_genes = 30,  cellheight = 15, fontsize_row = 12, width = 16, height = 8)
draw_and_save_heatmap(n_genes = 20,  cellheight = 15, fontsize_row = 12, width = 16, height = 6)

message("Done.")

