# Colleen Ahern
# Final manuscript heatmap from my Pchry DESeq2 analysis 

# all genes with CAZymes labeled with dbCAN annotation; tpm filter applied and padjmin ranked
# Also includes supplement heatmaps: all genes and CAZyme-only

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
# pcseqdata <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/featcount/pchry_countmat.txt", stringsAsFactors = FALSE)
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

#pcsampleinfo$libraryName <- gsub("/",".",pcsampleinfo$libraryName)
pcsampleinfo$libraryName_jgi <- gsub("/",".",pcsampleinfo$libraryName_jgi)
#pcsampleinfo$libraryName <- gsub("-",".",pcsampleinfo$libraryName)
pcsampleinfo$libraryName_jgi <- gsub("-",".",pcsampleinfo$libraryName_jgi)

# Make sure that the column names are the same as libraryName in the pcsampleinfo file
colnames(pccountdata) # using substr, you can change the characters of the colnames
# colnames(pccountdata) %in% pcsampleinfo$libraryName
colnames(pccountdata) %in% pcsampleinfo$libraryName_jgi
table(colnames(pccountdata)==pcsampleinfo$libraryName)
# pccountdata <- pccountdata[,as.character(pcsampleinfo$libraryName)]
pccountdata <- pccountdata[,as.character(pcsampleinfo$libraryName_jgi)]
# table(colnames(pccountdata)==pcsampleinfo$libraryName) # Now they match
table(colnames(pccountdata)==pcsampleinfo$libraryName_jgi) # Now they match
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

tpm <- tpm[(tpm$rem=="KEEP"),] # reduces from 12256 genes to 10843 (10894 for jgi)
tpm <- tpm[,-c(21)]
tpm$GeneID <- rownames(tpm)

# Import list of differentially regulated genes for heatmap
# Remove genes with a p-value = NA for all three comparitive conditions

#vary_43 <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_07242025.txt", header=T, check.names = FALSE) # My DESeq
vary_43 <- read.delim("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_jgi_08012025.txt", header=T, check.names = FALSE) # My DESeq
# vary_43 <- vary_43[,c(1:10)] 
# for jgi -- keep all comparisons

vary_43$rem <- 0
# for (i in 1:nrow(vary_43)) {
#   vary_43$rem[i] <- ifelse(is.na(vary_43$SACA_v_blank_padj[i]) == TRUE && is.na(vary_43$GFCA_v_blank_padj[i]) == TRUE && is.na(vary_43$cellulose_v_blank_padj[i]) == TRUE, "REMOVE", "KEEP")
# }
for (i in 1:nrow(vary_43)) {
  vary_43$rem[i] <- ifelse(is.na(vary_43$SACA_v_blank_padj[i]) == TRUE && is.na(vary_43$GFCA_v_blank_padj[i]) == TRUE && is.na(vary_43$cellulose_v_blank_padj[i]) == TRUE && is.na(vary_43$cellulose_v_SACA_padj[i]) == TRUE && is.na(vary_43$GFCA_v_SACA_padj[i]) == TRUE && is.na(vary_43$GFCA_v_cellulose_padj[i]) == TRUE, "REMOVE", "KEEP")
}

vary_43 <- vary_43[(vary_43$rem=="KEEP"),] # reduces from 11718 to 11690 (for JGI, from 11815 to 11791)
# vary_43 <- vary_43[,-11]
vary_43 <- vary_43[,-20]


# Merge tpm matrix with vary_43 matrix - to make a matrix that contains only genes that pass both of those filters
a <- merge(vary_43,tpm,"GeneID") # 10,822 genes (10874 for jgi)
# vary_43 <- a[,c(1:10)]
vary_43 <- a[,c(1:19)]

# Create a minimum p value column that takes the minimum p value from any of the three comparative conditions for each gene
# Then sort the matrix from lowest padjmin to largest
# vary_43$padjmin <- apply(vary_43[,c("SACA_v_blank_padj","GFCA_v_blank_padj","cellulose_v_blank_padj")],1,min, na.rm = TRUE)
vary_43$padjmin <- apply(vary_43[,c(colnames(vary_43)[grep("padj",colnames(vary_43))])],1,min, na.rm = TRUE)
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]
# vary_43 <- vary_43[,-c(4,7,10)]
# colnames(vary_43) <- c("GeneID", "SACA_v_Blank", "SACA_v_blank_padj", "GFCA_v_Blank", "GFCA_v_blank_padj", "Cellulose_v_Blank", "cellulose_v_blank_padj", "padjmin")

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

# for jgi - I want to have the DESeq results (tpm filtered) before and after merging with the annotated genes cus I don't want to
# lose non-annotated genes
write_csv(vary_43, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_beforeanno_08042025.csv")
# Add annotations to the data frame and reorder padjmin
vary_43 <- merge(vary_43, my_gene_col_44, "GeneID")
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]
# Uhhh actually I think i didn't lose any genes by merging, i guess the genome is really well annotated woohoo

for (i in 1:nrow(vary_43)) {
  vary_43$Name[i] <- paste(vary_43$Protein_name[i], " ", "(proteinID = ", vary_43$ProteinID[i], ")", sep = "")
}

# write_csv(vary_43, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary_43_07242025.csv")
# vary_43 <- read_csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary_43_07242025.csv")
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]
write_csv(vary_43, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_08042025.csv")
vary_43 <- read.csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_08042025.csv")

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
write_csv(vary_43caz, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_cazanno_08072025.csv")
vary_43caz <- read_csv("/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Tables/vary43_jgi_allconds_anno_cazanno_08072025.csv")
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

newnames <- lapply(
  colnames(vary_43[,c(2,4,6)]),
  function(x) bquote(.(x)))
newnames
head(vary_43)
nrow(vary_43)
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,12]
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

pc_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(2,4,6)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,12]), border_color=NA,breaks=seq(-4, 6, length.out=300))




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














