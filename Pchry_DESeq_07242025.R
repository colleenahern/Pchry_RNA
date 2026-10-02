# Colleen Ahern
# 07/24/2025
# Complete DESeq analysis on Pchry count matrix from featureCounts
# following tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
library("tximport")
library("readr")
BiocManager::install("tximportData")
library("tximportData")
library(tidyverse)
BiocManager::install("pasilla")
library(pasilla)

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
table(colnames(pccountdata)==pcsampleinfo$libraryName_jgi)
# pccountdata <- pccountdata[,as.character(pcsampleinfo$libraryName)]
pccountdata <- pccountdata[,as.character(pcsampleinfo$libraryName_jgi)]
# table(colnames(pccountdata)==pcsampleinfo$libraryName) # Now they match
table(colnames(pccountdata)==pcsampleinfo$libraryName_jgi) # Now they match
colnames(pccountdata) <- pcsampleinfo$sampleName

rownames(pcsampleinfo) <- pcsampleinfo$sampleName
pccoldata <- pcsampleinfo
pccts <- pccountdata

all(rownames(pccoldata) == colnames(pccts))

# With the count matrix and the sample info, construct a DESeqDataSet:
library("DESeq2")
pcdds <- DESeqDataSetFromMatrix(countData = pccts,
                                colData = pccoldata,
                                design = ~ groupName)
pcdds

# Pre-filter
keep <- rowSums(counts(pcdds)) > 0
pcdds <- pcdds[keep,]   # reduced number of genes from 12256 to 11718 (11815 for JGI)

# Set the factor levels 
pcdds$groupName <- factor(pcdds$groupName, levels = c("blank","cellulose","GFCA", "SACA"))

# Differential expression analysis: compare each substrate to glucose
pcdds <- DESeq(pcdds)
pcres <- results(pcdds)
pcres
pcres <- results(pcdds, contrast=c("groupName","SACA","blank"))
pcres <- results(pcdds, contrast=c("groupName","blank","SACA"))
pcres

# Now compare all of the different the pc substrates to each other: 6 in total
# SACA vs blank
pcres_SACA_blank <- results(pcdds, contrast=c("groupName","SACA","blank"))
pcres_SACA_blank
pcresmat_SACA_blank <- as.matrix(pcres_SACA_blank)
head(pcresmat_SACA_blank)
colnames(pcresmat_SACA_blank) <- paste("SACA_v_blank", colnames(pcresmat_SACA_blank), sep = "_")
pcresmat_SACA_blank <- as.data.frame(pcresmat_SACA_blank)
head(pcresmat_SACA_blank)
sum(pcresmat_SACA_blank$SACA_v_blank_padj < 0.05, na.rm=TRUE)

# GFCA vs blank
pcres_GFCA_blank <- results(pcdds, contrast=c("groupName","GFCA","blank"))
pcres_GFCA_blank
pcresmat_GFCA_blank <- as.matrix(pcres_GFCA_blank)
pcresmat_GFCA_blank
colnames(pcresmat_GFCA_blank) <- paste("GFCA_v_blank", colnames(pcresmat_GFCA_blank), sep = "_")
pcresmat_GFCA_blank <- as.data.frame(pcresmat_GFCA_blank)
head(pcresmat_GFCA_blank)
sum(pcresmat_GFCA_blank$GFCA_v_blank_padj < 0.05, na.rm=TRUE)

# cellulose vs blank
pcres_cellulose_blank <- results(pcdds, contrast=c("groupName","cellulose","blank"))
pcres_cellulose_blank
pcresmat_cellulose_blank <- as.matrix(pcres_cellulose_blank)
head(pcresmat_cellulose_blank)
colnames(pcresmat_cellulose_blank) <- paste("cellulose_v_blank", colnames(pcresmat_cellulose_blank), sep = "_")
pcresmat_cellulose_blank <- as.data.frame(pcresmat_cellulose_blank)
head(pcresmat_cellulose_blank)
sum(pcresmat_cellulose_blank$cellulose_v_blank_padj < 0.05, na.rm=TRUE)

# cellulose vs SACA
pcres_cellulose_SACA <- results(pcdds, contrast=c("groupName","cellulose","SACA"))
pcres_cellulose_SACA
pcresmat_cellulose_SACA <- as.matrix(pcres_cellulose_SACA)
head(pcresmat_cellulose_SACA)
colnames(pcresmat_cellulose_SACA) <- paste("cellulose_v_SACA", colnames(pcresmat_cellulose_SACA), sep = "_")
pcresmat_cellulose_SACA <- as.data.frame(pcresmat_cellulose_SACA)
head(pcresmat_cellulose_SACA)
sum(pcresmat_cellulose_SACA$cellulose_v_SACA_padj < 0.05, na.rm=TRUE)

# GFCA vs SACA
pcres_GFCA_SACA <- results(pcdds, contrast=c("groupName","GFCA","SACA"))
pcres_GFCA_SACA
pcresmat_GFCA_SACA <- as.matrix(pcres_GFCA_SACA)
head(pcresmat_GFCA_SACA)
colnames(pcresmat_GFCA_SACA) <- paste("GFCA_v_SACA", colnames(pcresmat_GFCA_SACA), sep = "_")
pcresmat_GFCA_SACA <- as.data.frame(pcresmat_GFCA_SACA)
head(pcresmat_GFCA_SACA)
sum(pcresmat_GFCA_SACA$GFCA_v_SACA_padj < 0.05, na.rm=TRUE)

# GFCA vs cellulose
pcres_GFCA_cellulose <- results(pcdds, contrast=c("groupName","GFCA","cellulose"))
pcres_GFCA_cellulose
pcresmat_GFCA_cellulose <- as.matrix(pcres_GFCA_cellulose)
head(pcresmat_GFCA_cellulose)
colnames(pcresmat_GFCA_cellulose) <- paste("GFCA_v_cellulose", colnames(pcresmat_GFCA_cellulose), sep = "_")
pcresmat_GFCA_cellulose <- as.data.frame(pcresmat_GFCA_cellulose)
head(pcresmat_GFCA_cellulose)
sum(pcresmat_GFCA_cellulose$GFCA_v_cellulose_padj < 0.05, na.rm=TRUE)

# Confirm that the order of genes in all 6 matrices is the same
all(rownames(pcresmat_SACA_blank) == rownames(pcresmat_GFCA_blank))
all(rownames(pcresmat_GFCA_blank) == rownames(pcresmat_cellulose_blank))
all(rownames(pcresmat_cellulose_blank) == rownames(pcresmat_cellulose_SACA))
all(rownames(pcresmat_cellulose_SACA) == rownames(pcresmat_GFCA_SACA))
all(rownames(pcresmat_GFCA_SACA) == rownames(pcresmat_GFCA_cellulose))

# original Pc_DESeq output from 07/24/2025
Pc_DESeq <- cbind(pcresmat_SACA_blank, pcresmat_GFCA_blank)
Pc_DESeq <- Pc_DESeq[, c(2,6,8,12)]
Pc_DESeq <- as.data.frame(Pc_DESeq)
Pc_DESeq$SACA_v_blank_Sig <- ifelse(Pc_DESeq$SACA_v_blank_padj <= 0.05,"TRUE","FALSE")
Pc_DESeq <- Pc_DESeq[, c(1,2,5,3,4)]
Pc_DESeq$GFCA_v_blank_Sig <- ifelse(Pc_DESeq$GFCA_v_blank_padj <= 0.05,"TRUE","FALSE")

Pc_DESeq <- cbind(Pc_DESeq, pcresmat_cellulose_blank)
Pc_DESeq <- Pc_DESeq[, c(1:6,8,12)]
Pc_DESeq$cellulose_v_blank_Sig <- ifelse(Pc_DESeq$cellulose_v_blank_padj <= 0.05,"TRUE","FALSE")

Pc_DESeq <- cbind(Pc_DESeq, pcresmat_cellulose_SACA)
Pc_DESeq <- Pc_DESeq[, c(1:9,11,15)]
Pc_DESeq$cellulose_v_SACA_Sig <- ifelse(Pc_DESeq$cellulose_v_SACA_padj <= 0.05,"TRUE","FALSE")

Pc_DESeq <- cbind(Pc_DESeq, pcresmat_GFCA_SACA)
Pc_DESeq <- Pc_DESeq[, c(1:12,14,18)]
Pc_DESeq$GFCA_v_SACA_Sig <- ifelse(Pc_DESeq$GFCA_v_SACA_padj <= 0.05,"TRUE","FALSE")

Pc_DESeq <- cbind(Pc_DESeq, pcresmat_GFCA_cellulose)
Pc_DESeq <- Pc_DESeq[, c(1:15,17,21)]
Pc_DESeq$GFCA_v_cellulose_Sig <- ifelse(Pc_DESeq$GFCA_v_cellulose_padj <= 0.05,"TRUE","FALSE")

Pc_DESeq$GeneID <- rownames(Pc_DESeq)
Pc_DESeq <- Pc_DESeq[,c(19,1:18)]

library(writexl)
#write_xlsx(Pc_DESeq, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_07242025.xlsx")
write_xlsx(Pc_DESeq, "/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_jgi_08012025.xlsx")


# 08/18/2025 addition - For GSEA, I want the p values, not just the padjusted. So let's make a Pc_DESeq_all file that has all outputs from DESeq
# so that this isn't a problem in the future
Pc_DESeq_all <- cbind(pcresmat_SACA_blank, pcresmat_GFCA_blank)
Pc_DESeq_all <- cbind(Pc_DESeq_all, pcresmat_cellulose_blank)
Pc_DESeq_all <- cbind(Pc_DESeq_all, pcresmat_cellulose_SACA)
Pc_DESeq_all <- cbind(Pc_DESeq_all, pcresmat_GFCA_SACA)
Pc_DESeq_all <- cbind(Pc_DESeq_all, pcresmat_GFCA_cellulose)

Pc_DESeq_all$GeneID <- rownames(Pc_DESeq_all)
Pc_DESeq_all <- Pc_DESeq_all[,c(37,1:36)]
write_csv(Pc_DESeq_all,"/Users/colleenahern/Documents/Pchry/DESeq2_analysis/Pc_DESeq_jgi_ALL_08182025.csv")
