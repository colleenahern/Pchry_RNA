# =============================================================================
# GO Biological Process enrichment: GF CA vs. Cellulose in P. chrysosporium
# Combined figure: (A) curated GFCA-vs-Cellulose GO dotplot
#                  (B) carbohydrate metabolism GO term overlap (GFCA/Cellulose vs. blank)
# C. Ahern
# =============================================================================

library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)
library(cowplot)
library(ape)

set.seed(42)

BASE    <- "/Users/colleenahern/Documents/Pchry/DESeq2_analysis"
DATE    <- format(Sys.Date(), "%m%d%Y")
FIG_DIR <- "/Users/colleenahern/Documents/Pchry/Figures"

# =============================================================================
# 1. Read in eggNOG-mapper annotations and map ProteinID -> GeneID
# =============================================================================
eggnog <- read_tsv(
  file.path(BASE, "Pchry_JGI_files/JGI_annotations/Phchr_eggnog.emapper.annotations"),
  comment = "##"
)

colnames(eggnog)[1] <- "ProteinID_raw"

eggnog <- eggnog %>%
  mutate(ProteinID = str_extract(ProteinID_raw, "(?<=jgi\\|Phchr4_2\\|)\\d+"))

gene_protein_map_raw <- read.gff(
  file.path(BASE, "Phchr4_2_GeneCatalog_20210419.gff3"),
  GFF3 = TRUE
) %>%
  filter(type == "gene")

gene_protein_map <- gene_protein_map_raw %>%
  dplyr::select(seqid, everything()) %>%
  mutate(attributes = str_split(attributes, ";")) %>%
  unnest_wider(where(is.list), names_sep = "") %>%
  dplyr::select(attributes1, attributes4, attributes5, attributes6) %>%
  mutate(across(everything(), ~ gsub("^(ID=|product_name=|proteinId=|transcriptId=)", "", .))) %>%
  set_names(c("GeneID", "Protein_name", "ProteinID", "TranscriptID"))

id_map <- gene_protein_map %>%
  dplyr::select(GeneID, ProteinID) %>%
  distinct()

eggnog_geneid <- eggnog %>%
  left_join(id_map, by = "ProteinID") %>%
  filter(!is.na(GeneID))

# =============================================================================
# 2. Build GO Biological Process TERM2GENE / TERM2NAME from eggNOG
# =============================================================================
go_long <- eggnog_geneid %>%
  dplyr::filter(GOs != "-", !is.na(GOs)) %>%
  dplyr::select(GeneID, GOs) %>%
  separate_rows(GOs, sep = ",") %>%
  dplyr::rename(go_id = GOs) %>%
  dplyr::distinct()

go_bp_terms <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_long$go_id),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
) %>%
  filter(ONTOLOGY == "BP")

go_bp_long <- go_long %>%
  filter(go_id %in% go_bp_terms$GOID)

term2gene <- go_bp_long %>%
  dplyr::rename(term = go_id, gene = GeneID) %>%
  dplyr::select(term, gene)

term2name <- go_bp_terms %>%
  dplyr::rename(term = GOID, name = TERM) %>%
  dplyr::select(term, name)

# =============================================================================
# 3. Read DESeq2 + CAZyme annotated results table
# =============================================================================
vary_43caz <- read_csv(file.path(BASE, "Tables/vary43_jgi_allconds_anno_cazanno_07222026.csv"))

# =============================================================================
# 4. Helper: build a ranked vector and run GSEA for a given stat column
# =============================================================================
run_go_gsea <- function(data, stat_col) {
  ranks <- data %>%
    filter(!is.na(.data[[stat_col]])) %>%
    dplyr::select(GeneID, all_of(stat_col)) %>%
    distinct(GeneID, .keep_all = TRUE) %>%
    arrange(desc(.data[[stat_col]]))
  
  vec <- ranks[[stat_col]]
  names(vec) <- as.character(ranks$GeneID)
  vec <- sort(vec, decreasing = TRUE)
  
  GSEA(
    geneList      = vec,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    minGSSize     = 5,
    maxGSSize     = 300,
    eps           = 0,
    by            = "fgsea"
  )
}

gsea_go            <- run_go_gsea(vary_43caz, "GFCA_v_cellulose_stat")
gsea_go_gfca_blank <- run_go_gsea(vary_43caz, "GFCA_v_blank_stat")
gsea_go_cell_blank <- run_go_gsea(vary_43caz, "cellulose_v_blank_stat")

go_res_gfca_blank <- as.data.frame(gsea_go_gfca_blank)
go_res_cell_blank <- as.data.frame(gsea_go_cell_blank)

# =============================================================================
# 5. Panel A -- curated GFCA vs. Cellulose GO dotplot
# =============================================================================
cellulose_terms <- c(
  "iron ion transmembrane transport",
  "iron ion transport",
  "intracellular iron ion homeostasis",
  "cellular response to iron ion starvation"
)

gfca_terms <- c(
  "glucan metabolic process",
  "carbohydrate biosynthetic process",
  "polysaccharide biosynthetic process",
  "glucose metabolic process",
  "glucan biosynthetic process",
  "carbohydrate metabolic process",
  "glucose catabolic process",
  "glycolytic process through fructose-6-phosphate",
  "glycolytic process through glucose-6-phosphate",
  "canonical glycolysis",
  "glucose catabolic process to pyruvate",
  "fungal-type cell wall beta-glucan biosynthetic process",
  "cell wall beta-glucan metabolic process",
  "fungal-type cell wall beta-glucan metabolic process",
  "glycolytic process"
)

go_narrative_curated <- gsea_go@result %>%
  filter(Description %in% c(cellulose_terms, gfca_terms)) %>%
  mutate(
    panel = ifelse(NES > 0, "GF CA-enriched", "Cellulose-enriched"),
    panel = factor(panel, levels = c("GF CA-enriched", "Cellulose-enriched")),
    Description = str_to_sentence(Description)
  ) %>%
  mutate(
    Description = case_when(
      panel == "Cellulose-enriched" ~ str_wrap(Description, width = 25),
      TRUE ~ str_wrap(Description, width = 35)
    ),
    Description = fct_reorder(Description, NES)
  )

nrow(go_narrative_curated)
table(go_narrative_curated$panel)

go_curated_plot <- ggplot(go_narrative_curated, aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
  geom_point() +
  facet_wrap(~ panel, scales = "free_y") +
  scale_color_gradient(low = "#377EB8", high = "#E41A1C") +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = "NES", y = NULL, size = "Count", color = "Adjusted p-value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, hjust = 0.5, color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    panel.spacing = unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

go_curated_plot

ggsave(plot = go_curated_plot,
       filename = file.path(FIG_DIR, paste0("Pchry_GOplot", DATE, ".pdf")),
       width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = go_curated_plot,
       filename = file.path(FIG_DIR, paste0("Pchry_GOplot", DATE, ".tiff")),
       width = 10, height = 10, units = "in", dpi = 200, bg = "white")


# =============================================================================
# 6. Panel B -- carbohydrate metabolism GO term overlap
# =============================================================================
carb_descendants <- as.list(GOBPOFFSPRING)[["GO:0005975"]]
carb_go_ids <- c("GO:0005975", carb_descendants)

gfca_carb_precise <- go_res_gfca_blank %>% filter(NES > 0, ID %in% carb_go_ids)
cell_carb_precise <- go_res_cell_blank %>% filter(NES > 0, ID %in% carb_go_ids)

shared_carb     <- intersect(gfca_carb_precise$Description, cell_carb_precise$Description)
gfca_only_carb  <- setdiff(gfca_carb_precise$Description, cell_carb_precise$Description)
cell_only_carb  <- setdiff(cell_carb_precise$Description, gfca_carb_precise$Description)

carb_overlap_df <- tibble(
  condition = rep(c("GF CA vs.\nNo substrate", "Cellulose vs.\nNo substrate"), each = 2),
  category  = rep(c("Shared between GF CA and Cellulose", "Unique to GF CA"), 2),
  count     = c(length(shared_carb), length(gfca_only_carb),
                length(shared_carb), length(cell_only_carb))
) %>%
  mutate(
    condition = factor(condition, levels = c("GF CA vs.\nNo substrate", "Cellulose vs.\nNo substrate")),
    category  = factor(category, levels = c("Unique to GF CA", "Shared between GF CA and Cellulose")),
    label     = ifelse(count == 0, "", as.character(count))
  )

carb_overlap_plot <- ggplot(carb_overlap_df, aes(x = condition, y = count, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4.2, color = "white") +
  scale_fill_manual(values = c("Unique to GF CA" = "#1D9E75", "Shared between GF CA and Cellulose" = "#534AB7")) +
  labs(x = NULL, y = "Number of significantly enriched\ncarbohydrate metabolism GO terms", fill = NULL) +
  theme_bw() +
  theme(
    axis.text        = element_text(size = 12, color = "black"),
    axis.title.y     = element_text(size = 12, color = "black"),
    axis.title.x     = element_text(size = 12, color = "black"),
    axis.ticks       = element_line(color = "black"),
    legend.text      = element_text(size = 12, color = "black"),
    legend.title     = element_text(size = 12, color = "black"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA)
  )

carb_overlap_plot

# =============================================================================
# 7. Combine and save
# =============================================================================
combined_go_plot <- plot_grid(
  go_curated_plot, carb_overlap_plot,
  labels = c("A)", "B)"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(2, 1)
)

combined_go_plot

ggsave(plot = combined_go_plot,
       filename = file.path(FIG_DIR, paste0("Pchry_GO_combined_v2", DATE, ".pdf")),
       width = 10, height = 12, units = "in", dpi = 300)
ggsave(plot = combined_go_plot,
       filename = file.path(FIG_DIR, paste0("Pchry_GO_combined_", DATE, ".tiff")),
       width = 10, height = 14, units = "in", dpi = 200, bg = "white")