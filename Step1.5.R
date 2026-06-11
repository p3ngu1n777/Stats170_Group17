# =============================================================================
#  Standalone Gene Annotation & Pathway Analysis (KEGG + GO)
#  For Human Lung Cell Atlas limma-voom Results
# =============================================================================
#
# USAGE:
#   1. Run your original HLCA limma-voom analysis to completion
#   2. Run this script in a NEW R session
#   3. Update PATHS section below to point to your OUTPUT_DIR
#   4. Source this script or run line-by-line
#
# OUTPUTS:
#   - Annotated results CSVs (Ensembl, KEGG, and GO terms)
#   - Top 20 gene lists (simple format for manual searches)
#   - Summary reports (focused on significant genes, p < 0.05)
#   - Publication-ready plots (Stacked bars for Up/Down regulation)
#

# =============================================================================
# CONFIGURATION & PATHS
# =============================================================================

# UPDATE THIS PATH to your OUTPUT_DIR from the original analysis
OUTPUT_DIR <- "/Users/nikokorvink/Stats170AB/HLCA_NB_GLM_results"

# Define paths to input files from your original analysis
INPUT_FILES <- list(
  global_age = file.path(OUTPUT_DIR, "global_age_results_limma.csv"),
  sig_global_age = file.path(OUTPUT_DIR, "sig_global_age_results_limma.csv"),
  per_celltype = file.path(OUTPUT_DIR, "per_celltype_age_results_limma.csv")
)

# Create output subdirectories
ANNOT_DIR <- file.path(OUTPUT_DIR, "ANNOTATION_RESULTS")
PLOTS_DIR <- file.path(ANNOT_DIR, "plots")
TOP20_DIR <- file.path(ANNOT_DIR, "top_20_genes")
dir.create(ANNOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TOP20_DIR, showWarnings = FALSE, recursive = TRUE)

message("\n", strrep("=", 80))
message(" GENE ANNOTATION & PATHWAY ANALYSIS PIPELINE")
message(strrep("=", 80))
message("\nOutput directory: ", ANNOT_DIR)

# =============================================================================
# 1. PACKAGES
# =============================================================================

message("\n[1/8] Loading packages...")

annotation_pkgs <- c("biomaRt", "clusterProfiler", "org.Hs.eg.db")
data_pkgs <- c("dplyr", "tibble", "tidyr", "ggplot2", "ggrepel", 
               "patchwork", "pheatmap", "scales", "RColorBrewer",
               "stringr", "forcats") 
all_pkgs <- c(annotation_pkgs, data_pkgs)

missing_pkgs <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(biomaRt); library(clusterProfiler); library(org.Hs.eg.db)
  library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(ggrepel); library(patchwork); library(pheatmap)
  library(scales); library(RColorBrewer); library(stringr); library(forcats)
})

message("✓ All packages loaded")

# =============================================================================
# 2. LOAD INPUT FILES
# =============================================================================

message("\n[2/8] Loading results from original analysis...")

# Check that files exist
missing_files <- names(INPUT_FILES)[!file.exists(unlist(INPUT_FILES))]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "),
       "\nCheck OUTPUT_DIR path in configuration section")
}

results_global <- read.csv(INPUT_FILES$global_age, stringsAsFactors = FALSE) %>% as_tibble()
results_sig <- read.csv(INPUT_FILES$sig_global_age, stringsAsFactors = FALSE) %>% as_tibble()

results_celltype <- NULL
if (file.exists(INPUT_FILES$per_celltype)) {
  results_celltype <- read.csv(INPUT_FILES$per_celltype, stringsAsFactors = FALSE) %>% as_tibble()
  message("  ✓ Loaded global results: ", nrow(results_global), " genes")
  message("  ✓ Loaded significant results: ", nrow(results_sig), " genes")
  message("  ✓ Loaded per-celltype results: ", nrow(results_celltype), " genes")
} else {
  message("  ✓ Loaded global results: ", nrow(results_global), " genes")
  message("  ✓ Loaded significant results: ", nrow(results_sig), " genes")
  message("  ⚠ Per-celltype file not found (skipping)")
}

# =============================================================================
# 3. GENE ANNOTATION FUNCTIONS
# =============================================================================

message("\n[3/8] Setting up annotation functions...")

# ── 3A: Get gene info from Ensembl ────────────────────────────────────────
get_gene_info <- function(gene_symbols, verbose = FALSE) {
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
  empty_return <- tibble(Gene = character(), Description = character(), 
                         GeneType = character(), EntrezID = character(), 
                         EnsemblID = character())
  
  if (length(gene_symbols) == 0) return(empty_return)
  
  tryCatch({
    if (verbose) message("    Querying Ensembl for ", length(gene_symbols), " genes...")
    ensembl <- tryCatch({
      useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                 mirror = "useast", verbose = FALSE)
    }, error = function(e1) {
      if (verbose) message("    ⚠ US East mirror down, trying Asia mirror...")
      tryCatch({
        useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                   mirror = "asia", verbose = FALSE)
      }, error = function(e2) {
        if (verbose) message("    ⚠ Asia mirror down, trying Main site...")
        useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                   mirror = "www", verbose = FALSE)
      })
    })
    
    gene_info <- getBM(
      attributes = c("external_gene_name", "description", "gene_biotype", 
                     "entrezgene_id", "ensembl_gene_id"),
      filters = "external_gene_name",
      values = gene_symbols,
      mart = ensembl,
      verbose = FALSE,
      useCache = TRUE
    ) %>%
      as_tibble() %>%
      dplyr::rename(Gene = external_gene_name, 
                    Description = description,
                    GeneType = gene_biotype,
                    EntrezID = entrezgene_id,
                    EnsemblID = ensembl_gene_id) %>%
      dplyr::mutate(Description = gsub("\\s*\\[.*?\\]$", "", Description))
    
    return(gene_info)
  }, error = function(e) {
    message("    ⚠ All Ensembl mirrors failed: ", conditionMessage(e))
    return(empty_return)
  })
}

# ── 3B: Map genes to KEGG pathways ────────────────────────────────────────
get_kegg_pathways <- function(gene_symbols, verbose = FALSE) {
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
  empty_return <- tibble(Gene = character(), KEGG_ID = character(), 
                         KEGG_Name = character(), p.adjust = numeric())
  
  if (length(gene_symbols) == 0) return(empty_return)
  
  tryCatch({
    if (verbose) message("    Mapping ", length(gene_symbols), " genes to KEGG...")
    gene_to_entrez <- suppressWarnings(
      bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
           OrgDb = org.Hs.eg.db, drop = FALSE) %>% as_tibble()
    )
    if (nrow(gene_to_entrez) == 0) return(empty_return)
    
    if (verbose) message("    Querying KEGG database...")
    kegg_result <- tryCatch({
      enrichKEGG(gene = gene_to_entrez$ENTREZID,
                 organism = "hsa", pvalueCutoff = 1, pAdjustMethod = "BH") %>%
        as.data.frame() %>% as_tibble() %>%
        dplyr::rename(KEGG_ID = ID, KEGG_Name = Description) %>%
        dplyr::select(geneID, KEGG_ID, KEGG_Name, p.adjust) %>%
        tidyr::separate_rows(geneID, sep = "/") %>%
        dplyr::rename(ENTREZID = geneID)
    }, error = function(e) {
      return(tibble(ENTREZID = character(), KEGG_ID = character(), 
                    KEGG_Name = character(), p.adjust = numeric()))
    })
    
    if (nrow(kegg_result) == 0) return(empty_return)
    
    kegg_final <- kegg_result %>%
      dplyr::left_join(gene_to_entrez, by = c("ENTREZID" = "ENTREZID")) %>%
      dplyr::select(Gene = SYMBOL, KEGG_ID, KEGG_Name, p.adjust) %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(Gene))
    
    return(kegg_final)
  }, error = function(e) {
    message("    ⚠ KEGG query failed: ", conditionMessage(e))
    return(empty_return)
  })
}

# ── 3C: Map genes to GO Biological Processes ──────────────────────────────
get_go_terms <- function(gene_symbols, verbose = FALSE) {
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
  empty_return <- tibble(Gene = character(), GO_ID = character(), 
                         GO_Name = character(), p.adjust_GO = numeric())
  
  if (length(gene_symbols) == 0) return(empty_return)
  
  tryCatch({
    if (verbose) message("    Mapping ", length(gene_symbols), " genes to GO Biological Processes...")
    gene_to_entrez <- suppressWarnings(
      bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
           OrgDb = org.Hs.eg.db, drop = FALSE) %>% as_tibble()
    )
    if (nrow(gene_to_entrez) == 0) return(empty_return)
    
    if (verbose) message("    Querying GO database...")
    go_result <- tryCatch({
      enrichGO(gene = gene_to_entrez$ENTREZID, OrgDb = org.Hs.eg.db, 
               ont = "BP", pvalueCutoff = 1, pAdjustMethod = "BH") %>%
        as.data.frame() %>% as_tibble() %>%
        dplyr::rename(GO_ID = ID, GO_Name = Description, p.adjust_GO = p.adjust) %>%
        dplyr::select(geneID, GO_ID, GO_Name, p.adjust_GO) %>%
        tidyr::separate_rows(geneID, sep = "/") %>%
        dplyr::rename(ENTREZID = geneID)
    }, error = function(e) {
      return(tibble(ENTREZID = character(), GO_ID = character(), 
                    GO_Name = character(), p.adjust_GO = numeric()))
    })
    
    if (nrow(go_result) == 0) return(empty_return)
    
    go_final <- go_result %>%
      dplyr::left_join(gene_to_entrez, by = c("ENTREZID" = "ENTREZID")) %>%
      dplyr::select(Gene = SYMBOL, GO_ID, GO_Name, p.adjust_GO) %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(Gene))
    
    return(go_final)
  }, error = function(e) {
    message("    ⚠ GO query failed: ", conditionMessage(e))
    return(empty_return)
  })
}

# ── 3D: Main annotation function ──────────────────────────────────────────
annotate_results <- function(results_df, gene_col = "Gene", verbose = FALSE) {
  gene_list <- unique(results_df[[gene_col]])
  
  if (verbose) message("  Getting descriptions...")
  gene_info <- get_gene_info(gene_list, verbose = verbose)
  
  if (verbose) message("  Getting KEGG pathway information...")
  kegg_path <- get_kegg_pathways(gene_list, verbose = verbose)
  
  if (verbose) message("  Getting GO Biological Process information...")
  go_path <- get_go_terms(gene_list, verbose = verbose)
  
  pathway_summary <- kegg_path %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      KEGG_Pathways = paste(KEGG_Name, collapse = "; "),
      KEGG_IDs = paste(KEGG_ID, collapse = "; "),
      N_Pathways = dplyr::n_distinct(KEGG_ID),
      .groups = "drop"
    )
  
  go_summary <- go_path %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      GO_Terms = paste(GO_Name, collapse = "; "),
      GO_IDs = paste(GO_ID, collapse = "; "),
      N_GO_Terms = dplyr::n_distinct(GO_ID),
      .groups = "drop"
    )
  
  annotated <- results_df %>%
    dplyr::left_join(gene_info %>% dplyr::select(Gene, Description, GeneType, EntrezID), 
                     by = gene_col) %>%
    dplyr::left_join(pathway_summary, by = gene_col) %>%
    dplyr::left_join(go_summary, by = gene_col) %>%
    dplyr::mutate(
      KEGG_Pathways = tidyr::replace_na(KEGG_Pathways, "No KEGG pathway"),
      N_Pathways = tidyr::replace_na(N_Pathways, 0),
      GO_Terms = tidyr::replace_na(GO_Terms, "No GO term"),
      N_GO_Terms = tidyr::replace_na(N_GO_Terms, 0)
    )
  
  return(annotated)
}

message("✓ Annotation functions ready")

# =============================================================================
# 4. ANNOTATE GLOBAL RESULTS
# =============================================================================

message("\n[4/8] Annotating global age results...")

results_global_annotated <- annotate_results(results_global, gene_col = "Gene", verbose = TRUE)

write.csv(results_global_annotated,
          file.path(ANNOT_DIR, "01_global_age_results_ANNOTATED.csv"),
          row.names = FALSE)

message("✓ Global results annotated")

# =============================================================================
# 5. ANNOTATE PER-CELLTYPE RESULTS (if available)
# =============================================================================

if (!is.null(results_celltype)) {
  message("\n[5/8] Annotating per-celltype results...")
  
  results_celltype_annotated <- annotate_results(results_celltype, gene_col = "Gene", 
                                                 verbose = TRUE)
  
  write.csv(results_celltype_annotated,
            file.path(ANNOT_DIR, "02_per_celltype_results_ANNOTATED.csv"),
            row.names = FALSE)
  
  message("✓ Per-celltype results annotated")
} else {
  message("\n[5/8] Skipping per-celltype annotation (no file found)")
  results_celltype_annotated <- NULL
}

# =============================================================================
# 6. EXTRACT TOP 20 GENES FOR MANUAL SEARCHES
# =============================================================================

message("\n[6/8] Extracting top 20 genes from each comparison...")

# ── 6A: Top 20 from global analysis ──────────────────────────────────────
top20_global <- results_global_annotated %>%
  dplyr::arrange(FDR) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(dplyr::any_of(c("Gene", "logFC", "FDR", "PValue", "Description", 
                                "N_Pathways", "KEGG_Pathways", "N_GO_Terms", "GO_Terms")))

write.csv(top20_global, file.path(TOP20_DIR, "01_top20_global_age.csv"), row.names = FALSE)
writeLines(top20_global %>% dplyr::pull(Gene), 
           file.path(TOP20_DIR, "01_top20_global_age_SIMPLE_LIST.txt"))

message("  ✓ Top 20 global genes extracted")

# ── 6B: Top 20 from significant results ──────────────────────────────────
if (nrow(results_sig) > 0) {
  top20_sig <- results_global_annotated %>%
    dplyr::filter(Gene %in% results_sig$Gene) %>%
    dplyr::arrange(FDR) %>%
    dplyr::slice_head(n = min(20, nrow(results_sig))) %>%
    dplyr::select(dplyr::any_of(c("Gene", "logFC", "FDR", "PValue", "Description", 
                                  "N_Pathways", "KEGG_Pathways", "N_GO_Terms", "GO_Terms")))
  
  write.csv(top20_sig, file.path(TOP20_DIR, "02_top20_significant_genes.csv"), row.names = FALSE)
  writeLines(top20_sig %>% dplyr::pull(Gene), 
             file.path(TOP20_DIR, "02_top20_significant_genes_SIMPLE_LIST.txt"))
  
  message("  ✓ Top 20 significant genes extracted")
}

# ── 6C: Top 20 per cell type (if available) ──────────────────────────────
if (!is.null(results_celltype_annotated)) {
  celltypes <- unique(results_celltype_annotated$CellType)
  for (ct in celltypes) {
    ct_safe <- gsub("[^a-zA-Z0-9_]", "_", ct)
    ct_data <- results_celltype_annotated %>%
      dplyr::filter(CellType == ct) %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::select(dplyr::any_of(c("Gene", "logFC", "FDR", "PValue", "Description", 
                                    "N_Pathways", "KEGG_Pathways", "N_GO_Terms", "GO_Terms")))
    
    if (nrow(ct_data) > 0) {
      write.csv(ct_data, file.path(TOP20_DIR, paste0("03_top20_", ct_safe, ".csv")), row.names = FALSE)
      writeLines(ct_data %>% dplyr::pull(Gene), 
                 file.path(TOP20_DIR, paste0("03_top20_", ct_safe, "_SIMPLE_LIST.txt")))
    }
  }
  message("  ✓ Top 20 genes per cell type extracted (", length(celltypes), " cell types)")
}

# =============================================================================
# 7. GENERATE SUMMARY REPORTS
# =============================================================================

message("\n[7/8] Generating summary reports...")

# Filter strictly for raw P-Value < 0.05 and assign biological direction
sig_raw_p <- results_global_annotated %>% 
  dplyr::filter(PValue < 0.05) %>%
  dplyr::mutate(Direction = ifelse(logFC > 0, "Up", "Down"))

# ── 7A: KEGG Pathway frequency (Split Counts) ────────────────────────────
if (nrow(sig_raw_p) > 0) {
  # Get top 30 pathways by total count
  top_kegg <- sig_raw_p %>%
    dplyr::filter(N_Pathways > 0) %>% 
    tidyr::separate_rows(KEGG_Pathways, sep = "; ") %>%
    dplyr::count(KEGG_Pathways, name = "total", sort = TRUE) %>%
    dplyr::slice_head(n = 30)
  
  # Split those top 30 into Up/Down subsets
  pathway_freq_sig <- sig_raw_p %>%
    dplyr::filter(N_Pathways > 0) %>% 
    tidyr::separate_rows(KEGG_Pathways, sep = "; ") %>%
    dplyr::filter(KEGG_Pathways %in% top_kegg$KEGG_Pathways) %>%
    dplyr::count(KEGG_Pathways, Direction) %>%
    dplyr::left_join(top_kegg, by = "KEGG_Pathways") %>%
    dplyr::arrange(desc(total), Direction)
  
  write.csv(pathway_freq_sig, file.path(ANNOT_DIR, "03_top30_KEGG_pathways_P05_Split.csv"), row.names = FALSE)
  message("  ✓ Significant KEGG frequencies calculated")
}

# ── 7B: GO Term frequency (Split Counts) ─────────────────────────────────
if (nrow(sig_raw_p) > 0) {
  # Get top 30 GO terms by total count
  top_go <- sig_raw_p %>%
    dplyr::filter(N_GO_Terms > 0) %>% 
    tidyr::separate_rows(GO_Terms, sep = "; ") %>%
    dplyr::count(GO_Terms, name = "total", sort = TRUE) %>%
    dplyr::slice_head(n = 30)
  
  # Split those top 30 into Up/Down subsets
  go_freq_sig <- sig_raw_p %>%
    dplyr::filter(N_GO_Terms > 0) %>% 
    tidyr::separate_rows(GO_Terms, sep = "; ") %>%
    dplyr::filter(GO_Terms %in% top_go$GO_Terms) %>%
    dplyr::count(GO_Terms, Direction) %>%
    dplyr::left_join(top_go, by = "GO_Terms") %>%
    dplyr::arrange(desc(total), Direction)
  
  write.csv(go_freq_sig, file.path(ANNOT_DIR, "03_top30_GO_Terms_P05_Split.csv"), row.names = FALSE)
  message("  ✓ Significant GO Biological Process frequencies calculated")
}

# ── 7C: Gene type distribution ──────────────────────────────────────────
if ("GeneType" %in% colnames(results_global_annotated)) {
  genetype_dist <- results_global_annotated %>%
    dplyr::count(GeneType, sort = TRUE) %>%
    dplyr::mutate(Proportion = round(n / nrow(results_global_annotated), 4),
                  Percentage = paste0(round(Proportion * 100, 1), "%"))
  write.csv(genetype_dist, file.path(ANNOT_DIR, "04_gene_type_distribution.csv"), row.names = FALSE)
}

# =============================================================================
# 8. PUBLICATION-READY PLOTS
# =============================================================================

message("\n[8/8] Creating plots...")

safe_save_plot <- function(plot_obj, base_name) {
  tryCatch({
    ggsave(file.path(PLOTS_DIR, paste0(base_name, ".pdf")), plot_obj, width = 10, height = 7, dpi = 300)
    ggsave(file.path(PLOTS_DIR, paste0(base_name, ".png")), plot_obj, width = 10, height = 7, dpi = 200)
  }, error = function(e) message("  ⚠ Failed to save ", base_name))
}

# ── PLOT 1: Annotated Volcano ───────────────────────────────────────────
p_volcano <- results_global_annotated %>%
  dplyr::arrange(FDR) %>%
  dplyr::mutate(
    Sig = factor(dplyr::case_when(
      FDR < 0.05 & logFC > 1 ~ "Up (FDR<0.05, |logFC|>1)",
      FDR < 0.05 & logFC < -1 ~ "Down (FDR<0.05, |logFC|>1)",
      FDR < 0.05 ~ "Significant (FDR<0.05)",
      TRUE ~ "Not significant"
    ), levels = c("Up (FDR<0.05, |logFC|>1)", "Down (FDR<0.05, |logFC|>1)", 
                  "Significant (FDR<0.05)", "Not significant")),
    TopGene = ifelse(dplyr::row_number() <= 10, Gene, NA)
  ) %>%
  ggplot(aes(logFC, -log10(PValue), colour = Sig)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values = c("#d62728", "#1f77b4", "#ff7f0e", "grey75")) +
  geom_hline(yintercept = -log10(0.05 / nrow(results_global_annotated)), linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 0.4, colour = "grey50") +
  geom_text_repel(aes(label = TopGene), size = 3, colour = "black", max.overlaps = 15, box.padding = 0.3) +
  labs(title = "Age-Associated Genes: Volcano Plot", subtitle = "Old (65+) vs Young (5-40) years",
       x = expression(log[2](FC)), y = expression(-log[10](p-value)), colour = "Classification") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

safe_save_plot(p_volcano, "01_volcano_annotated")

# ── PLOT 2: KEGG Enrichment (Stacked Bar, Up vs Down) ───────────────────
if (exists("pathway_freq_sig") && nrow(pathway_freq_sig) > 0) {
  
  # Select Top 15 pathways
  top15_kegg_names <- unique(pathway_freq_sig$KEGG_Pathways)[1:min(15, length(unique(pathway_freq_sig$KEGG_Pathways)))]
  
  plot_data_kegg <- pathway_freq_sig %>%
    dplyr::filter(KEGG_Pathways %in% top15_kegg_names) %>%
    dplyr::mutate(KEGG_Pathways = stringr::str_trunc(KEGG_Pathways, 45),
                  KEGG_Pathways = forcats::fct_reorder(KEGG_Pathways, total))
  
  # Separate dataframe for the total count labels
  label_data_kegg <- plot_data_kegg %>%
    dplyr::select(KEGG_Pathways, total) %>%
    dplyr::distinct()
  
  p_kegg <- ggplot(plot_data_kegg, aes(x = n, y = KEGG_Pathways, fill = Direction)) +
    geom_col(color = "black", linewidth = 0.2, alpha = 0.9) +
    geom_text(data = label_data_kegg, aes(x = total, y = KEGG_Pathways, label = total), 
              inherit.aes = FALSE, hjust = -0.2, size = 3) +
    scale_fill_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) + # Prevents text cutoff
    labs(title = "Most Common KEGG Pathways", 
         subtitle = "Among significant genes (Raw p < 0.05), split by expression direction",
         x = "Number of significant genes", y = NULL) +
    theme_bw(base_size = 11) + 
    theme(axis.text.x = element_text(size = 10))
  
  safe_save_plot(p_kegg, "02A_KEGG_pathways_bar_SPLIT")
}

# ── PLOT 3: GO Term Enrichment (Stacked Bar, Up vs Down) ────────────────
if (exists("go_freq_sig") && nrow(go_freq_sig) > 0) {
  
  # Select Top 15 GO terms
  top15_go_names <- unique(go_freq_sig$GO_Terms)[1:min(15, length(unique(go_freq_sig$GO_Terms)))]
  
  plot_data_go <- go_freq_sig %>%
    dplyr::filter(GO_Terms %in% top15_go_names) %>%
    dplyr::mutate(GO_Terms = stringr::str_trunc(GO_Terms, 45),
                  GO_Terms = forcats::fct_reorder(GO_Terms, total))
  
  label_data_go <- plot_data_go %>%
    dplyr::select(GO_Terms, total) %>%
    dplyr::distinct()
  
  p_go <- ggplot(plot_data_go, aes(x = n, y = GO_Terms, fill = Direction)) +
    geom_col(color = "black", linewidth = 0.2, alpha = 0.9) +
    geom_text(data = label_data_go, aes(x = total, y = GO_Terms, label = total), 
              inherit.aes = FALSE, hjust = -0.2, size = 3) +
    scale_fill_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Most Common GO Biological Processes", 
         subtitle = "Among significant genes (Raw p < 0.05), split by expression direction",
         x = "Number of significant genes", y = NULL) +
    theme_bw(base_size = 11) + 
    theme(axis.text.x = element_text(size = 10))
  
  safe_save_plot(p_go, "02B_GO_Terms_bar_SPLIT")
}

# ── PLOT 4: Heatmap of top genes ────────────────────────────────────────
top_genes_hm <- results_global_annotated %>%
  dplyr::arrange(FDR) %>% dplyr::slice_head(n = 30) %>% dplyr::pull(Gene)

hm_data <- results_global_annotated %>%
  dplyr::filter(Gene %in% top_genes_hm) %>%
  dplyr::arrange(match(Gene, top_genes_hm)) %>%
  dplyr::select(dplyr::any_of(c("Gene", "logFC", "FDR", "Description"))) %>%
  dplyr::mutate(Gene_label = paste0(Gene, "\n(FDR=", format(FDR, digits = 2), ")"))

hm_matrix <- matrix(hm_data$logFC, nrow = 1, dimnames = list("logFC", hm_data$Gene_label))
abs_max <- max(abs(hm_matrix))

pdf(file.path(PLOTS_DIR, "03_top30_genes_heatmap.pdf"), width = 14, height = 4)
pheatmap(hm_matrix, color = colorRampPalette(c("#1f77b4", "white", "#d62728"))(100),
         breaks = seq(-abs_max, abs_max, length.out = 101),
         cluster_cols = FALSE, cluster_rows = FALSE, fontsize_col = 7, fontsize_row = 10,
         main = "Top 30 Age-Associated Genes (logFC, Old vs Young)")
invisible(dev.off())

png(file.path(PLOTS_DIR, "03_top30_genes_heatmap.png"), width = 1400, height = 400, res = 100)
pheatmap(hm_matrix, color = colorRampPalette(c("#1f77b4", "white", "#d62728"))(100),
         breaks = seq(-abs_max, abs_max, length.out = 101),
         cluster_cols = FALSE, cluster_rows = FALSE, fontsize_col = 7, fontsize_row = 10,
         main = "Top 30 Age-Associated Genes (logFC, Old vs Young)")
invisible(dev.off())

message("✓ Plots generated")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

message("\n", strrep("=", 80))
message(" ANALYSIS COMPLETE")
message(strrep("=", 80))
message("\n📁 All files saved to: ", ANNOT_DIR)
message("\n", strrep("=", 80), "\n")