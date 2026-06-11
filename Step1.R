# =============================================================================
#  Human Lung Cell Atlas  —  limma-voom Engine (V14 — Temporal Contrasts)
# =============================================================================

if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
  Sys.setenv(R_MAX_VSIZE = "64Gb")
}

# ── 1. Packages ────────────────────────────────────────────────────────────────
bioc_pkgs <- c("zellkonverter", "SingleCellExperiment", "SummarizedExperiment", 
               "edgeR", "limma", "BiocParallel")
cran_pkgs <- c("dplyr", "tibble", "tidyr", "ggplot2", "patchwork", 
               "pheatmap", "scales", "RColorBrewer", "ggrepel", "Matrix")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("BPCells", quietly = TRUE)) {
  message("Installing BPCells…")
  install.packages("BPCells", repos = c("https://bnprks.r-universe.dev", "https://cloud.r-project.org"))
}

suppressPackageStartupMessages({
  library(BPCells); library(zellkonverter); library(SingleCellExperiment)
  library(SummarizedExperiment); library(edgeR); library(limma); library(dplyr)
  library(tibble); library(tidyr); library(ggplot2); library(patchwork)
  library(ggrepel); library(pheatmap); library(scales); library(RColorBrewer)
  library(Matrix); library(BiocParallel)
})
register(SerialParam())

# ── 2. Paths & Config ──────────────────────────────────────────────────────────
H5AD_PATH    <- "Stats170AB/dbb5ad81-1713-4aee-8257-396fbabe7c6e.h5ad"
OUTPUT_DIR   <- file.path(dirname(H5AD_PATH), "HLCA_NB_GLM_results")
PLOTS_DIR    <- file.path(OUTPUT_DIR, "plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

DONOR_COL    <- "donor_id"
CELLTYPE_COL <- "ann_level_3"
CELLCLASS_COL<- "ann_level_1"  
AGE_COL      <- "age_or_mean_of_age_range"
SEX_COL      <- "sex"
SMOKING_COL  <- "smoking_status"
DATASET_COL  <- "dataset"

# ── 3. Load metadata & extract gene symbols ────────────────────────────────────
message("\n[1/7] Loading metadata from h5ad…")
sce_meta <- readH5AD(H5AD_PATH, use_hdf5 = TRUE, reader = "R", skip_obsm = TRUE, skip_obsp = TRUE, verbose = FALSE)
actual_cols <- names(colData(sce_meta))
LIB_COL <- intersect(c("total_counts", "n_counts", "sum"), actual_cols)[1]

rd <- as.data.frame(rowData(sce_meta))
sym_col <- intersect(c("feature_name", "gene_name", "Symbol", "symbol", "name", "feature_symbol"), names(rd))[1]
if (!is.na(sym_col)) {
  gene_names <- as.character(rd[[sym_col]])
  gene_names[is.na(gene_names) | gene_names == ""] <- rownames(sce_meta)[is.na(gene_names) | gene_names == ""]
  gene_symbols_unique <- make.unique(gene_names)
} else {
  gene_symbols_unique <- rownames(sce_meta)
}

# ── 4. Derive AgeGroup & Overwrite Unknowns ────────────────────────────────────
message("[2/7] Deriving AgeGroup and Fixing 'Unknown/None' cell types…")

age_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(sce_meta[[AGE_COL]]))))
sce_meta$age_numeric <- age_numeric
sce_meta$AgeGroup <- cut(age_numeric, breaks = c(4, 40, 65, Inf),
                         labels = c("Young (5-40)", "Middle-aged (41-65)", "Old (65+)"), right = TRUE)

cell_types   <- as.character(sce_meta[[CELLTYPE_COL]])
cell_classes <- as.character(sce_meta[[CELLCLASS_COL]]) 

bad_labels     <- c("unknown", "none")
is_problematic <- tolower(cell_types) %in% bad_labels | is.na(cell_types)

cell_types[is_problematic] <- paste0("Unknown ", cell_classes[is_problematic])
sce_meta[[CELLTYPE_COL]] <- cell_types

valid_cells <- !is.na(sce_meta$AgeGroup)
sce_meta    <- sce_meta[, valid_cells]

message("      Retained ", ncol(sce_meta), " valid cells. Unknowns successfully grouped by lineage.")

# ── 5. Pseudo-bulk aggregation via BPCells ────────────────────────────────────
message("[3/7] Pseudo-bulk aggregation with BPCells…")
bp_mat <- open_matrix_anndata_hdf5(H5AD_PATH)
bp_mat <- bp_mat[, colnames(sce_meta)]

cd <- as.data.frame(colData(sce_meta))
group_ids <- paste(cd[[DONOR_COL]], cd[[CELLTYPE_COL]], sep = "__")
group_fac <- factor(group_ids)

indicator <- sparseMatrix(
  i = seq_len(ncol(bp_mat)), j = as.integer(group_fac), x = 1L,
  dims = c(ncol(bp_mat), nlevels(group_fac)), dimnames = list(NULL, levels(group_fac))
)

counts_pb <- as.matrix(bp_mat %*% indicator)
rownames(counts_pb) <- gene_symbols_unique

# ── 6. Build matching metadata ─────────────────────────────────────────────────
message("      Building pseudo-bulk metadata…")
meta_pb_raw <- cd %>%
  mutate(group_id = group_ids) %>%
  group_by(group_id) %>%
  summarise(
    !!DONOR_COL    := first(.data[[DONOR_COL]]),
    !!CELLTYPE_COL := first(.data[[CELLTYPE_COL]]),
    !!CELLCLASS_COL:= first(.data[[CELLCLASS_COL]]),
    !!SEX_COL      := first(.data[[SEX_COL]]),
    !!SMOKING_COL  := first(.data[[SMOKING_COL]]),
    !!DATASET_COL  := first(.data[[DATASET_COL]]),
    AgeGroup       = first(AgeGroup),
    lib_size       = sum(.data[[LIB_COL]], na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  slice(match(colnames(counts_pb), group_id)) %>%
  as.data.frame() %>%
  column_to_rownames("group_id")

rm(bp_mat, indicator, cd, group_ids, group_fac, sce_meta); gc(verbose = FALSE)

# ── 7. Filter genes & samples ─────────────────────────────────────────────────
message("[4/7] Filtering & Removing NAs…")
covariates <- c(CELLTYPE_COL, SEX_COL, SMOKING_COL, DATASET_COL, "AgeGroup", "lib_size")
keep_samps <- complete.cases(meta_pb_raw[, covariates]) & meta_pb_raw$lib_size > 0

min_samps <- max(5, floor(sum(keep_samps) * 0.05))
keep_genes <- rowSums(counts_pb >= 10) >= min_samps

counts_filt <- counts_pb[keep_genes, keep_samps]
meta_pb <- meta_pb_raw[keep_samps, ]
rm(counts_pb); gc(verbose = FALSE)

# ── 8. Design matrix & DGEList ────────────────────────────────────────────────
message("[5/7] Design matrix and DGEList…")
meta_pb <- meta_pb %>%
  mutate(
    AgeGroup  = relevel(droplevels(factor(AgeGroup)), ref = "Young (5-40)"),
    CellType  = droplevels(factor(.data[[CELLTYPE_COL]])),
    CellClass = factor(.data[[CELLCLASS_COL]]),
    Sex       = relevel(factor(.data[[SEX_COL]]),      ref = "male"),
    Smoking   = relevel(factor(.data[[SMOKING_COL]]), ref = "never"),
    Dataset   = factor(.data[[DATASET_COL]])
  )

design <- model.matrix(~ AgeGroup * CellType + Sex + Smoking + Dataset, data = meta_pb)
qr_des <- qr(design)
if (qr_des$rank < ncol(design)) design <- design[, qr_des$pivot[seq_len(qr_des$rank)]]

dge <- DGEList(counts = counts_filt, lib.size = meta_pb$lib_size, samples = meta_pb)
dge <- calcNormFactors(dge, method = "TMM")

# ── 9. Global models via limma-voom (TEMPORAL PHASES) ─────────────────────────
message("[6/7] Fitting global temporal models via limma-voom…")
v <- voom(dge, design, plot = FALSE)
fit_global <- eBayes(lmFit(v, design))

# Identify exact column names in the design matrix
oldest_coef <- grep("^AgeGroupOld[^:]+$", colnames(design), value = TRUE)
mid_coef    <- grep("^AgeGroupMiddle-aged[^:]+$", colnames(design), value = TRUE)

# 9A. Total Aging (Old vs Young)
results_age <- topTable(fit_global, coef = oldest_coef, number = Inf) %>%
  as_tibble(rownames = "Gene") %>%
  rename(FDR = adj.P.Val, PValue = P.Value) %>%
  mutate(Test = "Old_vs_Young")

sig_results <- filter(results_age, FDR < 0.05)
write.csv(results_age, file.path(OUTPUT_DIR, "global_age_results_limma.csv"), row.names = FALSE)
write.csv(sig_results, file.path(OUTPUT_DIR, "sig_global_age_results_limma.csv"), row.names = FALSE)

# 9B. Early Aging (Middle vs Young)
results_early <- topTable(fit_global, coef = mid_coef, number = Inf) %>%
  as_tibble(rownames = "Gene") %>%
  rename(FDR = adj.P.Val, PValue = P.Value) %>%
  mutate(Test = "Mid_vs_Young")
write.csv(results_early, file.path(OUTPUT_DIR, "global_early_aging_results.csv"), row.names = FALSE)

# 9C. Late Aging (Old vs Middle) 
# Calculate the mathematical contrast: (Old_vs_Young) - (Mid_vs_Young) = Old_vs_Mid
cont_vec <- rep(0, ncol(design))
cont_vec[which(colnames(design) == oldest_coef)] <- 1
cont_vec[which(colnames(design) == mid_coef)] <- -1

fit_late <- contrasts.fit(fit_global, cont_vec)
fit_late <- eBayes(fit_late)

results_late <- topTable(fit_late, coef = 1, number = Inf) %>%
  as_tibble(rownames = "Gene") %>%
  rename(FDR = adj.P.Val, PValue = P.Value) %>%
  mutate(Test = "Old_vs_Mid")
write.csv(results_late, file.path(OUTPUT_DIR, "global_late_aging_results.csv"), row.names = FALSE)

# ── 10. Per-cell-type stratified age models ────────────────────────────────────
message("[7/7] Per-cell-type stratified age models…")
run_per_ct_limma <- function(ct_label) {
  idx <- which(meta_pb$CellType == ct_label)
  if (length(idx) < 8) return(NULL)
  
  sub_meta <- meta_pb[idx, ] %>% mutate(AgeGroup = factor(AgeGroup, levels = c("Young (5-40)", "Middle-aged (41-65)", "Old (65+)")))
  if (!("Old (65+)" %in% sub_meta$AgeGroup) | !("Young (5-40)" %in% sub_meta$AgeGroup)) return(NULL)
  
  dm <- tryCatch(model.matrix(~ AgeGroup + Sex + Smoking + Dataset, data = sub_meta), error = function(e) NULL)
  if (is.null(dm)) return(NULL)
  qr_dm <- qr(dm)
  if (qr_dm$rank < ncol(dm)) dm <- dm[, qr_dm$pivot[seq_len(qr_dm$rank)]]
  
  cnt <- counts_filt[, idx, drop = FALSE]
  keep <- rowSums(cnt >= 5) >= 3
  if (nrow(cnt[keep, , drop = FALSE]) < 50) return(NULL)
  
  dge_ct <- DGEList(counts = cnt[keep, , drop = FALSE], lib.size = sub_meta$lib_size)
  fit_ct <- eBayes(lmFit(voom(calcNormFactors(dge_ct), dm, plot = FALSE), dm))
  
  age_c <- grep("AgeGroupOld", colnames(dm), value = TRUE)
  if (length(age_c) == 0) return(NULL)
  
  topTable(fit_ct, coef = age_c, number = Inf) %>%
    as_tibble(rownames = "Gene") %>% rename(FDR = adj.P.Val, PValue = P.Value) %>%
    mutate(CellType = ct_label)
}

ct_res_df <- bind_rows(Filter(Negate(is.null), lapply(levels(meta_pb$CellType), run_per_ct_limma)))
if (nrow(ct_res_df) > 0) write.csv(ct_res_df, file.path(OUTPUT_DIR, "per_celltype_age_results_limma.csv"), row.names = FALSE)

# ── 11. Sample Size & Variance Table ──────────────────────────────────────────
message("\nGenerating Sample Size Table to justify variance...")
sample_counts <- meta_pb %>%
  count(AgeGroup, name = "Number_of_Biological_Replicates") %>%
  mutate(Notes = ifelse(Number_of_Biological_Replicates == min(Number_of_Biological_Replicates),
                        "Lowest n = highest standard error (wider bars)", "Stable variance"))
write.csv(sample_counts, file.path(OUTPUT_DIR, "05_age_group_sample_sizes.csv"), row.names = FALSE)

# =============================================================================
#  PLOTS
# =============================================================================
message("Generating plots…")
logcpm_mat <- v$E 

safe_save <- function(plot_obj, base_name, w, h) {
  tryCatch({
    ggsave(file.path(PLOTS_DIR, paste0(base_name, ".pdf")), plot_obj, width = w, height = h)
    ggsave(file.path(PLOTS_DIR, paste0(base_name, ".png")), plot_obj, width = w, height = h, dpi = 200)
  }, error = function(e) message("      ⚠ Failed to save ", base_name))
}

# Plot 1: Volcano
if (nrow(results_age) > 0) {
  p_volcano <- results_age %>%
    mutate(Sig = factor(case_when(
      FDR < 0.05 & logFC > 1 ~ "Up (FDR<0.05, |logFC|>1)",
      FDR < 0.05 & logFC < -1 ~ "Down (FDR<0.05, |logFC|>1)",
      FDR < 0.05 ~ "Significant (FDR<0.05)", TRUE ~ "Not significant"
    ), levels = c("Up (FDR<0.05, |logFC|>1)", "Down (FDR<0.05, |logFC|>1)", "Significant (FDR<0.05)", "Not significant"))) %>%
    ggplot(aes(logFC, -log10(PValue), colour = Sig)) + geom_point(size = 0.6, alpha = 0.5) +
    scale_colour_manual(values = c("#d62728", "#1f77b4", "#ff7f0e", "grey75")) +
    geom_hline(yintercept = -log10(0.05 / nrow(results_age)), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 0.4) +
    geom_text_repel(data = . %>% arrange(FDR) %>% slice_head(n = 10), aes(label = Gene), size = 3, colour = "black") +
    labs(title = "Volcano - age effect (Old vs Young)", x = "log2 FC (Old vs Young)", y = expression(-log[10](p-value)), colour = NULL) +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
  safe_save(p_volcano, "01_volcano_age_main", 9, 7)
}

# Plot 2: Bar chart
if (nrow(ct_res_df) > 0 && sum(ct_res_df$FDR < 0.05) > 0) {
  bar_data <- ct_res_df %>% filter(FDR < 0.05) %>%
    mutate(Direction = ifelse(logFC > 0, "Up in Old", "Down in Old")) %>%
    count(CellType, Direction) %>% group_by(CellType) %>% mutate(total = sum(n)) %>% ungroup() %>%
    arrange(desc(total)) %>% mutate(CellType = factor(CellType, levels = unique(CellType)))
  
  p_bar <- ggplot(bar_data, aes(CellType, n, fill = Direction)) + geom_col(width = 0.7) +
    scale_fill_manual(values = c("Up in Old" = "#d62728", "Down in Old" = "#1f77b4")) + coord_flip() +
    labs(title = "Significant Age-Associated Genes per Cell Type", x = NULL, y = "Number of DE genes") +
    theme_bw(base_size = 12) + theme(legend.position = "top")
  safe_save(p_bar, "02_bar_DE_per_celltype", 9, max(5, length(unique(bar_data$CellType)) * 0.35 + 2))
}

# Plot 3: Annotated Heatmap
if (nrow(ct_res_df) > 0 && nrow(results_age) > 0) {
  top_genes_hm <- results_age %>% arrange(FDR) %>% slice_head(n = 40) %>% pull(Gene)
  
  hm_mat <- ct_res_df %>% 
    dplyr::filter(Gene %in% top_genes_hm) %>% 
    dplyr::select(Gene, CellType, logFC) %>%
    pivot_wider(names_from = CellType, values_from = logFC, values_fn = mean, values_fill = 0) %>%
    column_to_rownames("Gene") %>% 
    as.matrix()
  
  if (nrow(hm_mat) >= 2 && ncol(hm_mat) >= 2) {
    row_ord <- results_age %>% 
      dplyr::filter(Gene %in% rownames(hm_mat)) %>% 
      arrange(desc(logFC)) %>% 
      pull(Gene)
    
    hm_mat <- hm_mat[row_ord[row_ord %in% rownames(hm_mat)], , drop = FALSE]
    abs_m <- max(abs(hm_mat), na.rm = TRUE)
    
    # Creates the annotation dataframe mapping ann_level_3 to ann_level_1 (with ties broken)
    annot_df <- meta_pb %>%
      dplyr::select(CellType, CellClass) %>%
      group_by(CellType) %>%
      slice_head(n = 1) %>% 
      ungroup() %>%
      rename(Lineage = CellClass) %>%
      as.data.frame() %>%
      remove_rownames() %>%
      column_to_rownames("CellType")
    
    try({
      pdf(file.path(PLOTS_DIR, "03_heatmap_genes_x_celltype.pdf"), 
          width = max(10, ncol(hm_mat) * 0.5 + 3), 
          height = max(8, nrow(hm_mat) * 0.28 + 3))
      
      pheatmap(hm_mat, 
               color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), 
               breaks         = seq(-abs_m, abs_m, length.out = 101),
               annotation_col = annot_df, 
               cluster_rows   = FALSE, 
               cluster_cols   = TRUE, 
               fontsize_row   = 8, 
               angle_col      = 45, 
               main           = "Age logFC (Old vs Young): top genes x cell type")
      dev.off()
    })
    
    try({
      png(file.path(PLOTS_DIR, "03_heatmap_genes_x_celltype.png"), 
          width = max(10, ncol(hm_mat) * 0.5 + 3) * 100, 
          height = max(8, nrow(hm_mat) * 0.28 + 3) * 100, 
          res = 100)
      
      pheatmap(hm_mat, 
               color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), 
               breaks         = seq(-abs_m, abs_m, length.out = 101),
               annotation_col = annot_df, 
               cluster_rows   = FALSE, 
               cluster_cols   = TRUE, 
               fontsize_row   = 8, 
               angle_col      = 45, 
               main           = "Age logFC (Old vs Young): top genes x cell type")
      dev.off()
    })
  }
}

# Plot 4: Expression Trajectory Line Plot
if (nrow(results_age) > 0) {
  top_genes_line <- results_age %>% arrange(FDR) %>% slice_head(n = 9) %>% pull(Gene)
  top_genes_line <- intersect(top_genes_line, rownames(logcpm_mat))
  
  if (length(top_genes_line) > 0) {
    line_data <- as.data.frame(t(logcpm_mat[top_genes_line, , drop = FALSE])) %>%
      rownames_to_column("sample_id") %>%
      mutate(AgeGroup = meta_pb$AgeGroup) %>%
      pivot_longer(-c(sample_id, AgeGroup), names_to = "Gene", values_to = "logCPM") %>%
      group_by(Gene, AgeGroup) %>%
      summarise(
        mean_expr = mean(logCPM, na.rm = TRUE),
        se_expr   = sd(logCPM, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      mutate(
        Gene = factor(Gene, levels = top_genes_line),
        AgeGroup = factor(AgeGroup, levels = c("Young (5-40)", "Middle-aged (41-65)", "Old (65+)"))
      )
    
    p_line <- ggplot(line_data, aes(x = AgeGroup, y = mean_expr, group = 1)) +
      geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr), width = 0.15, color = "grey50", linewidth = 0.8) +
      geom_line(color = "#2166ac", linewidth = 1) +
      geom_point(size = 3, color = "#d62728") +
      facet_wrap(~ Gene, scales = "free_y") +
      labs(title = "Average Gene Expression Across Age Bins",
           subtitle = "Top 9 significant genes (Mean log2-CPM ± Standard Error)",
           x = NULL, y = "Mean log2-CPM") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "grey95"),
            strip.text = element_text(face = "bold"))
    
    safe_save(p_line, "04_lineplot_expression_trajectory", 9, 8)
  }
}

message("\n✓  Complete.  All outputs in: ", OUTPUT_DIR)