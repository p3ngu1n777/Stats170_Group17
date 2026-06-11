# =============================================================================
#  Human Lung Cell Atlas  —  Negative-Binomial GLM (rhdf5 Direct Loading)
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────
library(rhdf5)
library(HDF5Array)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(edgeR)
library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(BiocParallel)

# ── 1. Load the Atlas directly with rhdf5 ─────────────────────────────────────
setwd("C:/Users/p3ngu/OneDrive/桌面/Stats170AB")
H5AD_PATH <- "fulldata.h5ad"

message("Loading h5ad with rhdf5...")

# ── 1a. Load count matrix (X) as a sparse HDF5-backed matrix ──────────────────
# X is stored as CSC sparse matrix: data, indices, indptr
message("  Reading count matrix...")

# Open the file connection first
h5_file <- H5Fopen(H5AD_PATH)

# Read using the open connection
x_data    <- h5read(H5AD_PATH, "/X/data")
x_indices <- h5read(H5AD_PATH, "/X/indices")
x_indptr  <- h5read(H5AD_PATH, "/X/indptr")

n_cells <- length(h5read(H5AD_PATH, "/obs/_index"))
n_genes <- length(h5read(H5AD_PATH, "/var/_index"))

# Close when done reading
H5Fclose(h5_file)

# Build sparse matrix (dgCMatrix) — CSC format
library(Matrix)
counts_sparse <- sparseMatrix(
  i = x_indices + 1L,       # 0-based -> 1-based
  p = x_indptr,
  x = x_data,
  dims = c(n_genes, n_cells)
)
rm(x_data, x_indices, x_indptr)
gc()
message("  Count matrix loaded: ", n_genes, " genes x ", n_cells, " cells.")

# ── 1b. Load obs metadata ──────────────────────────────────────────────────────
message("  Reading obs metadata...")

decode_col <- function(x) {
  # Categorical columns are stored as codes + categories (factor encoding)
  if (is.list(x) && all(c("categories", "codes") %in% names(x))) {
    cats  <- x$categories
    codes <- as.integer(x$codes) + 1L  # 0-based -> 1-based
    return(cats[codes])
  }
  return(x)
}

obs_raw     <- h5read(H5AD_PATH, "/obs")
obs_decoded <- lapply(obs_raw, decode_col)

# Keep only columns with exactly n_cells rows
obs_decoded <- obs_decoded[sapply(obs_decoded, function(x) !is.list(x) && length(x) == n_cells)]

# Extract cell IDs then remove _index
cell_ids <- obs_decoded[["_index"]]
obs_decoded[["_index"]] <- NULL

obs_df <- do.call(DataFrame, obs_decoded)
rownames(obs_df) <- cell_ids
message("  obs metadata loaded: ", ncol(obs_df), " columns.")

# ── 1c. Load var (gene) metadata ───────────────────────────────────────────────
message("  Reading var metadata...")
var_raw     <- h5read(H5AD_PATH, "/var")
var_decoded <- lapply(var_raw, decode_col)
var_decoded <- var_decoded[sapply(var_decoded, function(x) !is.list(x) && length(x) == n_genes)]

gene_ids <- var_decoded[["_index"]]
var_decoded[["_index"]] <- NULL

var_df <- do.call(DataFrame, var_decoded)
rownames(var_df) <- gene_ids

# ── 1d. Assemble SingleCellExperiment ─────────────────────────────────────────
rownames(counts_sparse) <- gene_ids
colnames(counts_sparse) <- cell_ids

sce <- SingleCellExperiment(
  assays  = list(X = counts_sparse),
  rowData = var_df,
  colData = obs_df
)
rm(counts_sparse, obs_df, var_df, obs_raw, var_raw)
gc()
message("SCE assembled: ", ncol(sce), " cells, ", nrow(sce), " genes.")




# ── 2. Column Name Constants ───────────────────────────────────────────────────
DONOR_COL    <- "donor_id"
CELLTYPE_COL <- "ann_level_3"
AGE_COL      <- "age_or_mean_of_age_range"
SEX_COL      <- "sex"
SMOKING_COL  <- "smoking_status"
DATASET_COL  <- "dataset"

# ── 3. Derive AgeGroup ────────────────────────────────────────────────────────
message("Processing age metadata...")

if (!AGE_COL %in% names(colData(sce))) {
  stop(paste("Column", AGE_COL, "not found. Available:", 
             paste(names(colData(sce)), collapse = ", ")))
}

# Already numeric float — no string parsing needed
numeric_age <- as.numeric(sce[[AGE_COL]])

sce$AgeGroup <- cut(
  numeric_age,
  breaks = c(0, 29, 49, 64, Inf),
  labels = c("Young (<30)", "Adult (30-49)", "Middle (50-64)", "Older (65+)"),
  right  = TRUE
)

sce <- sce[, !is.na(sce$AgeGroup)]
message("AgeGroup derived: ", ncol(sce), " cells retained.")
print(table(sce$AgeGroup))

# ── 4. Pseudo-Bulk Aggregation ────────────────────────────────────────────────
message("Aggregating counts to pseudo-bulk (this may take a few minutes)...")

pb <- aggregateAcrossCells(
  sce,
  ids      = colData(sce)[, c(DONOR_COL, CELLTYPE_COL)],
  use.assay.type = "X"
)

pb$lib_size <- colSums(assay(pb, "sum"))
rm(sce)
gc()

# ── 5. Filtering ──────────────────────────────────────────────────────────────
keep_genes <- rowSums(assay(pb, "sum") >= 10) >= 5
pb_filt    <- pb[keep_genes, ]
pb_filt    <- pb_filt[, pb_filt$lib_size > 0]
message(nrow(pb_filt), " genes retained after filtering.")

# ── 6. Build DGEList & Design Matrix ──────────────────────────────────────────
counts_mat <- as.matrix(assay(pb_filt, "sum"))
meta_pb    <- as.data.frame(colData(pb_filt))

meta_pb <- meta_pb %>%
  mutate(
    AgeGroup = relevel(droplevels(factor(AgeGroup)), ref = "Adult (30-49)"),
    CellType = droplevels(factor(.data[[CELLTYPE_COL]])),
    Sex      = factor(.data[[SEX_COL]]),
    Smoking  = factor(.data[[SMOKING_COL]]),
    Dataset  = factor(.data[[DATASET_COL]])
  )

design <- model.matrix(~ AgeGroup * CellType + Sex + Smoking + Dataset, data = meta_pb)
dge    <- DGEList(counts = counts_mat, lib.size = meta_pb$lib_size)
dge    <- calcNormFactors(dge, method = "TMM")

# ── 7. Estimate Dispersion & Fit NB-GLM ───────────────────────────────────────
message("Estimating NB dispersions and fitting QL-GLM...")
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)

# ── 8. Hypothesis Tests ───────────────────────────────────────────────────────
run_test <- function(fit, coef_pattern, label) {
  target_coefs <- grep(coef_pattern, colnames(design), value = TRUE, perl = TRUE)
  if (length(target_coefs) == 0) return(NULL)
  res <- glmQLFTest(fit, coef = target_coefs)
  topTags(res, n = Inf)$table %>%
    as_tibble(rownames = "Gene") %>%
    mutate(Test = label)
}

results_age   <- run_test(fit, "^AgeGroup(?!.*:)", "AgeGroup_main")
results_int   <- run_test(fit, "AgeGroup.*:CellType|CellType.*:AgeGroup", "Interaction")
results_sex   <- run_test(fit, "^Sex", "Sex_main")
results_smoke <- run_test(fit, "^Smoking", "Smoking_main")

all_results <- bind_rows(results_age, results_int, results_sex, results_smoke)

# ── 9. Save Results ───────────────────────────────────────────────────────────
output_dir <- file.path(dirname(H5AD_PATH), "HLCA_NB_GLM_results")
dir.create(output_dir, showWarnings = FALSE)
write.csv(all_results, file.path(output_dir, "all_results.csv"), row.names = FALSE)
saveRDS(fit, file.path(output_dir, "glmQLFit.rds"))

# ── 10. Diagnostic & Summary Plots ────────────────────────────────────────────
pdf(file.path(output_dir, "HLCA_NB_GLM_diagnostics.pdf"), width = 12, height = 8)
plotBCV(dge)
plotQLDisp(fit)
if (!is.null(results_age)) {
  print(ggplot(results_age, aes(logFC, -log10(PValue), color = FDR < 0.05)) +
          geom_point(alpha = 0.5) + theme_minimal() + labs(title = "Age Main Effect"))
}
dev.off()

# ── 11. Cell-Type-Specific Age Effects ────────────────────────────────────────
message("Running stratified analysis per cell type...")

run_per_celltype <- function(ct_label, pb_se) {
  sub_se <- pb_se[, pb_se[[CELLTYPE_COL]] == ct_label]
  if (ncol(sub_se) < 8) return(NULL)
  sub_meta <- as.data.frame(colData(sub_se))
  if (length(unique(sub_meta$AgeGroup)) < 2) return(NULL)
  dm <- tryCatch(model.matrix(~ AgeGroup + Sex + Smoking + Dataset, data = sub_meta),
                 error = function(e) NULL)
  if (is.null(dm)) return(NULL)
  dge_ct <- DGEList(counts = assay(sub_se, "sum"), lib.size = sub_meta$lib_size)
  dge_ct <- calcNormFactors(dge_ct)
  dge_ct <- estimateDisp(dge_ct, dm)
  fit_ct <- glmQLFit(dge_ct, dm)
  age_c  <- grep("^AgeGroup", colnames(dm), value = TRUE)
  if (length(age_c) == 0) return(NULL)
  topTags(glmQLFTest(fit_ct, coef = age_c), n = Inf)$table %>%
    as_tibble(rownames = "Gene") %>% mutate(CellType = ct_label)
}

ct_results <- lapply(levels(meta_pb$CellType), run_per_celltype, pb_se = pb_filt)
write.csv(bind_rows(ct_results), file.path(output_dir, "per_celltype_age_results.csv"))

# ── 12. Session Info ──────────────────────────────────────────────────────────
sink(file.path(output_dir, "session_info.txt"))
print(sessionInfo())
sink()
message("Analysis Complete.")