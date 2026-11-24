# ============================================================================
# PAX6 CERVICAL CANCER: FIVE-GENE CO-EXPRESSION ANALYSIS
#
# Purpose
#   Reproduce the five-gene Spearman co-expression analysis, Figure 5,
#   Figure 6, Table 4, and the associated sensitivity analyses.
#
# Required input files (kept locally; not committed to the public repository)
#   data/TCGA_304_FiveGene_log2TPM.csv
#   data/TCGA_304_Aligned_Manifest.csv
#
# Expression input is ALREADY on the log2(TPM + 1) scale. This script does
# NOT transform it again.
#
# Run from the repository root:
#   source("PAX6_Coexpression_Analysis.R")
# or:
#   Rscript PAX6_Coexpression_Analysis.R
#
# Outputs are written to Results/, Figures/, and Tables/.
# ============================================================================

# ------------------------------ USER SETTINGS -------------------------------
EXPRESSION_FILE <- Sys.getenv(
  "PAX6_EXPRESSION_FILE",
  unset = file.path("data", "TCGA_304_FiveGene_log2TPM.csv")
)
MANIFEST_FILE <- Sys.getenv(
  "PAX6_MANIFEST_FILE",
  unset = file.path("data", "TCGA_304_Aligned_Manifest.csv")
)

EXPECTED_N <- 304L
BOOTSTRAP_REPLICATES <- 5000L
BOOTSTRAP_SEED <- 20261124L
SIGNIFICANCE_LEVEL <- 0.05
GENE_ORDER <- c("PAX6", "SOX1", "SOX2", "NEUROG2", "SIX3")
EXPRESSION_SCALE <- "log2(TPM + 1)"

# Previously reported Spearman values and BH q-values are used only for
# quality-control comparison. They are never substituted for calculated data.
EXPECTED_RESULTS <- data.frame(
  Gene_pair = c(
    "PAX6–SOX1", "PAX6–SOX2", "PAX6–NEUROG2", "PAX6–SIX3",
    "SOX1–SOX2", "SOX1–NEUROG2", "SOX1–SIX3",
    "SOX2–NEUROG2", "SOX2–SIX3", "NEUROG2–SIX3"
  ),
  Expected_rho = c(
    0.011, 0.296, 0.055, 0.274, 0.038,
    0.054, 0.081, 0.125, 0.160, 0.046
  ),
  Expected_q = c(
    0.845, 1.41e-6, 0.499, 6.41e-6, 0.561,
    0.499, 0.316, 0.074, 0.0172, 0.530
  ),
  stringsAsFactors = FALSE
)

# ------------------------------- DEPENDENCIES --------------------------------
required_packages <- c("ggplot2", "openxlsx")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))\nThen rerun this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
})

# ------------------------------- INPUT CHECKS --------------------------------
if (!file.exists(EXPRESSION_FILE)) {
  stop(
    "Expression matrix not found: ", EXPRESSION_FILE,
    "\nPlace TCGA_304_FiveGene_log2TPM.csv in the local data/ folder, or set ",
    "PAX6_EXPRESSION_FILE to its path.",
    call. = FALSE
  )
}
if (!file.exists(MANIFEST_FILE)) {
  stop(
    "Aligned manifest not found: ", MANIFEST_FILE,
    "\nPlace TCGA_304_Aligned_Manifest.csv in the local data/ folder, or set ",
    "PAX6_MANIFEST_FILE to its path.",
    call. = FALSE
  )
}

expression_data <- tryCatch(
  read.csv(EXPRESSION_FILE, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop("Could not read expression matrix: ", conditionMessage(e), call. = FALSE)
)
manifest <- tryCatch(
  read.csv(MANIFEST_FILE, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop("Could not read manifest: ", conditionMessage(e), call. = FALSE)
)

required_expression_columns <- c("sample_id", GENE_ORDER)
missing_expression_columns <- setdiff(required_expression_columns, names(expression_data))
if (length(missing_expression_columns) > 0L) {
  stop(
    "Expression matrix is missing column(s): ",
    paste(missing_expression_columns, collapse = ", "),
    call. = FALSE
  )
}
required_manifest_columns <- c(
  "sample_id", "patient_id", "sample_type", "workflow", "workflow_version"
)
missing_manifest_columns <- setdiff(required_manifest_columns, names(manifest))
if (length(missing_manifest_columns) > 0L) {
  stop(
    "Manifest is missing column(s): ",
    paste(missing_manifest_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(expression_data) != EXPECTED_N) {
  stop("Expected ", EXPECTED_N, " expression rows; found ", nrow(expression_data), ".", call. = FALSE)
}
if (nrow(manifest) != EXPECTED_N) {
  stop("Expected ", EXPECTED_N, " manifest rows; found ", nrow(manifest), ".", call. = FALSE)
}

for (object_name in c("expression_data", "manifest")) {
  object <- get(object_name)
  ids <- object$sample_id
  if (anyNA(ids) || any(trimws(as.character(ids)) == "")) {
    stop(object_name, "$sample_id contains missing or blank values.", call. = FALSE)
  }
  if (anyDuplicated(as.character(ids)) > 0L) {
    stop(object_name, " contains duplicate sample_id values.", call. = FALSE)
  }
}

# Check that the expression matrix and metadata manifest refer to exactly the
# same 304 samples, then align the manifest to the expression-matrix order.
expression_ids <- as.character(expression_data$sample_id)
manifest_ids <- as.character(manifest$sample_id)
if (!setequal(expression_ids, manifest_ids)) {
  stop(
    "Expression matrix and manifest do not contain the same sample_id set. ",
    "Do not proceed until the sample lists are reconciled.",
    call. = FALSE
  )
}
manifest <- manifest[match(expression_ids, manifest_ids), , drop = FALSE]
if (!identical(expression_ids, as.character(manifest$sample_id))) {
  stop("Manifest alignment failed.", call. = FALSE)
}

if (anyNA(manifest$patient_id) || any(trimws(as.character(manifest$patient_id)) == "")) {
  stop("Manifest contains missing or blank patient_id values.", call. = FALSE)
}
if (length(unique(as.character(manifest$patient_id))) != EXPECTED_N) {
  stop(
    "Manifest does not contain one unique patient_id per expression sample. ",
    "Review sample selection before analysis.",
    call. = FALSE
  )
}
if (!all(trimws(as.character(manifest$sample_type)) == "Primary Tumor")) {
  stop("Manifest includes records not labelled 'Primary Tumor'.", call. = FALSE)
}
if (length(unique(as.character(manifest$workflow))) != 1L ||
    !all(as.character(manifest$workflow) == "STAR - Counts")) {
  stop("Manifest does not consistently identify the STAR - Counts workflow.", call. = FALSE)
}
if (length(unique(as.character(manifest$workflow_version))) != 1L) {
  stop("Manifest contains more than one workflow_version.", call. = FALSE)
}

# Expression file name and supplied project context identify the values as
# already transformed log2(TPM + 1). Do not apply the transformation again.
for (gene in GENE_ORDER) {
  values <- suppressWarnings(as.numeric(as.character(expression_data[[gene]])))
  if (anyNA(values) || any(!is.finite(values))) {
    stop("Gene column ", gene, " contains missing or non-numeric values.", call. = FALSE)
  }
  if (any(values < 0)) {
    stop("Gene column ", gene, " contains negative values; expected log2(TPM + 1).", call. = FALSE)
  }
  expression_data[[gene]] <- values
}

if (any(vapply(expression_data[GENE_ORDER], function(x) length(unique(x)) < 2L, logical(1)))) {
  stop("At least one gene has no expression variation; correlation is undefined.", call. = FALSE)
}

# ------------------------------- OUTPUT FOLDERS ------------------------------
dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Tables", showWarnings = FALSE, recursive = TRUE)

# ---------------------------- HELPER FUNCTIONS -------------------------------
safe_correlation_test <- function(x, y, method = "spearman") {
  result <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  list(
    estimate = unname(as.numeric(result$estimate)),
    p_value = unname(as.numeric(result$p.value))
  )
}

bootstrap_spearman_ci <- function(x, y, B, seed) {
  set.seed(seed)
  n <- length(x)
  bootstrap_values <- replicate(B, {
    indices <- sample.int(n, size = n, replace = TRUE)
    estimate <- suppressWarnings(cor(x[indices], y[indices], method = "spearman"))
    if (is.finite(estimate)) as.numeric(estimate) else NA_real_
  })
  bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]
  if (length(bootstrap_values) < ceiling(0.95 * B)) {
    stop(
      "Fewer than 95% of bootstrap replicates were valid. Inspect the expression data.",
      call. = FALSE
    )
  }
  interval <- unname(stats::quantile(
    bootstrap_values, probs = c(0.025, 0.975), na.rm = TRUE,
    names = FALSE, type = 7
  ))
  list(
    lower = interval[1],
    upper = interval[2],
    valid_replicates = length(bootstrap_values)
  )
}

format_p_value <- function(x) {
  vapply(x, function(value) {
    if (is.na(value)) return(NA_character_)
    if (value < 0.001) formatC(value, format = "e", digits = 2)
    else sprintf("%.3f", value)
  }, character(1))
}

# -------------------------- PRIMARY SPEARMAN ANALYSIS ------------------------
pair_list <- combn(GENE_ORDER, 2L, simplify = FALSE)
spearman_rows <- vector("list", length(pair_list))

for (i in seq_along(pair_list)) {
  gene1 <- pair_list[[i]][1]
  gene2 <- pair_list[[i]][2]
  test_result <- safe_correlation_test(
    expression_data[[gene1]], expression_data[[gene2]], method = "spearman"
  )
  bootstrap_result <- bootstrap_spearman_ci(
    expression_data[[gene1]], expression_data[[gene2]],
    B = BOOTSTRAP_REPLICATES,
    seed = BOOTSTRAP_SEED + i
  )

  spearman_rows[[i]] <- data.frame(
    Gene1 = gene1,
    Gene2 = gene2,
    Gene_pair = paste(gene1, gene2, sep = "–"),
    N = EXPECTED_N,
    Spearman_rho = test_result$estimate,
    CI_lower = bootstrap_result$lower,
    CI_upper = bootstrap_result$upper,
    Raw_p_value = test_result$p_value,
    Bootstrap_requested = BOOTSTRAP_REPLICATES,
    Bootstrap_valid = bootstrap_result$valid_replicates,
    stringsAsFactors = FALSE
  )
}

spearman_results <- do.call(rbind, spearman_rows)
spearman_results$BH_adjusted_q <- p.adjust(spearman_results$Raw_p_value, method = "BH")
spearman_results$Significant_BH_005 <- spearman_results$BH_adjusted_q < SIGNIFICANCE_LEVEL

# ----------------------------- KENDALL TAU-B --------------------------------
kendall_rows <- vector("list", length(pair_list))
for (i in seq_along(pair_list)) {
  gene1 <- pair_list[[i]][1]
  gene2 <- pair_list[[i]][2]
  test_result <- safe_correlation_test(
    expression_data[[gene1]], expression_data[[gene2]], method = "kendall"
  )
  kendall_rows[[i]] <- data.frame(
    Gene1 = gene1,
    Gene2 = gene2,
    Gene_pair = paste(gene1, gene2, sep = "–"),
    N = EXPECTED_N,
    Kendall_tau_b = test_result$estimate,
    Raw_p_value = test_result$p_value,
    stringsAsFactors = FALSE
  )
}
kendall_results <- do.call(rbind, kendall_rows)
kendall_results$BH_adjusted_q <- p.adjust(kendall_results$Raw_p_value, method = "BH")
kendall_results$Significant_BH_005 <- kendall_results$BH_adjusted_q < SIGNIFICANCE_LEVEL

# ---------------------------- LEAVE-ONE-OUT ----------------------------------
# LOO results intentionally use sample indices only and do not export IDs.
loo_pairs <- list(c("PAX6", "SOX2"), c("PAX6", "SIX3"))
loo_detail_rows <- list()
loo_summary_rows <- list()

for (i in seq_along(loo_pairs)) {
  gene1 <- loo_pairs[[i]][1]
  gene2 <- loo_pairs[[i]][2]
  full_rho <- safe_correlation_test(
    expression_data[[gene1]], expression_data[[gene2]], method = "spearman"
  )$estimate
  loo_rho <- vapply(seq_len(EXPECTED_N), function(omitted_index) {
    retained <- setdiff(seq_len(EXPECTED_N), omitted_index)
    suppressWarnings(cor(
      expression_data[[gene1]][retained],
      expression_data[[gene2]][retained],
      method = "spearman"
    ))
  }, numeric(1))

  if (any(!is.finite(loo_rho))) {
    stop("Non-finite leave-one-out estimate for ", gene1, "–", gene2, ".", call. = FALSE)
  }

  loo_detail_rows[[i]] <- data.frame(
    Gene_pair = paste(gene1, gene2, sep = "–"),
    Omitted_sample_index = seq_len(EXPECTED_N),
    Spearman_rho_after_omission = loo_rho,
    stringsAsFactors = FALSE
  )
  loo_summary_rows[[i]] <- data.frame(
    Gene_pair = paste(gene1, gene2, sep = "–"),
    Full_data_Spearman_rho = full_rho,
    LOO_minimum_rho = min(loo_rho),
    LOO_maximum_rho = max(loo_rho),
    All_LOO_estimates_positive = all(loo_rho > 0),
    Number_of_omissions = length(loo_rho),
    stringsAsFactors = FALSE
  )
}
loo_details <- do.call(rbind, loo_detail_rows)
loo_summary <- do.call(rbind, loo_summary_rows)

# -------------------------- BENCHMARK VALIDATION -----------------------------
benchmark <- EXPECTED_RESULTS[match(spearman_results$Gene_pair, EXPECTED_RESULTS$Gene_pair), ]
rho_difference <- abs(spearman_results$Spearman_rho - benchmark$Expected_rho)
q_difference <- abs(spearman_results$BH_adjusted_q - benchmark$Expected_q)
rho_matches <- rho_difference <= 0.0011
q_matches <- q_difference <= pmax(1e-8, benchmark$Expected_q * 0.02)

if (!all(rho_matches) || !all(q_matches)) {
  comparison <- data.frame(
    Gene_pair = spearman_results$Gene_pair,
    Calculated_rho = spearman_results$Spearman_rho,
    Reference_rho_rounded = benchmark$Expected_rho,
    Calculated_q = spearman_results$BH_adjusted_q,
    Reference_q_rounded = benchmark$Expected_q,
    Rho_within_tolerance = rho_matches,
    q_within_tolerance = q_matches
  )
  print(comparison, row.names = FALSE)
  stop(
    "Calculated Spearman results do not match the documented benchmark. ",
    "Do not use manuscript outputs until sample selection and data provenance are reviewed.",
    call. = FALSE
  )
}

expected_significant <- c("PAX6–SOX2", "PAX6–SIX3", "SOX2–SIX3")
actual_significant <- spearman_results$Gene_pair[spearman_results$Significant_BH_005]
if (!setequal(actual_significant, expected_significant)) {
  stop("The FDR-significant pair set differs from the documented three pairs.", call. = FALSE)
}

# ------------------------------- FIGURE 5 ------------------------------------
pax6_partners <- GENE_ORDER[-1]
scatter_data <- do.call(rbind, lapply(pax6_partners, function(partner) {
  data.frame(
    PAX6 = expression_data$PAX6,
    Partner_expression = expression_data[[partner]],
    Partner = partner,
    stringsAsFactors = FALSE
  )
}))
scatter_data$Partner <- factor(scatter_data$Partner, levels = pax6_partners)

pax6_results <- spearman_results[spearman_results$Gene1 == "PAX6", , drop = FALSE]
annotation_data <- data.frame(
  Partner = factor(pax6_results$Gene2, levels = pax6_partners),
  label = paste0(
    "Spearman rho = ", sprintf("%.3f", pax6_results$Spearman_rho),
    "\n95% bootstrap CI [", sprintf("%.3f", pax6_results$CI_lower), ", ",
    sprintf("%.3f", pax6_results$CI_upper), "]",
    "\nBH q = ", format_p_value(pax6_results$BH_adjusted_q)
  ),
  stringsAsFactors = FALSE
)

figure5 <- ggplot(scatter_data, aes(x = PAX6, y = Partner_expression)) +
  geom_point(alpha = 0.65, size = 1.6) +
  facet_wrap(~ Partner, scales = "free_y", ncol = 2) +
  geom_label(
    data = annotation_data,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.10,
    size = 2.8, label.size = 0.2
  ) +
  labs(
    title = "PAX6-centred co-expression with selected developmental regulators",
    subtitle = "TCGA-CESC primary tumour samples (N = 304)",
    x = "PAX6 expression [log2(TPM + 1)]",
    y = "Partner gene expression [log2(TPM + 1)]",
    caption = paste(
      "Two-sided Spearman correlations; BH correction across ten unique gene pairs.",
      "No linear regression lines are displayed."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  filename = file.path("Figures", "Figure5_PAX6_Coexpression_Scatterplots.tiff"),
  plot = figure5, device = "tiff", width = 10, height = 8,
  units = "in", dpi = 600, compression = "lzw", bg = "white"
)
ggsave(
  filename = file.path("Figures", "Figure5_PAX6_Coexpression_Scatterplots.png"),
  plot = figure5, device = "png", width = 10, height = 8,
  units = "in", dpi = 600, bg = "white"
)

# ------------------------------- FIGURE 6 ------------------------------------
rho_matrix <- diag(length(GENE_ORDER))
rownames(rho_matrix) <- GENE_ORDER
colnames(rho_matrix) <- GENE_ORDER
q_matrix <- matrix(
  NA_real_, nrow = length(GENE_ORDER), ncol = length(GENE_ORDER),
  dimnames = list(GENE_ORDER, GENE_ORDER)
)

for (i in seq_len(nrow(spearman_results))) {
  gene1 <- spearman_results$Gene1[i]
  gene2 <- spearman_results$Gene2[i]
  rho_matrix[gene1, gene2] <- spearman_results$Spearman_rho[i]
  rho_matrix[gene2, gene1] <- spearman_results$Spearman_rho[i]
  q_matrix[gene1, gene2] <- spearman_results$BH_adjusted_q[i]
  q_matrix[gene2, gene1] <- spearman_results$BH_adjusted_q[i]
}

grid_data <- expand.grid(
  Gene1 = GENE_ORDER,
  Gene2 = rev(GENE_ORDER),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid_data$rho <- mapply(
  function(row_gene, column_gene) rho_matrix[row_gene, column_gene],
  grid_data$Gene2, grid_data$Gene1
)
grid_data$q <- mapply(
  function(row_gene, column_gene) q_matrix[row_gene, column_gene],
  grid_data$Gene2, grid_data$Gene1
)
grid_data$label <- sprintf("%.2f", grid_data$rho)
significant_cells <- !is.na(grid_data$q) & grid_data$q < SIGNIFICANCE_LEVEL
grid_data$label[significant_cells] <- paste0(grid_data$label[significant_cells], "*")
grid_data$Gene1 <- factor(grid_data$Gene1, levels = GENE_ORDER)
grid_data$Gene2 <- factor(grid_data$Gene2, levels = rev(GENE_ORDER))

if (nrow(grid_data) != 25L || sum(is.na(grid_data$rho)) != 0L) {
  stop("Figure 6 matrix construction failed validation.", call. = FALSE)
}

figure6 <- ggplot(grid_data, aes(x = Gene1, y = Gene2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = label), color = "black", size = 4) +
  scale_x_discrete(limits = GENE_ORDER, drop = FALSE) +
  scale_y_discrete(limits = rev(GENE_ORDER), drop = FALSE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1), name = "Spearman rho"
  ) +
  coord_fixed() +
  labs(
    title = "Pairwise co-expression among PAX6 and selected developmental regulators",
    subtitle = "TCGA-CESC primary tumour samples (N = 304)",
    x = NULL, y = NULL,
    caption = paste(
      "* BH-adjusted q < 0.05 across ten unique off-diagonal comparisons.",
      "Mirrored cells represent the same association; diagonal values equal 1.00."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  filename = file.path("Figures", "Figure6_FiveGene_Correlation_Heatmap.tiff"),
  plot = figure6, device = "tiff", width = 8.5, height = 7.5,
  units = "in", dpi = 600, compression = "lzw", bg = "white"
)
ggsave(
  filename = file.path("Figures", "Figure6_FiveGene_Correlation_Heatmap.png"),
  plot = figure6, device = "png", width = 8.5, height = 7.5,
  units = "in", dpi = 600, bg = "white"
)

# -------------------------------- TABLE 4 ------------------------------------
manuscript_table <- data.frame(
  `Gene pair` = spearman_results$Gene_pair,
  N = spearman_results$N,
  `Spearman rho` = paste0(
    sprintf("%.3f", spearman_results$Spearman_rho),
    ifelse(spearman_results$Significant_BH_005, "*", "")
  ),
  `95% bootstrap CI` = paste0(
    "[", sprintf("%.3f", spearman_results$CI_lower), ", ",
    sprintf("%.3f", spearman_results$CI_upper), "]"
  ),
  `Raw p-value` = format_p_value(spearman_results$Raw_p_value),
  `BH-adjusted q-value` = format_p_value(spearman_results$BH_adjusted_q),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Export numerical results and sensitivity analyses.
write.csv(
  spearman_results, file.path("Results", "Spearman_All10_FullPrecision.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  kendall_results, file.path("Results", "Kendall_TauB_Sensitivity.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  loo_summary, file.path("Results", "LeaveOneOut_Sensitivity_Summary.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  loo_details, file.path("Results", "LeaveOneOut_Sensitivity_Details.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  manuscript_table, file.path("Tables", "Table4_Coexpression_Manuscript.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  spearman_results, file.path("Tables", "Table4_Coexpression_FullPrecision.csv"),
  row.names = FALSE, na = ""
)

# Input QC summary contains aggregate metadata only; no sample IDs are exported.
qc_summary <- data.frame(
  Metric = c(
    "Expression samples", "Unique patients", "Manifest sample type",
    "Workflow", "Workflow version", "Expression scale",
    "Expression matrix rows aligned to manifest", "Unique gene-pair tests",
    "BH-significant Spearman associations", "Bootstrap resamples per pair"
  ),
  Value = c(
    nrow(expression_data), length(unique(manifest$patient_id)),
    paste(unique(as.character(manifest$sample_type)), collapse = "; "),
    paste(unique(as.character(manifest$workflow)), collapse = "; "),
    paste(unique(as.character(manifest$workflow_version)), collapse = "; "),
    EXPRESSION_SCALE, "Yes", nrow(spearman_results),
    paste(actual_significant, collapse = "; "), BOOTSTRAP_REPLICATES
  ),
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path("Results", "Input_QC_Summary.csv"), row.names = FALSE)

# ---------------------- TIMES NEW ROMAN EXCEL WORKBOOK -----------------------
workbook <- createWorkbook()
sheets <- c("Table 4", "Full precision", "Kendall tau-b", "LOO summary", "LOO details", "Input QC", "Notes")
for (sheet in sheets) addWorksheet(workbook, sheet)

caption <- paste(
  "Table 4. Pairwise transcriptomic co-expression among PAX6 and selected",
  "developmental regulators in TCGA-CESC primary tumour samples (N = 304)."
)
writeData(workbook, "Table 4", caption, startRow = 1, startCol = 1, colNames = FALSE)
writeData(workbook, "Table 4", manuscript_table, startRow = 3, startCol = 1)
footnote_row <- nrow(manuscript_table) + 5L
writeData(
  workbook, "Table 4",
  paste(
    "* BH-adjusted q < 0.05 across ten unique pairwise comparisons.",
    "Confidence intervals are pointwise percentile-bootstrap intervals based on 5,000 resamples."
  ),
  startRow = footnote_row, startCol = 1, colNames = FALSE
)
writeData(workbook, "Full precision", spearman_results)
writeData(workbook, "Kendall tau-b", kendall_results)
writeData(workbook, "LOO summary", loo_summary)
writeData(workbook, "LOO details", loo_details)
writeData(workbook, "Input QC", qc_summary)

notes <- data.frame(
  Item = c(
    "Cohort", "Expression input", "Expression scale", "Primary correlation",
    "Multiple testing", "Bootstrap confidence intervals", "Kendall sensitivity",
    "Leave-one-out sensitivity", "Significance threshold", "Interpretation"
  ),
  Description = c(
    "304 TCGA-CESC primary tumour samples, one sample per patient.",
    basename(EXPRESSION_FILE),
    "Input values are already log2(TPM + 1); no additional transformation is applied.",
    "Two-sided Spearman rank correlation for ten unique gene pairs.",
    "Benjamini-Hochberg correction across all ten unique comparisons; mirrored heatmap cells are not separate tests.",
    paste0("5,000 non-parametric bootstrap resamples per pair; base seed = ", BOOTSTRAP_SEED,
           "; pointwise percentile 95% intervals.",
           " Bootstrap intervals are reproducible for this script version and R RNG implementation."),
    "Two-sided Kendall tau-b correlations; BH correction across ten unique pairs.",
    "Each of 304 samples omitted individually for PAX6-SOX2 and PAX6-SIX3; sample identifiers are not exported.",
    "BH-adjusted q < 0.05.",
    "Co-expression associations do not establish direct regulation, physical interaction or causality."
  ),
  stringsAsFactors = FALSE
)
writeData(workbook, "Notes", notes)

header_style <- createStyle(
  fontName = "Times New Roman", fontSize = 9, textDecoration = "bold",
  fontColour = "#FFFFFF", fgFill = "#34495E",
  halign = "center", valign = "center", wrapText = TRUE
)
body_style <- createStyle(
  fontName = "Times New Roman", fontSize = 9,
  valign = "center", wrapText = TRUE
)
title_style <- createStyle(
  fontName = "Times New Roman", fontSize = 10, textDecoration = "bold",
  wrapText = TRUE
)
note_style <- createStyle(
  fontName = "Times New Roman", fontSize = 9, textDecoration = "italic",
  wrapText = TRUE
)

sheet_data <- list(
  "Table 4" = list(data = manuscript_table, header_row = 3L),
  "Full precision" = list(data = spearman_results, header_row = 1L),
  "Kendall tau-b" = list(data = kendall_results, header_row = 1L),
  "LOO summary" = list(data = loo_summary, header_row = 1L),
  "LOO details" = list(data = loo_details, header_row = 1L),
  "Input QC" = list(data = qc_summary, header_row = 1L),
  "Notes" = list(data = notes, header_row = 1L)
)

for (sheet in names(sheet_data)) {
  dat <- sheet_data[[sheet]]$data
  header_row <- sheet_data[[sheet]]$header_row
  addStyle(
    workbook, sheet, header_style,
    rows = header_row, cols = seq_len(ncol(dat)),
    gridExpand = TRUE, stack = TRUE
  )
  if (nrow(dat) > 0L) {
    addStyle(
      workbook, sheet, body_style,
      rows = (header_row + 1L):(header_row + nrow(dat)),
      cols = seq_len(ncol(dat)), gridExpand = TRUE, stack = TRUE
    )
  }
  setColWidths(workbook, sheet, cols = seq_len(ncol(dat)), widths = "auto")
  freezePane(workbook, sheet, firstActiveRow = header_row + 1L)
}
addStyle(workbook, "Table 4", title_style, rows = 1, cols = 1, stack = TRUE)
addStyle(workbook, "Table 4", note_style, rows = footnote_row, cols = 1, stack = TRUE)
mergeCells(workbook, "Table 4", cols = 1:ncol(manuscript_table), rows = 1)
mergeCells(workbook, "Table 4", cols = 1:ncol(manuscript_table), rows = footnote_row)

saveWorkbook(workbook, file.path("Tables", "Table4_Coexpression.xlsx"), overwrite = TRUE)

# ------------------------------ SESSION INFO ---------------------------------
writeLines(capture.output(sessionInfo()), file.path("Results", "R_sessionInfo.txt"))

# ---------------------------------- SUMMARY ----------------------------------
cat("\n========================================\n")
cat("PAX6 CO-EXPRESSION ANALYSIS COMPLETED\n")
cat("========================================\n")
cat("Expression input:", EXPRESSION_FILE, "\n")
cat("Manifest input:", MANIFEST_FILE, "\n")
cat("Samples:", EXPECTED_N, "\n")
cat("Expression scale:", EXPRESSION_SCALE, "(used as supplied; no re-transformation)\n")
cat("Manifest: 304 Primary Tumor samples; STAR - Counts workflow.\n")
cat("Unique gene-pair tests:", nrow(spearman_results), "\n")
cat("BH family: ten unique comparisons\n")
cat("Bootstrap resamples per pair:", BOOTSTRAP_REPLICATES, "\n")
cat("Significant Spearman pairs:", paste(actual_significant, collapse = ", "), "\n")
cat("\nSpearman results:\n")
print(
  spearman_results[, c("Gene_pair", "N", "Spearman_rho", "CI_lower", "CI_upper", "Raw_p_value", "BH_adjusted_q")],
  row.names = FALSE
)
cat("\nAll benchmark checks passed.\n")
cat("Outputs saved in Results/, Figures/, and Tables/.\n")
cat("========================================\n")
