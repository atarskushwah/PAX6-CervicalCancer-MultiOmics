# ==============================================================================
# Correlation Analysis of PAX6 and Candidate Genes (SOX1, SOX2, NEUROG2, SIX3)
# Dataset: TCGA-CESC (Cervical Squamous Cell Carcinoma & Endocervical Adenocarcinoma)
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0: Install and Load Required Packages
# ------------------------------------------------------------------------------
required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "tidyverse", 
                       "patchwork", "ggcorrplot", "car", "broom")

# Install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "TCGAbiolinks" || pkg == "SummarizedExperiment") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(patchwork)
library(ggcorrplot)
library(car)
library(broom)

# ------------------------------------------------------------------------------
# STEP 1: Download/Import the Cervical Cancer Expression Dataset
# ------------------------------------------------------------------------------
# Note: If automatic download fails, you can load a pre-downloaded matrix here:
# expr_data <- read.csv("your_expression_matrix.csv") 

query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
tcga_data <- GDCprepare(query)

# ------------------------------------------------------------------------------
# STEP 2 & 3: Clean the Dataset and Extract Gene Expression Values
# ------------------------------------------------------------------------------
# Filter for primary tumor samples only (exclude normal/control samples)
tumor_samples <- tcga_data[, tcga_data$sample_type == "Primary Tumor"]

# Extract clinical metadata for later stratification
# Using any_of() prevents the script from crashing if standard stage variables are missing/renamed
clinical_meta <- colData(tumor_samples) %>%
  as.data.frame() %>%
  select(
    any_of(c(
      "barcode", 
      "primary_diagnosis", 
      "tumor_stage", 
      "ajcc_pathologic_stage", 
      "figo_stage"
    ))
  )

# Extract TPM (Transcripts Per Million) unstranded values
# Log2 transformation: log2(TPM + 1)
expr_matrix <- assay(tumor_samples, "tpm_unstrand")
expr_matrix_log2 <- log2(expr_matrix + 1)

# Get gene mapping to HGNC symbols
gene_info <- rowData(tumor_samples) %>% as.data.frame()
gene_mapping <- data.frame(
  ensembl_id = rownames(expr_matrix_log2),
  gene_name = gene_info$gene_name
)

# Target genes
target_genes <- c("PAX6", "SOX1", "SOX2", "NEUROG2", "SIX3")

# Identify correct Ensembl IDs for target genes
target_ids <- gene_mapping %>% 
  filter(gene_name %in% target_genes) %>%
  # In case of multiple Ensembl IDs mapping to one symbol, take the first to avoid duplicates
  distinct(gene_name, .keep_all = TRUE) 

# Subset the expression matrix
sub_expr <- expr_matrix_log2[target_ids$ensembl_id, ]
rownames(sub_expr) <- target_ids$gene_name

# Transpose to have samples as rows and genes as columns
df_expr <- as.data.frame(t(sub_expr))
df_expr$barcode <- rownames(df_expr)

# Merge with clinical data, remove missing values and duplicates
df_clean <- df_expr %>%
  left_join(clinical_meta, by = "barcode") %>%
  # Remove duplicated samples based on patient barcode
  mutate(patient_id = substr(barcode, 1, 12)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  # Ensure all 5 genes are present (drop any rows with NAs in gene columns)
  drop_na(all_of(target_genes))

N_samples <- nrow(df_clean)
cat("Final number of analyzed tumor samples (N):", N_samples, "\n")

# ------------------------------------------------------------------------------
# STEP 4: Check Expression Distributions
# ------------------------------------------------------------------------------
df_long <- df_clean %>%
  select(all_of(target_genes)) %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression")

# Density Plot
p_density <- ggplot(df_long, aes(x = Expression, fill = Gene)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Gene Expression Distributions (Density)", x = "log2(TPM + 1)")

# Boxplot
p_boxplot <- ggplot(df_long, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression Distributions (Boxplot)", y = "log2(TPM + 1)")

ggsave("Distribution_Density.png", p_density, width = 8, height = 5, dpi = 300)
ggsave("Distribution_Boxplot.png", p_boxplot, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# STEP 5, 6, 7 & 8: Correlation Analyses & Multiple Testing Correction
# ------------------------------------------------------------------------------
candidates <- c("SOX1", "SOX2", "NEUROG2", "SIX3")
results_list <- list()

for (gene in candidates) {
  # Pearson correlation
  pearson_res <- cor.test(df_clean$PAX6, df_clean[[gene]], method = "pearson")
  
  # Spearman correlation
  spearman_res <- cor.test(df_clean$PAX6, df_clean[[gene]], method = "spearman", exact = FALSE)
  
  results_list[[gene]] <- data.frame(
    Gene_pair = paste0("PAX6–", gene),
    N = N_samples,
    Pearson_r = pearson_res$estimate,
    Pearson_p = pearson_res$p.value,
    Spearman_rho = spearman_res$estimate,
    Spearman_p = spearman_res$p.value,
    CI_lower = pearson_res$conf.int[1],
    CI_upper = pearson_res$conf.int[2]
  )
}

results_df <- bind_rows(results_list)

# Apply Benjamini-Hochberg FDR correction on Spearman p-values
results_df$BH_adjusted_Spearman_p <- p.adjust(results_df$Spearman_p, method = "BH")

# Function to categorize correlation strength based on specified framework
categorize_strength <- function(r) {
  abs_r <- abs(r)
  case_when(
    abs_r < 0.20 ~ "Very weak",
    abs_r >= 0.20 & abs_r < 0.40 ~ "Weak",
    abs_r >= 0.40 & abs_r < 0.60 ~ "Moderate",
    abs_r >= 0.60 & abs_r < 0.80 ~ "Strong",
    abs_r >= 0.80 ~ "Very strong"
  )
}

results_df <- results_df %>%
  mutate(
    Direction = ifelse(Spearman_rho > 0, "Positive", "Negative"),
    Strength = categorize_strength(Spearman_rho),
    `95%_CI` = paste0("[", round(CI_lower, 3), ", ", round(CI_upper, 3), "]")
  ) %>%
  select(Gene_pair, N, Pearson_r, Pearson_p, Spearman_rho, Spearman_p, 
         BH_adjusted_Spearman_p, `95%_CI`, Direction, Strength)

write.csv(results_df, "Correlation_Results_Table.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 9 & 10: Generate Publication-Quality Scatter Plots
# ------------------------------------------------------------------------------
create_scatter <- function(data, cand_gene, res_df) {
  
  # Extract specific stats for annotation
  stats <- res_df %>% filter(Gene_pair == paste0("PAX6–", cand_gene))
  
  label_text <- sprintf(
    "N = %d\nSpearman rho = %.3f\np = %.2e\nFDR = %.2e",
    stats$N, stats$Spearman_rho, stats$Spearman_p, stats$BH_adjusted_Spearman_p
  )
  
  ggplot(data, aes(x = PAX6, y = .data[[cand_gene]])) +
    geom_point(alpha = 0.6, color = "darkblue", size = 1.5) +
    geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.3) +
    labs(x = "PAX6 Expression log2(TPM+1)", y = paste(cand_gene, "Expression log2(TPM+1)")) +
    theme_classic(base_size = 14) +
    annotate("text", x = min(data$PAX6), y = max(data[[cand_gene]]), 
             label = label_text, hjust = 0, vjust = 1, size = 4)
}

plot_sox1 <- create_scatter(df_clean, "SOX1", results_df)
plot_sox2 <- create_scatter(df_clean, "SOX2", results_df)
plot_neurog2 <- create_scatter(df_clean, "NEUROG2", results_df)
plot_six3 <- create_scatter(df_clean, "SIX3", results_df)

# Combine using patchwork
combined_plot <- (plot_sox1 | plot_sox2) / (plot_neurog2 | plot_six3) +
  plot_annotation(tag_levels = 'A')

ggsave("Figure1_ScatterPlots.pdf", combined_plot, width = 12, height = 10, dpi = 300)
ggsave("Figure1_ScatterPlots.png", combined_plot, width = 12, height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# STEP 11: Generate Correlation Matrix/Heatmap
# ------------------------------------------------------------------------------
corr_matrix <- cor(df_clean[, target_genes], method = "spearman")
p_val_matrix <- cor_pmat(df_clean[, target_genes], method = "spearman", exact = FALSE)

p_heatmap <- ggcorrplot(corr_matrix, 
                        method = "square", 
                        type = "lower",
                        p.mat = p_val_matrix, 
                        insig = "blank",
                        lab = TRUE,
                        colors = c("#6D9EC1", "white", "#E46726"),
                        title = "Spearman Correlation Heatmap")

ggsave("Figure2_Heatmap.pdf", p_heatmap, width = 7, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# STEP 12 & 13: Sensitivity Analysis & Outlier Detection
# ------------------------------------------------------------------------------
# Sensitivity Analysis: Comparing Pearson vs Spearman is inherently done 
# in the results_df. We now check for influential outliers using Cook's distance.

outlier_results <- list()
df_clean_no_outliers <- df_clean

for (gene in candidates) {
  model <- lm(as.formula(paste(gene, "~ PAX6")), data = df_clean)
  cooksd <- cooks.distance(model)
  
  # Threshold rule of thumb: 4 / N
  influential <- as.numeric(names(cooksd)[(cooksd > (4 / N_samples))])
  
  if (length(influential) > 0) {
    df_temp <- df_clean[-influential, ]
    spearman_sens <- cor.test(df_temp$PAX6, df_temp[[gene]], method = "spearman", exact = FALSE)
    
    outlier_results[[gene]] <- data.frame(
      Gene_pair = paste0("PAX6–", gene),
      Outliers_Removed = length(influential),
      New_N = nrow(df_temp),
      Orig_Spearman_rho = results_df$Spearman_rho[results_df$Gene_pair == paste0("PAX6–", gene)],
      Sens_Spearman_rho = spearman_sens$estimate,
      Sens_Spearman_p = spearman_sens$p.value
    )
  } else {
    outlier_results[[gene]] <- data.frame(
      Gene_pair = paste0("PAX6–", gene),
      Outliers_Removed = 0, New_N = N_samples,
      Orig_Spearman_rho = NA, Sens_Spearman_rho = NA, Sens_Spearman_p = NA
    )
  }
}
outlier_df <- bind_rows(outlier_results)
write.csv(outlier_df, "Sensitivity_Analysis_Outliers.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 14: Subgroup Stratification by Histological Subtype
# ------------------------------------------------------------------------------
# Ensure groups have a sufficient sample size (e.g., N > 20)
if ("primary_diagnosis" %in% colnames(df_clean)) {
  subgroup_analysis <- df_clean %>%
    group_by(primary_diagnosis) %>%
    filter(n() > 20) %>%
    summarise(
      N = n(),
      SOX1_rho = cor(PAX6, SOX1, method = "spearman"),
      SOX2_rho = cor(PAX6, SOX2, method = "spearman"),
      NEUROG2_rho = cor(PAX6, NEUROG2, method = "spearman"),
      SIX3_rho = cor(PAX6, SIX3, method = "spearman"),
      .groups = "drop"
    )
  
  write.csv(subgroup_analysis, "Subgroup_Analysis_Histology.csv", row.names = FALSE)
} else {
  message("Clinical metadata 'primary_diagnosis' not found; skipping subgroup analysis.")
}