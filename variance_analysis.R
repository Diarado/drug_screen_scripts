# Set working directory
setwd("D:/Jiahe/columbia/mr_analysis/narnea_results/top250_1000btstrap")

# Load required libraries
library(stats)
library(qvalue)

###########################################
# Part 0: Generate Dummy PC1 Data
###########################################
# Generate dummy PC1 values for each round and plate
set.seed(42)  # For reproducibility

# Generate PC1 for each plate in each round
PC1_by_plate <- list()
for(round_num in 1:12) {
  PC1_by_plate[[round_num]] <- list(
    rep1 = rnorm(96),    # 96 wells for plate 1
    rep2 = rnorm(97),    # 97 wells for plate 2
    rep3 = rnorm(95)     # 95 wells for plate 3
  )
}

# Generate PC1 for each complete run (all plates combined)
PC1_by_run <- list()
for(round_num in 1:12) {
  PC1_by_run[[round_num]] <- rnorm(288)  # 288 total wells per round
}

###########################################
# Part 1: Data Loading and Processing
###########################################

# Function to load count data for a specific round
load_count_data <- function(round_num) {
  round_str <- sprintf("%02d", round_num)
  file_path <- file.path(paste0("round", round_str), "narnea", 
                         paste0("drs_round", round_str, "_raw_counts.tsv"))
  
  tryCatch({
    cnt_data <- read.table(file_path, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
    return(cnt_data)
  }, error = function(e) {
    cat(sprintf("Error loading data for round %s: %s\n", round_str, e$message))
    return(NULL)
  })
}

# Function to split data into replicates
split_replicates <- function(cnt_data) {
  if (is.null(cnt_data)) return(NULL)
  
  list(
    rep1 = cnt_data[, 1:96],
    rep2 = cnt_data[, 97:193],
    rep3 = cnt_data[, 194:288]
  )
}

###########################################
# Part 2: PC1 Correlation Analysis
###########################################

# Function to calculate correlations with PC1
calculate_pc1_correlations <- function(expression_data, pc1) {
  # Calculate Spearman correlation for each gene
  results <- t(apply(expression_data, 1, function(x) {
    test_result <- cor.test(x, pc1, method = "spearman")
    c(correlation = test_result$estimate,
      p_value = test_result$p.value)
  }))
  
  # Convert to data frame and add gene names
  results_df <- data.frame(
    gene = rownames(expression_data),
    correlation = results[, "correlation"],
    p_value = results[, "p_value"]
  )
  
  # Calculate FDR
  results_df$q_value <- qvalue(results_df$p_value)$qvalues
  
  # Filter for significant genes (FDR < 0.05)
  significant_genes <- results_df[results_df$q_value < 0.05, ]
  
  return(list(
    all_results = results_df,
    significant_genes = significant_genes
  ))
}

###########################################
# Part 3: Gene List Comparison
###########################################

# Function to compare gene lists
compare_gene_lists <- function(variable_genes_lists, pc1_genes) {
  lapply(variable_genes_lists, function(var_genes) {
    common <- intersect(var_genes, pc1_genes)
    list(
      common_genes = common,
      count = length(common),
      percent_overlap = length(common) / length(var_genes) * 100
    )
  })
}

###########################################
# Part 4: Mean Expression Analysis
###########################################

# Function to analyze mean expression
analyze_mean_expression <- function(expression_data, gene_list) {
  # Calculate mean expression vector
  mean_expr <- colMeans(expression_data[gene_list, , drop = FALSE])
  
  # Calculate correlations with mean expression
  results <- t(apply(expression_data[gene_list, , drop = FALSE], 1, function(x) {
    test_result <- cor.test(x, mean_expr, method = "spearman")
    c(correlation = test_result$estimate,
      p_value = test_result$p.value)
  }))
  
  # Create results data frame
  results_df <- data.frame(
    gene = gene_list,
    correlation = results[, "correlation"],
    p_value = results[, "p_value"]
  )
  
  # Calculate FDR
  results_df$q_value <- qvalue(results_df$p_value)$qvalues
  
  return(results_df)
}

###########################################
# Part 5: Main Analysis Pipeline
###########################################

# Load previously identified variable genes
variable_genes <- list(
  combined_100 = read.table("common_genes_combined_100.txt", stringsAsFactors = FALSE)$V1,
  rep1_100 = read.table("common_genes_rep1_100.txt", stringsAsFactors = FALSE)$V1,
  rep2_100 = read.table("common_genes_rep2_100.txt", stringsAsFactors = FALSE)$V1,
  rep3_100 = read.table("common_genes_rep3_100.txt", stringsAsFactors = FALSE)$V1
)

# Initialize results storage
pc1_correlations_by_plate <- list()
pc1_correlations_by_run <- list()
mean_expr_analysis <- list()

# Process each round
for(round_num in 1:12) {
  # Load data
  cnt_data <- load_count_data(round_num)
  if (is.null(cnt_data)) next
  
  # Split into replicates
  replicates <- split_replicates(cnt_data)
  
  # Calculate PC1 correlations for each plate
  for(rep_name in names(replicates)) {
    pc1_plate <- PC1_by_plate[[round_num]][[rep_name]]  # Assuming PC1_by_plate exists
    correlations <- calculate_pc1_correlations(replicates[[rep_name]], pc1_plate)
    pc1_correlations_by_plate[[paste0("round", round_num, "_", rep_name)]] <- correlations
  }
  
  # Calculate PC1 correlations for combined run
  pc1_run <- PC1_by_run[[round_num]]  # Assuming PC1_by_run exists
  correlations_run <- calculate_pc1_correlations(cnt_data, pc1_run)
  pc1_correlations_by_run[[paste0("round", round_num)]] <- correlations_run
  
  # Analyze mean expression for variable genes
  for(list_name in names(variable_genes)) {
    mean_expr_results <- analyze_mean_expression(cnt_data, variable_genes[[list_name]])
    mean_expr_analysis[[paste0("round", round_num, "_", list_name)]] <- mean_expr_results
  }
}

###########################################
# Part 6: Results Summary and Export
###########################################

# Find common significant genes across plates
significant_genes_by_plate <- lapply(pc1_correlations_by_plate, function(x) x$significant_genes$gene)
common_significant_plate <- Reduce(intersect, significant_genes_by_plate)

# Find common significant genes across runs
significant_genes_by_run <- lapply(pc1_correlations_by_run, function(x) x$significant_genes$gene)
common_significant_run <- Reduce(intersect, significant_genes_by_run)

# Compare variable genes with PC1-correlated genes
plate_comparisons <- compare_gene_lists(variable_genes, common_significant_plate)
run_comparisons <- compare_gene_lists(variable_genes, common_significant_run)

# Save results
saveRDS(list(
  pc1_correlations_plate = pc1_correlations_by_plate,
  pc1_correlations_run = pc1_correlations_by_run,
  mean_expr_analysis = mean_expr_analysis,
  common_significant_plate = common_significant_plate,
  common_significant_run = common_significant_run,
  plate_comparisons = plate_comparisons,
  run_comparisons = run_comparisons
), "gene_analysis_results.rds")

# Generate summary report
sink("analysis_summary.txt")
cat("Analysis Summary\n")
cat("================\n\n")
cat(sprintf("Common significant genes across plates: %d\n", length(common_significant_plate)))
cat(sprintf("Common significant genes across runs: %d\n", length(common_significant_run)))
cat("\nOverlap with variable genes:\n")
for(list_name in names(variable_genes)) {
  cat(sprintf("\n%s:\n", list_name))
  cat(sprintf("  Plate analysis overlap: %d genes (%.2f%%)\n", 
              plate_comparisons[[list_name]]$count,
              plate_comparisons[[list_name]]$percent_overlap))
  cat(sprintf("  Run analysis overlap: %d genes (%.2f%%)\n", 
              run_comparisons[[list_name]]$count,
              run_comparisons[[list_name]]$percent_overlap))
}
sink()