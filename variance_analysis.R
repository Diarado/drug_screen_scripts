# Set working directory
setwd("D:/Jiahe/columbia/mr_analysis/narnea_results/top250_1000btstrap")

# Initialize lists to store top genes from each round
# One list for combined data, three lists for individual replicates
top_genes_combined <- list()
top_genes_rep1 <- list()
top_genes_rep2 <- list()
top_genes_rep3 <- list()

# Loop through rounds 1 to 12
for(i in 1:12) {
  # Format round number (01, 02, ..., 12)
  round_num <- sprintf("%02d", i)
  cur_round <- paste0("round", round_num)
  
  # Construct file path
  file_path <- file.path(cur_round, "narnea", paste0("drs_", cur_round, "_raw_counts.tsv"))
  
  tryCatch({
    # Read count file
    cnt_data <- read.table(file_path, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
    
    # Split data into replicates
    cnt_rep1 <- cnt_data[, 1:96]        # Replicate 1
    cnt_rep2 <- cnt_data[, 97:193]      # Replicate 2
    cnt_rep3 <- cnt_data[, 194:288]     # Replicate 3
    
    # Calculate variances and get top 1000 genes for each replicate
    # Replicate 1
    gene_variances1 <- apply(cnt_rep1, MARGIN=1, var)
    top_1000_genes1 <- names(sort(gene_variances1, decreasing=TRUE)[1:1000])
    top_genes_rep1[[round_num]] <- top_1000_genes1
    
    # Replicate 2
    gene_variances2 <- apply(cnt_rep2, MARGIN=1, var)
    top_1000_genes2 <- names(sort(gene_variances2, decreasing=TRUE)[1:1000])
    top_genes_rep2[[round_num]] <- top_1000_genes2
    
    # Replicate 3
    gene_variances3 <- apply(cnt_rep3, MARGIN=1, var)
    top_1000_genes3 <- names(sort(gene_variances3, decreasing=TRUE)[1:1000])
    top_genes_rep3[[round_num]] <- top_1000_genes3
    
    # Calculate variance for combined data (original approach)
    gene_variances_combined <- apply(cnt_data, MARGIN=1, var)
    top_1000_genes_combined <- names(sort(gene_variances_combined, decreasing=TRUE)[1:1000])
    top_genes_combined[[round_num]] <- top_1000_genes_combined
    
    cat(sprintf("Processed round %s successfully\n", round_num))
    
  }, error = function(e) {
    cat(sprintf("Error processing round %s: %s\n", round_num, e$message))
  })
}

# Find common genes across all rounds for each replicate
common_genes_rep1 <- Reduce(intersect, top_genes_rep1)
common_genes_rep2 <- Reduce(intersect, top_genes_rep2)
common_genes_rep3 <- Reduce(intersect, top_genes_rep3)
common_genes_combined <- Reduce(intersect, top_genes_combined)

# Save results for each replicate
write.table(sort(common_genes_rep1),
            file="common_genes_rep1_100.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(sort(common_genes_rep2),
            file="common_genes_rep2_100.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(sort(common_genes_rep3),
            file="common_genes_rep3_100.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(sort(common_genes_combined),
            file="common_genes_combined_100.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Print summary statistics
cat("\nNumber of common genes in each analysis:\n")
cat(sprintf("Replicate 1: %d genes\n", length(common_genes_rep1)))
cat(sprintf("Replicate 2: %d genes\n", length(common_genes_rep2)))
cat(sprintf("Replicate 3: %d genes\n", length(common_genes_rep3)))
cat(sprintf("Combined data: %d genes\n", length(common_genes_combined)))

# Find genes common across all replicates
all_common_genes <- Reduce(intersect, list(common_genes_rep1, 
                                           common_genes_rep2, 
                                           common_genes_rep3))

write.table(sort(all_common_genes),
            file="common_genes_across_all_reps_100.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

cat(sprintf("\nGenes common across all replicates: %d\n", length(all_common_genes)))
 
# 
# This is quite interesting.  I would like to further explore the genes that characterize these lists, and compare with the prior work on 1st PC you did a while ago.
# 
# 
# 
# My memory is that you did the 1st PC two ways – per plate and per run (including all three plates together), similarly to how you generated the most variable genes in this exercise.  First, I would like you to identify genes that correlate with the 1st PC (using Spearman’s correlation) and pass FDR less than 0.05.  Do this for each plate with the 1st PC that you calculated for that plate, and then identify if there are any genes that appear on every list from every plate.  Then do the same thing per run, where you include all three plates in one matrix, identify genes that correlate with the 1st PC for that run across all three plates (again, using Spearman’s with FDR less than 0.05), and then see if there are any genes in common across the lists for every run.
# 
# 
# 
# In addition to seeing what are on these lists, I would be interested in knowing if any of the genes on the 4 lists you generated of highly variable genes (top 100 and top 1000 for each run or each plate) are on either of the 2 lists you are generating in the above paragraph.  So in other words, all 4 lists of highly variable genes will be compared to the 1st PC gene list by plate, and then again to the 1st PC gene list by run, and identify how many (and which) genes are in common in each comparison.
# 
# 
# 
# Assuming we use the highly variable genes in the drug screen, I am interested in the best way to remove the effect of these genes.  One way would be to generate a mean gene expression vector and remove the effect of this vector, but I would like to study the mean gene expression vector a bit first.  So for each of the 4 highly variable gene lists, please create a mean gene expression vector for each well of the drug screen (i.e. just linearly combine the value of each gene on the list to generate a single value for each well).  Then determine what the Spearman’s correlation is between each gene on the list and the mean gene expression vector of all genes on the list; calculate this correlation using all wells in the drug screen (including p-value, and FDR adjusted p-value).  Do this for each of the 4 highly variable lists.  What I’m trying to do here is determine how well a mean gene expression vector is accounting for the variability of each gene that makes it up.
# 

