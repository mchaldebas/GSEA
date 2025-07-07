# -----------------------------------------------------------------------------
# Required Libraries
# -----------------------------------------------------------------------------
# Ensure you have these installed: install.packages(c("readr", "readxl", "ggplot2", "dplyr", "ggrepel", "GSEABase", "tools"))
library(readr)       # For reading CSV files efficiently
library(readxl)      # To read the old excel file if needed, or new ones
library(ggplot2)
library(dplyr)
library(ggrepel)
library(GSEABase)
library(tools)       # For toTitleCase

# -----------------------------------------------------------------------------
# Main Analysis Function
# -----------------------------------------------------------------------------

#' Run Gene Set Enrichment Analysis for Continuous Data
#'
#' This function takes a list of continuous values per gene and performs an
#' enrichment analysis against a provided gene set collection (GMT file). It uses
#' the Wilcoxon Rank-Sum test to determine if the values for genes within a set
#' are significantly different from genes outside the set.
#'
#' @param gene_data_path Path to the input data file (CSV or TSV). The file must
#'   contain at least a gene column and a value column.
#' @param gmt_file_path Path to the GMT file containing gene sets.
#' @param gene_col A string with the name of the column containing gene symbols
#'   in your input data file. Default is "GENE".
#' @param value_col A string with the name of the column containing the continuous
#'   numeric values. Default is "VALUE".
#' @param set_min_size Minimum number of genes in a gene set to be included. Default is 30.
#' @param set_max_size Maximum number of genes in a gene set to be included. Default is 500.
#' @param p_value_threshold The adjusted p-value cutoff for significance. Default is 0.05.
#' @param top_n_labels The number of top pathways to label (for both up- and down-regulated). Default is 10.
#'
#' @return A list containing two elements:
#'   1. `results`: A data frame with the full analysis results for each gene set.
#'   2. `plot`: A ggplot object (volcano plot) visualizing the results.
#'

run_continuous_gsea <- function(gene_data_path,
                                gmt_file_path,
                                gene_col = "GENE",
                                value_col = "VALUE",
                                set_min_size = 30,
                                set_max_size = 500,
                                p_value_threshold = 0.05,
                                top_n_labels = 10) {
  
  # --- 1. Load and Prepare Data ---
  
  # Load gene sets from GMT file
  cat("Loading gene sets from:", gmt_file_path, "\n")
  gene_sets <- getGmt(gmt_file_path)
  pathWithGenes <- geneIds(gene_sets)
  
  # Filter gene sets by size
  original_set_count <- length(pathWithGenes)
  pathWithGenes <- pathWithGenes[
    sapply(pathWithGenes, length) >= set_min_size & sapply(pathWithGenes, length) <= set_max_size
  ]
  cat("Filtered gene sets by size (", set_min_size, "-", set_max_size, "): ",
      length(pathWithGenes), " of ", original_set_count, " sets remaining.\n", sep = "")
  
  # Load continuous gene data
  cat("Loading gene values from:", gene_data_path, "\n")
  # Using read_csv for flexibility, assuming CSV. Add logic for other types if needed.
  gene_values <- read_csv(gene_data_path, show_col_types = FALSE) %>%
    select(all_of(c(gene_col, value_col))) %>%
    distinct(.data[[gene_col]], .keep_all = TRUE) %>%
    rename(GENE = all_of(gene_col), VALUE = all_of(value_col)) # Standardize column names
  
  # Get the universe of genes from the filtered pathways for background comparison
  universe_genes <- unique(unname(unlist(pathWithGenes)))
  gene_values <- gene_values %>% filter(GENE %in% universe_genes)
  cat("Found", nrow(gene_values), "genes from your data that are in the provided gene sets.\n")
  
  # --- 2. Perform Statistical Test for each Gene Set ---
  
  results_list <- lapply(names(pathWithGenes), function(gene_set_name) {
    gene_set <- pathWithGenes[[gene_set_name]]
    
    # Define genes in the set and out of the set (the background)
    in_set_mask <- gene_values$GENE %in% gene_set
    
    # Get the continuous values for both groups
    values_in_set <- gene_values$VALUE[in_set_mask]
    values_out_of_set <- gene_values$VALUE[!in_set_mask]
    
    # Only proceed if both groups have enough data points (at least 2)
    if (length(values_in_set) < 2 || length(values_out_of_set) < 2) {
      return(NULL)
    }
    
    # Perform Wilcoxon Rank-Sum test (non-parametric alternative to t-test)
    test_result <- wilcox.test(values_in_set, values_out_of_set)
    
    # Calculate the effect size (difference in means)
    mean_in_set <- mean(values_in_set, na.rm = TRUE)
    mean_out_of_set <- mean(values_out_of_set, na.rm = TRUE)
    mean_difference <- mean_in_set - mean_out_of_set
    
    # Return a one-row data frame
    data.frame(
      GeneSet = gene_set_name,
      Mean_in_Set = mean_in_set,
      Mean_Difference = mean_difference,
      p_value = test_result$p.value,
      N_in_Set = length(values_in_set)
    )
  })
  
  # Combine list of data frames into a single data frame
  results <- bind_rows(results_list)
  cat("Analysis complete for", nrow(results), "gene sets.\n")
  
  # --- 3. Adjust P-values and Prepare for Plotting ---
  
  results <- results %>%
    filter(!is.na(p_value)) %>%
    mutate(
      adj_p_value = p.adjust(p_value, method = "BH"),
      log_p_value = -log10(adj_p_value)
    ) %>%
    # Clean up pathway names for readability
    mutate(
      Clean_GeneSet = gsub("HALLMARK_|GOBP_|KEGG_|REACTOME_", "", GeneSet),
      Clean_GeneSet = gsub("_", " ", Clean_GeneSet),
      Clean_GeneSet = tolower(Clean_GeneSet),
      Clean_GeneSet = toTitleCase(Clean_GeneSet)
    ) %>%
    arrange(adj_p_value)
  
  # --- 4. Generate Volcano Plot ---
  
  # Identify top pathways to label
  top_enriched <- results %>%
    filter(Mean_Difference > 0) %>%
    slice_head(n = top_n_labels)
  
  top_depleted <- results %>%
    filter(Mean_Difference < 0) %>%
    slice_head(n = top_n_labels)
  
  top_pathways_to_label <- bind_rows(top_enriched, top_depleted)
  
  # Create the plot
  volcano_plot <- ggplot(results, aes(x = Mean_Difference, y = log_p_value)) +
    geom_point(aes(color = ifelse(adj_p_value < p_value_threshold, Mean_Difference, NA)),
               size = 2, alpha = 0.7) +
    scale_color_gradient2(
      low = "blue", mid = "grey80", high = "red",
      midpoint = 0, name = "Mean Difference\n(Significant Sets)",
      na.value = "grey80"
    ) +
    labs(
      x = "Mean Difference (In-Set vs. Out-of-Set)",
      y = "-log10(Adjusted p-value)",
      title = "Gene Set Enrichment Volcano Plot",
      subtitle = paste("Based on", basename(gmt_file_path))
    ) +
    theme_minimal(base_size = 14) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_value_threshold), linetype = "dotted", color = "black") +
    geom_text_repel(
      data = top_pathways_to_label,
      aes(label = Clean_GeneSet),
      size = 3.8,
      max.overlaps = 15,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50"
    )
  
  # --- 5. Return Results ---
  return(list(results = results, plot = volcano_plot))
}

# --- Example Usage ---

# Define file paths
# NOTE: Make sure your GMT file path is correct!
gmt_file <- "data/h.all.v2024.1.Hs.symbols.gmt.txt"
gene_data_file <- "my_gene_values.csv" # The file we just created

# Run the analysis by calling the function
# We specify the column names from our sample file
analysis_output <- run_continuous_gsea(
  gene_data_path = gene_data_file,
  gmt_file_path = gmt_file,
  gene_col = "GeneSymbol",         # Column with gene names in our file
  value_col = "ExpressionLevel"    # Column with continuous values in our file
)

# The function returns a list. You can access the results and the plot.

# 1. View the top results table
print(head(analysis_output$results))

# 2. Display the plot
print(analysis_output$plot)

# 3. Save the plot to a file
ggsave(
  "continuous_gsea_volcano.pdf",
  plot = analysis_output$plot,
  width = 10,
  height = 8
)
cat("Plot saved to continuous_gsea_volcano.pdf\n")