library(Seurat)
library(dplyr)
library(scPred)
library(magrittr)
library(smotefamily)
library(ggplot2)
library(MLmetrics)
library(pROC)

###############################################################################
## Helper functions
###############################################################################

cell_type_map <- c(
  "CD14..Monocytes" = "CD14+ Monocyte",
  "CD19..B.cells" = "CD19+ B",
  "CD34..cells" = "CD34+",
  "CD4..CD25..Regulatory.T.cells" = "CD4+/CD25 T Reg",
  "CD4..CD45RA..CD25..Naïve.T.cells" = "CD4+/CD45RA+/CD25- Naive T",
  "CD4..CD45RO..Memory.T.cells" = "CD4+/CD45RO+ Memory",
  "CD4..Helper.T.cells" = "CD4+ T Helper2",
  "CD56..NK.cells" = "CD56+ NK",
  "CD8..Cytotoxic.T.cells" = "CD8+ Cytotoxic T"
)

models_list <- c("AdaBag", "glmboost","knn", "lda", "mlpML", "pcaNNet", "rf", "svmRadial")

plot_value_counts <- function(seurat_obj, column_name) {
  # Get the data from the specified column
  data_column <- seurat_obj@meta.data[[column_name]]
  
  # Create a data frame with the counts of each unique value
  value_counts <- as.data.frame(table(data_column))
  
  # Rename the columns for clarity
  colnames(value_counts) <- c("Value", "Count")
  
  # Plot the data
  ggplot(value_counts, aes(x = Value, y = Count)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5) + 
    theme_minimal() +
    xlab("Cell Types") +
    ylab("Genes") +
    ggtitle(paste("Number of genes with expression in ", column_name)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

display_cell_types <- function(seurat_obj) {
  # View the column names in the metadata
  colnames(seurat_obj@meta.data)
  
  # Replace 'column_name' with the name of your column of interest
  unique_values <- unique(seurat_obj@meta.data$cell_type)
  
  # Display the unique values
  print(unique_values)
}

save_metrics <- function(obj, model_name, filename) {
  # Assuming 'actual' and 'predicted' are your vectors of true labels and predictions
  actual <- obj$cell_type  # Replace with actual labels vector
  predicted <- obj$scpred_prediction  # Replace with predicted labels vector
  
  # Add overall accuracy (this is the same for all classes, so we can do it once)
  accuracy <- Accuracy(y_true = actual, y_pred = predicted)
  
  actual <- as.factor(actual)
  predicted <- as.factor(predicted)
  
  # Get all unique classes
  classes <- levels(actual)
  
  # Initialize an empty data frame to store metrics
  metrics_df <- data.frame()
  
  # Loop through each class to calculate recall
  for (cl in classes) {
    recall <- Sensitivity(y_true = actual, y_pred = predicted, positive = cl)
    precision <- Precision(y_true = actual, y_pred = predicted, positive = cl)
    f1_score <- F1_Score(y_true = actual, y_pred = predicted, positive = cl)
    
    # Append the metrics to the data frame
    metrics_df <- rbind(metrics_df, data.frame(
      Class = cl,
      Model = m,
      Recall = recall,
      Precision = precision,
      F1_Score = f1_score
    ))
  }
  
  # Add accuracy to all rows (since it is global, not class-specific)
  metrics_df$Accuracy <- accuracy
  
  # Define the file path
  file_path <- paste("../results/",filename,".csv", sep="")
  
  results_dir <- "../results"
  # Create the folder if it doesn't exist
  if (!file.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
    cat("Folder created at:", results_dir, "\n")
  }
  
  # Write the data to CSV (append mode)
  if (!file.exists(file_path)) {
    write.csv(metrics_df, file = file_path, row.names = FALSE)
  } else {
    write.table(metrics_df, file = file_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

create_ggplot_heatmap <- function(heatmap_data) {
  
  # Ensure the row names are set correctly
  rownames(heatmap_data) <- heatmap_data[,1]
  heatmap_data <- heatmap_data[,-1] # Remove the first column (assuming it contains row names)
  
  # Convert to long format for ggplot2
  heatmap_long <- melt(heatmap_data)
  
  # Create heatmap using ggplot2
  ggplot(heatmap_long, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    xlab("Cell Types") +
    ylab("Markers") +
    ggtitle("Heatmap of Cell Types and Markers") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = value), color = "black", size = 3)
}

load_query_dataset <- function() {
  pbmc.data <- Read10X(data.dir = "../data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc68k", min.cells = 3, min.features = 200)
  
  cell_types_data <- read.delim("../data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices/68k_pbmc_barcodes_annotation.tsv")
  pbmc@meta.data$cell_type <- cell_types_data$celltype[match(rownames(pbmc@meta.data), cell_types_data$barcodes)]
  
  return(pbmc)
}
