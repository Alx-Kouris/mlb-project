###############################################################################
## Running the below code will load the FANTOM5 dataset as reference,
## train the models and use them to classify the 10x PBMC dataset.
###############################################################################

source("utils.R")

###############################################################################
## Read reference dataset
###############################################################################

data <- read.csv("../data/preprocessed/gene_expressions.csv", row.names=1)
reference <- CreateSeuratObject(counts = data)

cell_type <- gsub("\\.\\d+$", "", rownames(reference[[]]))
reference$cell_type <- cell_type

reference@meta.data$cell_type <- cell_type_map[reference@meta.data$cell_type]

head(reference[[]])
plot_value_counts(reference, "cell_type")

###############################################################################
## Process reference dataset
###############################################################################

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

reference <- getFeatureSpace(reference, "cell_type")

DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

###############################################################################
## Read query dataset
###############################################################################

pbmc <- load_query_dataset()

plot_value_counts(pbmc, "cell_type")

# Remove cell types that don't exist in the reference dataset
new_cell_types <- unique(cell_type_map)
pbmc <- subset(pbmc, subset = cell_type %in% new_cell_types)

pbmc <- NormalizeData(pbmc)

###############################################################################
## Train models and classify
###############################################################################

for (m in models_list) {
  cat("Training model ", m, ".\n")
  reference <- trainModel(reference, model = m)
  get_classifiers(reference)
  get_scpred(reference)
  
  pbmc <- scPredict(pbmc, reference)
  crossTab(pbmc, "cell_type", "scpred_prediction")
  DimPlot(pbmc, group.by = "scpred_prediction", reduction = "scpred")
  
  save_metrics(pbmc, m, "pbmc_over_fantom5")
}
