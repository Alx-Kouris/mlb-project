###############################################################################
## Running the below code will load the augmented FANTOM5 dataset as reference,
## train the models and use them to classify the 10x PBMC dataset.
###############################################################################

source("utils.R")

###############################################################################
## Read augmented reference dataset
###############################################################################

# read augmented dataset
data <- read.csv("../data/preprocessed/gene_expressions_augmented.csv", row.names=1)
reference_aug <- CreateSeuratObject(counts = data)

cell_type <- gsub("\\.\\d+$", "", rownames(reference_aug[[]]))
reference_aug$cell_type <- cell_type

reference_aug@meta.data$cell_type <- cell_type_map[reference_aug@meta.data$cell_type]

head(reference_aug[[]])
plot_value_counts(reference_aug, "cell_type")

###############################################################################
## Process reference dataset
###############################################################################

reference_aug <- reference_aug %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(reference_aug, group.by = "cell_type", label = TRUE, repel = TRUE)

reference_aug <- getFeatureSpace(reference_aug, "cell_type")

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
  reference_aug <- trainModel(reference_aug, model = m)
  get_classifiers(reference_aug)
  get_scpred(reference_aug)
  
  pbmc <- scPredict(pbmc, reference_aug)
  crossTab(pbmc, "cell_type", "scpred_prediction")
  DimPlot(pbmc, group.by = "scpred_prediction", reduction = "scpred")
  
  save_metrics(pbmc, m, "pbmc_over_aug_fantom5")
}
