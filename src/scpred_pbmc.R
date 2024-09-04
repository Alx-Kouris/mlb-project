###############################################################################
## Running the below code will load the 10X PBMC dataset as reference,
## split into 3 parts (validation, train, test) and then use the train part
## for classifying the PBMC test set, the PBMC validation set and finally
## the FANTOM5 original dataset
###############################################################################

source("utils.R")
library(caret)

###############################################################################
## Train-test split of query datatest
###############################################################################

pbmc <- load_query_dataset()

# Step 1: Extract the metadata
metadata <- pbmc@meta.data

# Step 2: Create a stratified split
set.seed(42)  # For reproducibility
query_index <- createDataPartition(metadata$cell_type, p = 0.5, list = FALSE)

# Step 3: Split the Seurat object
pbmc_query <- subset(pbmc, cells = colnames(pbmc)[query_index])
pbmc_model <- subset(pbmc, cells = colnames(pbmc)[-query_index])

# Step 4: Split again
metadata <- pbmc_model@meta.data
train_index <- createDataPartition(metadata$cell_type, p = 0.5, list = FALSE)

pbmc_train <- subset(pbmc_model, cells = colnames(pbmc_model)[train_index])
pbmc_test <- subset(pbmc_model, cells = colnames(pbmc_model)[-train_index])

plot_value_counts(pbmc_train, "cell_type")
plot_value_counts(pbmc_test, "cell_type")

###############################################################################
## Process reference dataset
###############################################################################

pbmc_train <- pbmc_train %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(pbmc_train, group.by = "cell_type", label = TRUE, repel = TRUE)

pbmc_train <- getFeatureSpace(pbmc_train, "cell_type")

###############################################################################
## Classify test dataset
###############################################################################

pbmc_test <- NormalizeData(pbmc_test)

for (m in models_list) {
  cat("Training model ", m, ".\n")
  pbmc_train <- trainModel(pbmc_train, model = m)
  get_classifiers(pbmc_train)
  get_scpred(pbmc_train)
  
  pbmc_test <- scPredict(pbmc_test, pbmc_train)
  
  save_metrics(pbmc_test, m, "pbmc_test")
}

###############################################################################
## Classify validation dataset
###############################################################################

pbmc_query <- NormalizeData(pbmc_query)

for (m in models_list) {
  cat("Training model ", m, ".\n")
  pbmc_train <- trainModel(pbmc_train, model = m)
  get_classifiers(pbmc_train)
  get_scpred(pbmc_train)
  
  pbmc_query <- scPredict(pbmc_query, pbmc_train)
  
  save_metrics(pbmc_query, m, "pbmc_validation")
}

###############################################################################
## Classify FANTOM5 dataset
###############################################################################

data <- read.csv("gene_expressions.csv", row.names=1)
reference <- CreateSeuratObject(counts = data)

cell_type <- gsub("\\.\\d+$", "", rownames(reference[[]]))
reference$cell_type <- cell_type

reference@meta.data$cell_type <- cell_type_map[reference@meta.data$cell_type]

reference <- NormalizeData(reference)
for (m in models_list) {
  cat("Training model ", m, ".\n")
  pbmc_train <- trainModel(pbmc_train, model = m)
  get_classifiers(pbmc_train)
  get_scpred(pbmc_train)
  
  reference <- scPredict(reference, pbmc_train)
  
  save_metrics(reference, m, "fantom5_over_pbmc")
}