source("utils.R")

###############################################################################
## Example
###############################################################################

reference <- scPred::pbmc_1
query <- scPred::pbmc_2

plot_value_counts(reference, "cell_type")
plot_value_counts(query, "cell_type")

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

query <- NormalizeData(query)

DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

reference <- getFeatureSpace(reference, "cell_type")

for (m in models_list) {
  cat("Training model ", m, ".\n")
  reference <- trainModel(reference, model = 'svmRadial')
  
  get_probabilities(reference) %>% head()
  get_scpred(reference)
  
  query <- scPredict(query, reference)
  DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
  
  crossTab(query, "cell_type", "scpred_prediction")
  
  save_metrics(query, m, "example_run")
}



