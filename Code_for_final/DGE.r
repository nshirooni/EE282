library(Seurat)
library(MAST)
library(zellkonverter)
library(dplyr)
library(ggplot2)

  # Load the AnnData object
object.sce <- readH5AD('/pub/nshiroon/scrna/merged_cleaned.h5ad')
  
  # Convert to Seurat object
object <- Seurat::as.Seurat(object.sce, counts = "counts", data = "X")
  
  # Rename assay if needed
object <- RenameAssays(object, originalexp = "RNA")
  
  # Clean up the AnnData object
rm(object.sce)

data <- object
saveRDS(data,'/pub/nshiroon/scrna/merged_cleaned.rds')

#Set identities to unbiased leiden clusters for DGE 
Idents(data) <- data@meta.data$leiden

#Visualize low-dimensional space
DimPlot(data, label = TRUE)

#Plot Important features denoting to Motor Neurons, CN3/4/6 respectively
FeaturePlot(data, features = c("Isl1","Tmem255a","Vip","Mical2"))

#Run DGE using MAST and save output
markers <- FindAllMarkers(data, test.use = 'MAST')
saveRDS(markers, 'markers.rds')

#Pull out top 3 DEGs per cluster
top3_genes <- markers%>% 
  group_by(cluster) %>% 
  top_n(3, avg_log2FC) %>% 
  pull(gene)
#Remove potential duplicates
top3_genes <- top3_genes[!duplicated(top3_genes)]
#Plot results
dot_plot = DotPlot(
  data,
  assay = 'RNA',
  top3_genes,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale = TRUE,
  scale.by = "radius",
)