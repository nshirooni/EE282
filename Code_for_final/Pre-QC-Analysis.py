#Load in packages
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import anndata as ad
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython

#Read in data
adata = sc.read_10x_h5("/pub/nshiroon/scrna/sc/e11__3_4_i_b.h5")

#Calculates mitochondrial percentage as well as light filtering
def calculate_mito_percentage(adata, mito_prefixes=('^MT-', '^mt-', '^Mito-')):
    import re
    
    mito_genes = adata.var_names[adata.var_names.str.contains('|'.join(mito_prefixes), regex=True)]
    
    adata.var['mt'] = adata.var_names.isin(mito_genes)
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, log1p=True, percent_top= [20])
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    return adata

#Ensure unique genes and run defined function for mitochondrial percentage and light filtering
adata.var_names_make_unique()
calculate_mito_percentage(adata) 

#Save raw counts to a layer 
adata.layers["counts"] = adata.X.copy()

#Initial downstream analysis
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')
genes_to_remove = ['eGFP','tdTomato']
adata.var['highly_variable'] = ~adata.var.index.isin(genes_to_remove) & adata.var['highly_variable']
sc.pl.highly_variable_genes(adata, show=False)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution= .5)
sc.pl.umap(adata, color='leiden')

#Save data for R transfer
adata.write('pre_QC.h5ad')

# Convert to R for ggplot 2 analysis
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
ro.r('''
#Load packages
    library(zellkonverter)
    library(Seurat)
    library(ggplot2)
    library(viridis) 

#Load the AnnData object
  object.sce <- readH5AD('pre_QC.h5ad')
  
  # Convert to Seurat object
  object <- Seurat::as.Seurat(object.sce, counts = "counts", data = "X")
  
  # Rename assay
  object <- RenameAssays(object, originalexp = "RNA")
  # Clean up the AnnData object
  rm(object.sce)

#Ensure cells are identified by unbiased cluters
seurat_obj <- object
Idents(seurat_obj) <- seurat_obj@meta.data$leiden

    # Plot 1: Histogram of Total Counts
    p1 <- ggplot(seurat_obj@meta.data, aes(x = total_counts)) +
      geom_histogram(bins = 50, fill = "blue", color = "black") +
      labs(title = "Distribution of Total Counts per Cell",
           x = "Total UMI Counts",
           y = "Number of Cells") +
      theme_minimal()

    # Plot 2: Violin Plot for % Mitochondrial Genes
    p2 <- VlnPlot(seurat_obj, features = 'pct_counts_mt')

    # Plot 3: Scatter Plot of Total Counts vs Total Genes
    p3 <- ggplot(seurat_obj@meta.data, aes(x = total_counts, y = n_genes_by_counts, color = pct_counts_mt)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_viridis_c(option = "plasma") +
      labs(title = "Total Counts vs Total Genes",
           x = "Total UMI Counts",
           y = "Number of Genes Detected",
           color = "% Mitochondrial Genes") +
      theme_minimal() +
      geom_hline(yintercept = 200, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 500, linetype = "dashed", color = "gray50")

    #Seurat plots for visualizing metrics in low-dimensional space
    p4 <- FeaturePlot(seurat_obj, features = 'total_counts')
    p5 <- FeaturePlot(seurat_obj, features = 'n_genes_by_counts')
    p6 <- FeaturePlot(seurat_obj, features = 'pct_counts_mt')

    # save plots
ggsave("histogram.png", p1)
ggsave("violin_plot.png", p2)
ggsave("scatter_plot.png", p3)
ggsave("total_counts.png", p4)
ggsave("n_genes_by_counts.png", p5)
ggsave("pct_counts_mt.png", p6)
''')

#For interactive usage within a jupyter notebook
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#For first three plots
# Create a figure with subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Read and display images
images = ['histogram.png', 'violin_plot.png', 'scatter_plot.png']
for i, img_path in enumerate(images):
    img = mpimg.imread(img_path)
    axs[i].imshow(img)
    axs[i].axis('off')  # Turn off axis numbers

plt.tight_layout()
plt.show()

#For last three plots
# Create a figure with subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Read and display images
images = ['total_counts.png', 'n_genes_by_counts.png', 'pct_counts_mt.png']
for i, img_path in enumerate(images):
    img = mpimg.imread(img_path)
    axs[i].imshow(img)
    axs[i].axis('off')  # Turn off axis numbers

plt.tight_layout()
plt.show()

#Save unfiltered data prior to QC
adata.write('/pub/nshiroon/scrna/pre_QC_E11/e11__SC_i_a.h5ad')