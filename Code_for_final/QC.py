#Load packages
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import median_abs_deviation as mad
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
import celltypist
from celltypist import models
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=50,
    facecolor="white",
    frameon=False,
)

#Read in data
adata = sc.read_h5ad('/pub/nshiroon/scrna/pre_QC_E11/e11__SC_i_a.h5ad')

#Define functions for outlier detection QC based on Median Absolute Deviation 
def mad_outlier(adata, metric, nmads, upper_only = False):
    M = adata.obs[metric]
    
    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
    
    return (M > np.median(M) + nmads * mad(M))
def auto_QC(adata):
            outliers = mad_outlier(adata, 'log1p_total_counts', 5) +\
            mad_outlier(adata, 'log1p_n_genes_by_counts', 5) +\
            mad_outlier(adata, 'pct_counts_in_top_20_genes', 5) +\
            mad_outlier(adata, 'pct_counts_mt', 3, upper_only= True)
            adata = adata[~outliers]               
            return adata

#Load Packages
ro.r('''
library(Seurat)
library(scater)
library(scDblFinder)
library(BiocParallel)
''')

#Function for scDblFinder
ro.r('''
process_anndata <- function(data_mat) {
    set.seed(123)
    sce <- scDblFinder(SingleCellExperiment(list(counts=data_mat)))
    doublet_score <- sce$scDblFinder.score
    doublet_class <- sce$scDblFinder.class
    return(list(doublet_score=doublet_score, doublet_class=doublet_class))
}
''')

#Use scanpy matrix as input for R function and save doublets to scanpy metadata for analysis
adata.X = adata.layers["counts"].copy()
data_mat = adata.X.T  

#Call the R function
result = ro.r['process_anndata'](data_mat)

doublets = result['doublet_class']

adata.obs['doublets'] = doublets

#Filter doublets and call QC function 
adata = adata[adata.obs['doublets'] == 'singlet']
adata = auto_QC(adata) 

#Create function for annotation of cells
def annotate_cells(adata):
    model = models.Model.load(model = 'Developing_Mouse_Brain.pkl')
    predictions = celltypist.annotate(adata, model = 'Developing_Mouse_Brain.pkl', majority_voting = True)
    adata = predictions.to_adata()
    adata.obs['celltypes'] = adata.obs['majority_voting']
    sc.pl.umap(adata, color = 'celltypes')
    return adata

#Downstream analysis + Annotations
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor='seurat_v3', layer='counts')
genes_to_remove = ['eGFP','tdTomato']
adata.var['highly_variable'] = ~adata.var.index.isin(genes_to_remove) & adata.var['highly_variable']
sc.pl.highly_variable_genes(adata, show=False)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution= .5)
sc.pl.umap(adata, color='leiden')
adata = annotate_cells(adata) 

#Save QC'd data
adata.write('/pub/nshiroon/scrna/clean_E11/e11__3_4_i_b.h5ad')

# Convert AnnData to R (R DataFrame is required for ggplot2)
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
ro.r('''
    library(zellkonverter)
    library(Seurat)
    library(ggplot2)
    library(viridis)  # For color palettes

# Load the AnnData object
    object.sce <- readH5AD('/pub/nshiroon/scrna/clean_E11/e11__3_4_i_b.h5ad')
  
  # Convert to Seurat object
  object <- Seurat::as.Seurat(object.sce, counts = "counts", data = "X")
  
  # Rename assay if needed
  object <- RenameAssays(object, originalexp = "RNA")
  # Clean up the AnnData object
  rm(object.sce)

seurat_obj <- object
#Ensure identities for cells are based on unbiased leiden clusters
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
    # Plot 4: Seurat Plot for celltypes
    p4 <- DimPlot(seurat_obj, group.by = 'celltypes')

    #Save plots
ggsave("histogram.png", p1)
ggsave("violin_plot.png", p2)
ggsave("scatter_plot.png", p3)
ggsave("celltypes.png", p4)
''')

#Plots for visualization in jupyter notebook
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Create a figure with subplots
fig, axs = plt.subplots(1, 3, figsize=(15,5))

# Read and display images
images = ['histogram.png', 'violin_plot.png', 'scatter_plot.png']
for i, img_path in enumerate(images):
    img = mpimg.imread(img_path)
    axs[i].imshow(img)
    axs[i].axis('off')  # Turn off axis numbers

plt.tight_layout()
plt.show()