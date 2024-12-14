# Directory containing the files
dir <- '/dfs6/pub/nshiroon/RSEM_results/3_v_neg'

# List all files in the directory that match the pattern '*.genes.results'
files <- list.files(path = dir, pattern = "\\.genes\\.results$", full.names = TRUE)

# Extract sample names by removing the '.genes.results' suffix from filenames
sample_names <- sub("\\.genes\\.results$", "", basename(files))

# Assign sample names to the files
names(files) <- sample_names

# Print to verify
print(files)


library(tximport)
#load in rsem output 
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

txi.rsem$length[txi.rsem$length == 0] <- 1


# Define conditions manually
condition_map <- c(
  "cn3_A" = "CN3",
  "cn3_B" = "CN3",
  "cn3_neg_A" = "neg",
  "cn3_neg_B" = "neg"
)

# Assign conditions based on the mapping
sample_metadata <- data.frame(
  sample_id = sample_names,
  condition = condition_map[sample_names]
)

# View the metadata
print(sample_metadata)

#Load in data in DESeq2 format
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi.rsem, sample_metadata, ~condition)

#Convert ensmbl IDs to gene names
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=rownames(dds),mart= mart)
G_list = G_list[!G_list$external_gene_name=="",]
new_ra = rownames(dds)
new_ra[match(G_list$ensembl_gene_id,new_ra)] = G_list$external_gene_name
new_ra[duplicated(new_ra)] = rownames(dds)[duplicated(new_ra)]
rownames(dds) = new_ra

#Run DGE and pull results
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "CN3", "neg"))
res

# Select top 50 significant genes
top_genes <- rownames(res[order(res$padj),])[1:50]

# Extract normalized counts
library(pheatmap)
normalized_counts <- counts(dds, normalized = TRUE)
top_counts <- normalized_counts[top_genes, ]

# Scale and plot heatmap
pheatmap(
  scale(t(top_counts)),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Top 50 Differentially Expressed Genes"
)