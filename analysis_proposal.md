## Analysis Proposal for EE282 Final Project
Author: Nima Shirooni

# Introduction
Motor neurons are key celltypes in the transmission of signals from the brain and spinal cord to target muscles/glands.
The developmental landscape of motor neurons are yet to be fully discovered, particularly more work has been done on development 
of motor neurons in the spinal cord compared to motor neurons in the developing brainstem. While the biological aspect 
of development is a key motivation there is also a medical component where rare congenital neurologic diseases are caused by 
disrupted or abberant inneravation of motor neurons. This is often due to differentially affected motor neuron subpopulations
in which these subpopulations and potential candidate gene markers for them have not yet been fully defined. This motivated
the usage of single-cell RNA-sequencing technology to explore the hetereogeneity of motor neurons. Matthew Rose, PI of the Rose Lab,
combined multiple mouse genetic reporter lines with intersectional temporal (embryonic days E9.5 to E18.5) and spatial transcriptomics (single cell/nuclei RNA-seq, Bulk RNA-seq, and Slide-Seq) 
to isolate and compare eight distinct mouse MN populations: the three oculomotor nuclei (CN3, CN4, CN6) and the other primary MN types (CN5, CN7, CN10, CN12 and spinal MNs). 
We propose to use the following data to create a transcriptomic atlas of the developing mouse brainstem with a focus on motor neurons. While the goal is ambitious and is proposed as a long-term goal, 
the short-term goal and the scope of the proposal is to study the motor neurons that mediate eye movement in the brainstem, the oculomotor neurons. These datasets are yet to be 
published but this work would be towards a publisable analysis and release of these datasets to any interested researcher.

# Goals
The key datasets for analysis will be all single cell/nuclei RNA-sequencing datasets from E9.5-E18.5, and 
bulk RNA-seq datasets from E11.5. The first goal is to integrate all sc/snRNA-sequencing datasets
and extract clusters related to the ocular motor neurons. The second goal is to analyze each ocular motor neuron in isolation and detect 
the various subpopulations per ocular motor neurons and there resulting differentially expressed markers (genes) and verify
expression of the in-situ expression using the Allen Brain Atlas Developmental (Allen Institute, 2004) Mouse in-situ resource. The third goal is to find 
"hidden markers" in the Bulk RNA-seq as well as full-length transcriptome RNA-seq data in which due to sequencing depth of
single-cell data, some markers will not be visible or very lowly expressed, in this case we will lose heterogeneity but pick up 
more genetic markers. 

# Methods
The first goal will involve dealing with the preprocessing (Quality Control), cell type identification, and batch integration 
of the sc/snRNA-sequencing data in which this data is formatted as a gene by cell counts matrix. 
For this goal both Python and R will be utilized through their 
respective bioinformatics analysis platforms, Scanpy (Wolf et al, 2018) and Seurat (Hao et al, 2023).
The datasets were generated utilizing 10X Droplet Sequencing in which a microfluidics chamber with oil releases cells from suspension
as well as mRNA capture beads to capture cells in a droplet and sequence their mRNA's. Three key factors are indicated in quality control,
Ambient RNA, a phenomenon in which cells do not survive dissociation and prematurely release mRNA prior to droplet capture. Doublets in which
multiple cells may be captured in a droplet. Lastly low-quality or non-viable cells which should be excluded from analysis. Ambient RNA decontamination 
will be done using CellBender (Fleming et al, 2023) which is a command-line tool in estimating ambient RNA contamination. Seurat will then input the matrix generated from
CellBender for further quality control. Threee key summary statistics will be plotted for visualization of cell quality, these three are, the amount of 
genes expressed per cell, the total counts per cell, and the mitochondrial expression percentage per cell (cell viability). These will be plotted individually
via the ViolinPlot function in Seurat to show distributions. To visualize them jointly ggplot2 will be used to create a scatter plot. Cells that 
show low genes expressed per cell, low total counts that are not biologically relevant will be removed. Further cells with high mitochondrial percentages
will be removed. Next a doublet finder algorithm, scDblFinder (Germain et al, 2022),  will be used to account for doublets, while one could assume doublets would be double the 
distribution of typical counts per cell, this is not always the case and where a doubletfinder shines. The data will then be normalized to account for 
sequencing depth, then the top 100-2000 variable genes (depending on dataset) will then be selected for downstream analysis. Dimensionality reduction will then
be applied, in which, high dimensional data is then simplified down into a low-dimensional 2D representation of clusters in which cells with similar transcriptional profiles will
cluster together. The low-dimensional space is then visualized via an algorithm in Seurat named UMAP which can be visualized using a dimensionality reduction plot, and a gene 
of interest can be viewed via feature plot in Seurat, a gene search with Isl1 will show clusters related to motor neurons, as Isl1 is differentially and highly expressed in 
motor neuron populations. To label cells CellTypist (Dominguez et al, 2022) a package in Python will be applied. This package has various reference datasets, such as developing mouse brain,
this can then be used to use the lab data as a query to the CellTypist developing mouse brain reference dataset to then label the query dataset based on the reference dataset.
To create a transcriptomic atlas all samples must be merged together in one low-dimensional space. This can be challenging due to what is known
as the "batch effect" this occurs when two similar samples are processed differently. An example of this is single cell and single nuclei, the same celltypes between them 
may cluster separately due to the sourcing of mRNA causing this batch effect. To combat this we propose using a deep learning algorithm, scDREAMER (Shree et al, 2023), in brevity 
this algorithm takes datasets of different conditions such as single cell and single nuclei and tries to create a low-dimensional space invariant of that condition.
After atlas creation the cells relating to the ocular motor neurons will then be extracted and analyzed separately. scDREAMER will then be used on the ocular motor neurons
to find unique subclusters and those subclusters will then be used for differential gene expression testing in R using MAST (Finak et al, 2015) to find unique markers that can then be 
verified for anatomical expression in the Allen Brain Atlas. Differentially expressed genes will also be visualized using pheatmap (Kolde, 2018) in R. 
Lastly Bulk and full-length transcriptome sequencing datasets have yet to be aligned to the lab reference genome. Using an alignment tool in the command line, 
CellRanger, these datasets will be aligned and then used for analysis. Full-length transcriptome data
will be analyzed in similar fashion to the single cell data. Bulk-Seq data will be compared against each ocular motor neuron through differential gene expression testing 
using DeSEQ2 (Love et al, 2014) in R and visualized with pheatmap. 

# Conclusion
The proposed analysis is integral to the understanding of the vast diversity and heterogeneity of motor neurons in the brainstem, in this case the ocular motor neurons. 
This work will also help in the future analyses of the other motor neuron populations as well. Initial work has already been set to establish this pipeline and the 
direction of it is clear, making this a feasible task. Once a pipeline is established the majority of the rest of the work is in it's biological discovery through
differential gene expression testing and use of the Allen Brain Atlas. Utilizing the HPC3 to streamline this pipeline will prove invaluable and usage of the 
free-gpu cluster will be especially helpful for deep-learning applications. This project will help establish an effective pipeline for any transcriptomic work
done beyond this proposal and eventually be used to publish a full thorough analysis of brainstem motor neurons in development.

# References 
Allen Institute for Brain Science (2004). Allen Mouse Brain Atlas [Developing Mouse Brain]. Available from mouse.brain-map.org.Allen Institute for Brain Science (2011).
Allen Reference Atlas – Mouse Brain [brain atlas]. Available from atlas.brain-map.org.

Ajita Shree*, Musale Krushna Pavan*, and Hamim Zafar. "scDREAMER for atlas-level integration of single-cell datasets using deep generative model paired with adversarial classifier." Nature Communications 14.1 (2023): 7781. doi: https://doi.org/10.1038/s41467-023-43590-8

Dominguez Conde et al., Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022).

Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, Srivastava A, Molla G, Madad S, Fernandez-Granda C, Satija R (2023). “Dictionary learning for integrative, multimodal and scalable single-cell analysis.” Nature Biotechnology. doi:10.1038/s41587-023-01767-y, https://doi.org/10.1038/s41587-023-01767-y.

Kolde R (2018). pheatmap: Pretty Heatmaps. R package version 1.0.

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.

MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data G Finak, A McDavid, M Yajima, J Deng, V Gersuk, AK Shalek, CK Slichter et al Genome biology 16 (1), 278

Stephen J Fleming, Mark D Chaffin, Alessandro Arduini, Amer-Denis Akkad, Eric Banks, John C Marioni, Anthony A Phillipakis, Patrick T Ellinor, and Mehrtash Babadi. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. Nature Methods, 2023. https://doi.org/10.1038/s41592-023-01943-7

Wolf, F., Angerer, P. & Theis, F. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0

