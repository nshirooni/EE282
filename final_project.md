# Exploring the transcriptomic profile of Ocular Motor Neurons during development

## Introduction

Motor neurons are key cell types in the transmission of signals from the brain and spinal cord to target muscles/glands. The developmental landscape of motor neurons are yet to be discovered, particularly more work has been done on development of motor neurons in the spinal cord compared to motor neurons in the brainstem. While the biological aspect of development is a key motivation there is also a medical component where rare congenital neurologic diseases are caused by disrupted or aberrant innervation of motor neurons. This is often due to differentially affected motor neuron subpopulations in which these subpopulations and potential candidate gene markers for them have yet to be fully defined. This motivated both the use of bulk RNA sequencing and single-cell RNA-sequencing technology to explore the genetic profiles of motor neurons and subsequently the heterogeneity of motor neurons. Matthew Rose, PI of the Rose Lab, combined multiple mouse genetic reporter lines with intersectional temporal (embryonic days E9.5 to E18.5) and spatial transcriptomics (single cell/nuclei RNA-seq, Bulk RNA-seq, and Slide-Seq) to isolate and compare eight distinct mouse MN populations: the three oculomotor nuclei (CN3, CN4, CN6) and the other primary MN types (CN5, CN7, CN10, CN12 and spinal MNs). We propose to use the following data to create a transcriptomic atlas of the developing mouse brainstem with a focus on motor neurons. While the goal is ambitious and is proposed as a long-term goal, the short-term goal and the scope of the proposal is to study the motor neurons that mediate eye movement in the brainstem, the oculomotor neurons. The study will highlight developing ocular motor neurons (CN3,4,6) at E11.5.

## Methods

### Bulk RNA-Seq

Bulk RNA-seq was generated through FACS sorting by utilizing a reporter gene eGFP-Isl1 in which Isl1 is differentially expressed in motor neurons comparatively to other cell types in the brainstem. Two tubes were collected for each motor neuron at E11.5, a positive and negative tube. Replicates were also established for each set of eGFP positive cells and eGFP negative cells. While the positive cells would be the motor neuron of interest expressing Isl1-eGFP the negative cells were collected for the purposes of differential gene expression. For these datasets all that was available were the FASTQ files from Illumina sequencing as well as the .fa and .gtf file for the reference genome created by the lab.  The FASTQs and reference genome were pulled down with scp/cp, the code is included as **pull_down.sh** in the code folder.

The next step is to align the reads from Bulk-Seq to a reference genome using an aligner, in this case the aligner of choice is STAR (Dobin et al, 2013) . There are other aligners as well, but STAR specifically will align to the whole genome rather than just the transcriptome as some other pseudo aligners might do. The first step in the workflow is to create a genome index based on the reference genome sequences which are provided as a fasta file and the annotated reference genome which is provided as a GTF file. This step is only necessary to run once and will then work for all interested FASTQ files for alignment.

Before alignment there is another index that must be created and this is for the transcript quantification tool of interest, RSEM (Li et al, 2011), which will provide the counts for differential gene expression. The same files are provided to RSEM to create the index. The following code to create the indices for both platforms is included as **prepare_indices.sh** in the code folder.

Once the indices are created alignment and quantification can begin. For alignment 92G of memory is allotted and 16 cores for parallelization. The FASTQ files are paired reads and such that there will be two input FASTQs to STAR. An issue with FASTQs often is that the sequences might be too long and that the files themselves need to be QC'd and trimmed, STAR, has the capability to do soft trims while doing the alignment causing this to not be a concern. To get the necessary output to get transcript quantification, --quantMode in STAR is set to TranscriptomeSAM to provide the BAM (Binary Alignment Map) that is used as input into RSEM's calculate expression function. Once alignment has occurred for a file it will then input the BAM to rsem-calculate-expression function to then generate counts for differential gene expression. This pipeline for alignment and transcript quantification can be found in the code folder as **align_and_quantify.sh**.

RSEM provides .results files which can be used as input into DESeq2 (Love et al, 2014) which is a platform to do differential gene expression. Ocular motor neurons were compared to their negative samples for regional differential expression as well as compared to each other. The top differential genes were then plotted using pheatmap (Kolde et al, 2018) in R. This was done using the code in the code folder, **Bulk.R**. The full analysis utilizes 4 scripts across the command line and R for completion.

### Single Cell RNA-Seq

The single cell datasets that are handled for this project are from E11.5 and were collected in the same FACS sorting process as bulk-seq. These are already fully aligned and ready files. The datasets were collected utilizing 10X technology through a platform known as droplet-sequencing. Droplet-sequencing occurs in a microfluidics chamber in which cells from dissociation and mRNA capture beads are released and it is expected for an oil droplet to capture a bead and a cell. The bead is then responsible for capturing the 3' end of mRNA's. This is not always a perfect scenario as at times cells can lyse before being captured releasing mRNA prematurely which can get captured in other droplets and cause contamination, a phenomenon known as ambient RNA. It is also possible for multiple cells to get captured in a droplet causing doublets. A final concern is the possibility of low-sequenced cells low quality cells. Ambient mRNA was previously handled in the data utilizing CellBender (Fleming et al, 2023), a command-line tool which estimates ambient RNAs from empty droplets (background) to then assess the contamination within cells. Doublets will be handled via an R package known as scDblFinder (Germain et al, 2022), and lastly low-quality cells are assessed through MAD (Median Absolute Deviation) filtering in which low quality cells are statistically handled(Huemos et al, 2023). To ensure a high-quality QC process both Scanpy (Wolf et al 2019) and Seurat (Hao et al, 2023) are utilized and maintained through a conda environment with the package anndata2ri created by the lab who created Scanpy (Theis et al, 2018) which allows for utilization of both packages across languages.

Cells are first assessed through an initial quality control check by visualizing summary statistics as shown in Figure 1 through the script **Pre-QC-Analysis.py**, primarily quality is assessed by total gene complexity, total counts, and mitochondrial percentage, as high mitochondrial cells might indicate stressed and dying cells not viable for analysis. Histograms for figure 1 were created by using ggplot2 (Wickham 2016) in R through the anndata2ri package conversion of Python to R. Next each sample will run through a quality control pipeline to remove doublets, low quality cells, as well as annotate cells using the package CellTypist (Dominguez et al, 2022) and utilize the developing mouse brain model to label cells, the output from this compared to the initial check is shown in Figure 1. The histogram showing total counts has cleaned up and follows a more normal distribution, the percent mitochondria is much lower, and any cells with high mitochondria, low gene complexity, and low total counts are removed. This utilizes the script **QC.py**. Once all datasets have been sufficiently filtered and annotated all samples for ocular motor neurons will be merged. Data will then be normalized to account for sequencing depth followed by variable gene selection for identifying marker genes, PCA for dimensionality reduction, then nearest neighbors and clustering will be applied to visualize cell types and cell states on low-dimensional space via UMAP. Ocular motor neurons will then be identified by using markers selected from bulk-seq for identification. The dataset will then be transferred into Seurat for differential gene expression testing using MAST (McDavid et al, 2024) to identify heterogeneity and strong candidate markers across ocular motor neurons utilizing the **DGE.R** script in the code folder.The full analysis utilizes 3 scripts from both python and R for a complete analysis.

![](/Figures_for_final/Figure_1.png)

**Figure 1: (A).** Pre-QC analysis of a single-cell sample, first plot represents a histogram of total counts per cell, the second plot shows the percentage of mitochondrial content per unbiased cluster, the third plot shows the intersection of quality control metrics in which the x axis represents total counts, the y axis represents the total number of expressed genes per cell, and the color gradient represents mitochondrial content per cell. **(B).** The same plots after QC showing tighter distributions. Histogram and inter-sectional plot were completed using ggplot2.

## Results

Aligned and quantified bulk-seq samples correlating to ocular motor neurons ran through several differential gene expression tests such as 3 vs negative, 4 vs negative, 6 vs negative, as well as each ocular motor neuron against each other. Figure 2 is an example of the differential gene expression results acquired through DESeq2 in which CN3 is tested against surrounding negative cells in the tissue in which CN3 is located. The lab has identified previously genes such as *Tmem255a* and *Palmd* which are located in the results showing that the alignment and quantification worked as successful. Another indicator of success is the clear distinction of *eGFP* in the samples against the negative samples.

![](/Figures_for_final/Figure_2.png)

**Figure 2:** Differential gene expression of CN3 vs *Isl1-eGFP* negative cells surrounding CN3

The supplementary figures show the results for all other differential gene expression tests. With the differential gene expression completed on the bulk-seq showing candidate markers for each ocular motor neuron population. It is now important to scope through further heterogeneity by utilizing single cell RNA-seq. Figure 3 shows the merged ocular motor neuron datasets with unbiased clusters as well as gene expression pulled from the bulk-RNA sequencing data. *Isl1* shows all motor neurons as well as identifying negative cells, here it would be expected to only have 3 motor neurons, since 6 is close to 7, 7 is also in the sample. This is why single cell can be critical as being able to visualize heterogeneity is a key facet into data exploration in transcriptomics. *Tmem255a*, *Vip*, *Mical2*, correlate respectively to 3,4,6. An example of the positives of having bulk data is due to the increase in sequencing depth compared to single cell. *Mical2* shows clear differential expression in CN6 compared to the other ocular motor neurons but is hidden due to sequencing depth. Both approaches can significantly complement one another in terms of exploration. It can also be noted that *Vip* and *Tmem255a* are shown to be differentially expressed even within their respective motor neuron types showing motor neuron heterogeneity.

![](/Figures_for_final/Figure_3.png)

**Figure 4: (A).** Unbiased Clustering of CN3/4/6 samples, **(B).** FeaturePlot of *Isl1* to show motor neurons, **(C).** FeaturePlot of *Tmem255a* to show CN3, **(D).** FeaturePlot of *Vip* to show CN4, **(E).** FeaturePlot of *Mical2* to show CN6.

To probe differential gene expression and explore motor neuron heterogeneity the single cell merged dataset is run on the MAST protocol to identify genes significant per cluster. The resulting dot plot is shown in Figure 4. The comparison is made on the unbiased Leiden clusters to show potential heterogeneity amongst the motor neurons themselves as motor neurons such as CN3 innervate multiple ocular muscles. As shown motor neurons such as CN3 which are divided into cluster ID 18,8,4 have unique expression showing heterogeneity at the same level.

![](/Figures_for_final/Figure_4.png)

**Figure 5**: DotPlot showing top 3 differentially expressed genes per unbiased clusters that are noted in Figure 4A

## Discussion

The usage of both bulk and single cell RNA-sequencing provide insights into the transcriptomic landscape of developing ocular motor neurons. Bulk provides low resolution (heterogeneity) but higher sequencing depth and deeper insight into total gene expression while single cell provides higher resolution but lower sequencing depth. An example can be seen from Figure 3 in which the *Mical2* which is a gene picked up on the bulk-seq as a top 50 differentially expressed gene is not as strongly picked up on the single-cell compared to what the bulk has suggested. The single cell however provides a layer of heterogeneity that cannot be captured in bulk and that can be assessed in Figure 3 and Figure 4 as the UMAP suggests that there are multiple clusters per motor neuron such as in CN3 which is composed of 3 clusters 18,8, and 4 which in Figure 4 shows unique gene expression per cluster. This aligns with what is known with motor neurons in which CN3 particularly innervates 5 ocular muscles, it would be expected if one were to look at E18.5 rather than a younger E11.5 to potentially see 5 populations.

By assessing the full transcriptomic profile of these motor neurons, reporter mice can be purchased to study certain congenital disorders related to motor neurons. To get a strong sense of what genes could be targets, the lab has done in-situ studies to validate some of these markers such as the *Tmem255a* which is also seen in Figure 3c. There is still much work to build a list of potential markers for each subpopulation of a motor neuron. Seeing differential gene expression picked up in both bulk/scRNA-seq provides higher confidence in targets that can be used first for in-situ hybridization experiments for confirmation of localized gene expression and then as a reporter mice for further study.

A first step in providing further confidence is to assess if these genes are shown on the Allen Brain Atlas (Allen Institute for Brain Science, 2004) in their in-situ studies. *Isl1* gene expression can be used as the reference as shown in Figure 5a and compared to a potential target gene of choice, *Sema5a*, show in Figure 5b which was pulled out of the DGE results showing differential gene expression in CN4, part of cluster 17 compared to all other clusters. As shown the boxed region in the *Isl1*/*Sema5a* is CN4.

![](/Figures_for_final/Figure_5.png)

**Figure 5: (A).** *Isl1* stain on Allen Brain Atlas website, red box correlates to physical location of CN4 via *Isl1* gene expression, **(B).** *Sema5a* stain on Allen Brain Atlas website, red box correlates to a region of CN4 in which *Sema5a* can be seen expressed, images from <http://atlas.brain-map.org>.

This shows how the single-cell data highlights that gene expression is heterogeneous even within the same celltype, Figure 6, shows the expression of *Sema5a* which is not localized through all of CN4 itself. For further studies to identify subpopulations it would be best to look at the most mature age in the data which is E18.5 and perform the same analysis as shown. However, it is interesting to identify the early transcriptomic profile that constitutes these ocular motor neurons.

![](/Figures_for_final/Figure_6.png)

**Figure 6:** FeaturePlot of *Sema5a* with red box correlating to cluster of CN4

## References:

Allen Institute for Brain Science. 2011. Allen Reference Atlas -- Mouse Brain [brain atlas]. Available from <http://atlas.brain-map.org>.

Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29: 15--21. <doi:10.1093/bioinformatics/bts635>.

Dominguez Conde C, et al. 2022. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376: eabl5197.

Germain P, Lun A, Garcia Meixide C, Macnair W, Robinson M. 2022. Doublet identification in single-cell sequencing data using scDblFinder. F1000Research. <doi:10.12688/f1000research.73600.2>.

Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, Srivastava A, Molla G, Madad S, Fernandez-Granda C, Satija R. 2023. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol. <doi:10.1038/s41587-023-01767-y>.

Heumos L, Schaar AC, Lance C, et al. 2023. Best practices for single-cell analysis across modalities. Nat Rev Genet. <doi:10.1038/s41576-023-00586-w>. Kolde R. 2018. pheatmap: Pretty Heatmaps. R package version 1.0.

Li B, Dewey CN. 2011. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12: 323. <doi:10.1186/1471-2105-12-323>.

Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15: 550. <doi:10.1186/s13059-014-0550-8>.

McDavid A, Finak G, Yajima M. 2024. MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.32.0. <https://github.com/RGLab/MAST/>.

Fleming SJ, Chaffin MD, Arduini A, Akkad AD, Banks E, Marioni JC, Phillipakis AA, Ellinor PT, Babadi M. 2023. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. Nat Methods. <doi:10.1038/s41592-023-01943-7>.

Wickham H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19: 15. <doi:10.1186/s13059-017-1382-0>.

## Supplementary Plots

### CN4 vs CN4_eGFP_negative

![](/Figures_for_final/4_v_neg.png)

### CN6 vs CN6_eGFP_negative

![](/Figures_for_final/Figures/6_v_neg.png)

### CN4 vs CN3

![](/Figures_for_final/3v4.png)

### CN3 vs CN6

![](/Figures_for_final/3v6.png)

### CN4 vs CN6

![](/Figures_for_final/4v6.png)
