# Analysis Proposal for EE282 Final Project

Author: Nima Shirooni

## Introduction 
Motor neurons are key celltypes in the brain in which they are nerve cells responsible 
for transmission of signals from the brain and spinal cord to target muscles. The developmental
understanding of motor neurons is still nascent. Particularly, spinal cord motor neurons are 
more well established in the literature compared to brainstem motor neurons. This provides a 
key motivation into the investigation of the development of brainstem motor neurons. Not 
only are they developmentally interesting from a biological perspective but from a medical 
perspective, inappropriate development can lead to rare congential neurologic diseases.
Since much are not known about these brainstem motor neurons, single cell RNA-seq can provide
interesting insight into the heterogeneity of these cells. The project has led to the efforts of
providing various key datasets that will be used for motor neuron exploration. These datasets are
bulk RNA-sequencing, single cell/nuclei RNA-sequencing, and spatial RNA-sequencing of the developing 
mouse embryo brainstems from embryonic days E9.5-E18.5. Data was collected by my PI, Matthew Rose, during
his post-doctoral studies and is now an integral set of data for the lab to analyze. These datasets
pose a vast and fascinating exploration of brainstem developmental data that can guide future experiments.

## Goals 
The large scope of the project is to take all of the existing datasets and create a transcriptomic 
atlas of motor neuron development across embryonic ages E9.5-E18.5, while this is the big picture,
a smaller goal is to learn more about 3 particular motor neurons involved in eye movement, CN3,4, and 6
which are all involved in ocular movement. Learning deeper into these motor neurons would provide answers such
as potential subpopulations, differentially expressed markers (genes) in these subpopulations, and 
differentially expressed markers between each motor neuron. These goals are particularly suited for 
single cell/nuclei datasets. Bulk-sequencing was also collected at age E11.5, while single-cell RNA
sequencing provides cellular hetereogeneity it suffers in sequencing depth, genes we identify as differentially 
expressed but rather low sequencing depth, Bulk data can potentially bridge the sequencing depth gap and be used to
support potentially "hidden" genes in the single-cell data. Additionally spatial data using a technique called 
Slide-Sequencing was obtained at E11.5, an interesting analyses that can be used with this data is to take the 
single cell data and then predict the spatial locations on the slide-sequencing data, this can help confirm cellular 
hetereogenity of identified motor neuron subpopulations as well as show the validity of celltype annotations in the
single cell data for all other brainstem celltypes. Given these goals, for the project, I propose a full analysis 
of CN3,4,6 and integrating bulk and spatial RNA-seq to support my analyses. 

# Analyses

