Transcriptional Analysis of the Correlation Between Traumatic Brain Injury and Dementia: Investigating the Impact of APOE Genotype

Overview:

This repository contains the datasets and code required to replicate the figures in the 20.440 final project write up by Sarah Sobol and Lauren Ramlan. At a high level, the repository contains the necessary datasets and code to perform a transcriptional analysis of the correlation between traumatic brain injury (TBI) and demenita with a focus on investigating the impact of APOE4 genotype. The data for this analysis comes from the Aging, Dementia and Traumatic Brain Injury Study (Miller JA, et al., 2017), which was a detailed characterization of brain samples from control and TBI exposed donors for an aged population-based cohort from the Adult Changes in Though (ACT) study. The datasets are contained in the Data folder, the code required to process the data and generate figures is contained in the Code folder, and figures are contained in the Figures folder of this repository, and excel files of the results from GO term enrichment and DAVID Functional Annotation analysis are contained in the 
. As an overview, the respository contains the necessary datasets and code to produce volcano plots comparing gene expression in donors with both TBI and dementia to donors with TBI and no dementia as well as between donors with TBI and dementia to donors without TBI and with dementia. These volcano plots can be generated for samples for all brain regions combined as well as each of four brain regions (hippocampus, temportal cortex, parietal cortex, and forebrain white matter) stratified by whether or not the sample came from a donor who had the APOE4 allele. Furthermore, the repository contains code to generate venn diagrams looking at how many genes overlap between the TBI + dementia vs. TBI + no dementia and TBI + dementia vs. no TBI + dementia comparisons for each stratification. of brain region and APOE genotype. The repo

produce a volcano plot of the differentially expressed genes between aged donors that either are healthy controls (have not experienced a traumatic brain injury and do not have dementia) or donors that have experienced a traumatic brain injury and have dementia.  The goal of this analysis and figure generation is to determine which genes are differentially expressed between the two groups of donors and whether those genes are upregulated or downregulated.  The RNA-sequencing datafrom this study was used to compare gene expression between two of the groups of patients (control and those with both TBI and dementia). This repository takes the FPKM values obtained for the genes analyzed on each of the patient samples and first creates a subset of this dataframe that only contains genes that had a FPKM value greater than 0 for at least 50% of the patient samples. Then this dataframe transposes the data frame with rows of genes and columns of samples to be rows of samples and columns of genes and then merges this dataframe of gene expression for each sample with the diagnostic group of that patient sample (No TBI + No Dementia, No TBI + Dementia, TBI + No Dementia, or TBI + Dementia). Using these dataframes, the code performs a limma-trend analysis (Ritchie ME, et al., 2015) with the contrasting groups being TBI + Dementia (Y_Dementia) and No TBI + No Dementia (N_No_Dementia). Using the adjusted p values for each of the genes (corrected with Benjamini-Hochberg procedure) and the log fold change the genes are plotted on a volcano plot to determine which genes are upregulated, downregulated, or not significant. The code also stores the gene_id numbers of the significantly differentially expressed genes and obtains the gene_symbol associated with that gene_id and labels the datapoint on the plot.

Data:

All of the datasets contained in the data folder are from the Aging, Dementia and TBI Study (Miller JA, et al. 2017). These datasets contain information for 107 donors and over 50,000 genes. To generate the dataset, samples were collected fromthe hippocampus, temporal cortex, and parietal cortex (both grey and white matter) in 55 aged donors with TBI and their matched controls (a total of 107 donors after quality control) and RNA-sequencing was performed on the 377 samples. The sequencing results were then aligned and aggregated at the gene levels using the RSEM algorithm resulting in fpkm values. Furthermore, information about each of the donors such as whether they experienced a traumatic brain injury, dementia status, race, and age is contained in the datasets.

fpkm_table_normalized.csv

This dataset contains a (row, column) matrix of fpkm values obtained for each (gene, sample) from the RSEM analysis of the RNA sequencing data. The fpkm values in this dataset were normalized by correcting for RNA quality (RIN) and batch effects and scaling such that the log2(fpkm) across the entire dataset remained unchanged after normalization. The first row of the dataset contain the RNA-sequencing profiles of the samples and the first column contains the unique identifiers for each of the genes.

columns-samples.csv

This dataset contains additional information about the samples that were profiled with RNA sequencing such as the donor identification number that corresponds with the sample and the brain region that sample came from.

rows-genes.csv

This dataset contains that information about the genes measured by RNA-sequencing for each of the samples such as the gene symbol and name in addition to the unique gene identification number.

donor_information.csv

This dataset contains additional information about the donors that the samples came from. With this dataset the donor identification number is linked to information such as the age of the donor, the donor's race, whether the donor experienced a traumatic brain injury, whether the donor has the apo_e4 allele, and whether the donor has dementia.

Folder Structure:

There are 3 folders in this repository: Data, Code, and Figures

The Data folder contained the 4 csv files that will be needed to run the code to generate the figure
- fpkm_table_normalized.csv
- columns-samples.csv
- rows-genes.csv
- donor_information.csv

The Code folder contains a singular .Rmd file entitled Differential_Gene_Expression.Rmd which will process the datasets into the necessary format and use those dataframes to perform a limma trend analysis (Ritchie ME, et al., 2015) and generate a volcano plot that shows differential gene expression between the control group (no TBI + no dementia) and the group with both TBI and dementia. This code will also store the subset of genes that are differentially expressed as well as whether they are upregulated or down regulated and include the gene information associated with that gene_id.

The Figures folder contains a singular .png file that is a volcano plot demonstrating that genes that are differentially expressed between donors with traumatic brain injury and dementia and donors that are healthy controls. The genes that are upregulated in the tbi + dementia donors relative to the healthy control donors are in clue and the downregulated genes are in red with the not significant genes in grey. The differentially expressed genes are also labeled with the gene symbol at the data point.

Installation:

To run the code, R and RStudio need to be downloaded. Then the .Rmd file in the code folder and all of the datasets need to be downloaded and stored in the same file on ones computer. An R project needs to be created in the folder that the .Rmd file and the datasets are stored in. To run the code one also needs to install the following packages:
- limma (one may have to install BioCManager to install limma if not already installed)
- ggplot2
- dplyr
Once the packages are installed and if the datasets are in the file path with the R project and .Rmd file the code should run from there to generate the volcano plot.

Citations:

2016 Allen Institute for Brain Science. Aging, Dementia and TBI study. Available fromL http://aging.brain-map.org/
Miller JA et al. (2017) Neuropathological and transcriptomic characteristics of the aged brain, Elife doi: 10.7554/eLife.31126

Ritchie ME, et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research doi: 10.1093/nar/gkv007
