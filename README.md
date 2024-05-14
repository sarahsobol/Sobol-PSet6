**Transcriptional Analysis of the Correlation Between Traumatic Brain Injury and Dementia: Investigating the Impact of APOE Genotype**

**Overview:**

This repository contains the datasets and code required to replicate the figures in the 20.440 final project write up by Sarah Sobol and Lauren Ramlan. At a high level, the repository contains the necessary datasets and code to perform a transcriptional analysis of the correlation between traumatic brain injury (TBI) and demenita with a focus on investigating the impact of APOE4 genotype. The data for this analysis comes from the Aging, Dementia and Traumatic Brain Injury Study (Miller JA, et al., 2017), which was a detailed characterization of brain samples from control and TBI exposed donors for an aged population-based cohort from the Adult Changes in Though (ACT) study. The datasets are contained in the Data folder, the code required to process the data and generate figures is contained in the Code folder, and figures are contained in the Figures folder of this repository, and excel files of the results from GO term enrichment and DAVID Functional Annotation analysis are contained in the Gene_ontology_results folder. As an overview, the respository contains the necessary datasets and code to produce volcano plots comparing gene expression in donors with both TBI and dementia to donors with TBI and no dementia as well as between donors with TBI and dementia to donors without TBI and with dementia. These volcano plots can be generated for samples for all brain regions combined as well as each of four brain regions (hippocampus, temportal cortex, parietal cortex, and forebrain white matter) stratified by whether or not the sample came from a donor who had the APOE4 allele. These volcano plots are generated based on a differential expression analysis using the limma package in R. Furthermore, the repository contains code to generate venn diagrams looking at how many genes overlap between the TBI + dementia vs. TBI + no dementia and TBI + dementia vs. no TBI + dementia comparisons for each stratification. of brain region and APOE genotype. The repository also contains the code necessary to create chord diagrams to visually link identified differentially expressed genes with pathway annotations. Finally, the repository contains the necessary datasets and code to perform principal component analysis (PCA) using dataframes of RNA-sequencing results filtered to contain only genes identified as differentially expressed in the TBI + dementia group of donors.

**Data:**

All of the datasets contained in the data folder are from the Aging, Dementia and TBI Study (Miller JA, et al. 2017). These datasets contain information for 107 donors and over 50,000 genes. To generate the dataset, samples were collected from the hippocampus, temporal cortex, parietal cortex, and forebrain white matter in 55 aged donors with TBI and their matched controls (a total of 107 donors after quality control) and RNA-sequencing was performed on the 377 samples. The sequencing results were then aligned and aggregated at the gene levels using the RSEM algorithm resulting in FPKM values. Furthermore, information about each of the donors such as whether they experienced a traumatic brain injury, dementia status, race, and age is contained in the datasets.

fpkm_table_normalized.csv

This dataset contains a (row, column) matrix of fpkm values obtained for each (gene, sample) from the RSEM analysis of the RNA sequencing data. The fpkm values in this dataset were normalized by correcting for RNA quality (RIN) and batch effects and scaling such that the log2(fpkm) across the entire dataset remained unchanged after normalization. The first row of the dataset contain the RNA-sequencing profiles of the samples and the first column contains the unique identifiers for each of the genes.

columns-samples.csv

This dataset contains additional information about the samples that were profiled with RNA sequencing such as the donor identification number that corresponds with the sample and the brain region that sample came from.

rows-genes.csv

This dataset contains that information about the genes measured by RNA-sequencing for each of the samples such as the gene symbol and name in addition to the unique gene identification number.

donor_information.csv

This dataset contains additional information about the donors that the samples came from. With this dataset the donor identification number is linked to information such as the age of the donor, the donor's race, whether the donor experienced a traumatic brain injury, whether the donor has the apo_e4 allele, and whether the donor has dementia.

**Folder Structure:**

There are 4 folders in this repository: Data, Code, Figures, and Gene_ontology_results


The Data folder contains 4 csv files that will be needed to execute the code and create figures and lists of differentially expressed genes. This data can all be downloaded from the Allen Institute for Brain Science on the Aging, Dementia, and TBI Study website under the download tab.
- fpkm_table_normalized.csv
- columns-samples.csv
- rows-genes.csv
- donor_information.csv


The Code folder contains 4 different scripts to generate both outputs of lists of differentially expressed genes and figures.
- data_preprocessing.py
  
This scipt runs in python and uses the fpkm_table_normalizes.csv, donor_information.csv, and columns-samples.csv as inputs. It takes the FPKM values and merges them with the information from the donor_information and columns-samples datasets to create two separate datasets: one with the FPKM values only where the rows are genes and the columns are sample IDs and one with the FPKM values and information about which "group" the sample belongs to where the rows are sample IDs and the columns are genes. The groups indicate the APOE genotype, TBI history, and dementia status of the donor the sample comes from. This script also subsets the dataframe by the brain region the sample is from. Furthermore, the script filters the datasets to only include genes where over 50% of the FPKM values for the samples in each brain region are non-zero. It also removes patients who have N/A as their apoe4 allele status. The script returns two dataframe types for each brain region as well as all brain regions combined: FPKM values with genes as rows and sample IDs as columns and FPKM values + group with sample IDs as rows and genes as columns.
  
- Differential_Gene_Expression_Analysis.Rmd
  
This script runs in R and uses the outputs of the data_preprocessing.py script as well as the rows-genes.csv dataframe as inputs. It uses the packages limma, ggplot2, dplyr, and VennDiagram to perform a limma differential gene expression analysis for the contrasts of TBI + Dementia vs. TBI + No Dementia and TBI + Dementia vs. No TBI + Dementia for both samples from the APOE4 positive and negative sample groups. It uses this differential gene expression analysis to store the upregulated, downregulated, and nominally significant genes for each contrast as well as the genes that a commonly upregulated, downregulated, and nominally significant between the two types of contrasts for both the APOE4 positive and negative populations. It also uses these differentially expressed genes to plot volcano plots for each contrast and venn diagrams to show the overlap in genes between the contrasts stratified by brain region and APOE genotype. This script also outputs two data files: apoe4_all_df.csv and no_apoe4_all_df.csv. These data files contain the FPKM values with the group for all samples (sample ID as rows and gene as columns) subset so that only genes identified as commonly differentially expressed are included for the APOE4+ and APOE4- analyses respectively.

- Gene_Ontology.Rmd
  
This script runs in R and uses the inputs of apoe4_pathway.csv and non_apoe4_pathway.csv, which are data files that contain information from the DAVID Functional Annotation Tool where genes are paired with the pathway they are annotated for. The code uses the library circlize to plot a chord diagram connecting genes with the pathways they are annoted with.
  
- pca_analysis.py
  
This script runs in python and uses the inputs of apoe4_all_df and no_apoe4_all_df, which are outputs of the Differential_Gene_Expression_Analysis.Rmd code script. It uses these inputs to perform a PCA analysis and plot scatterplots and bar graphs of variance explained by the number of principle components. The PCA analysis is performed for the subset of patients with traumatic brain injury separated by whether they have the APOE4 allele or not. The PCA is performed for the TBI+ and APOE4 positive or negative donor groups using each of the data sets. The Mahalanobis distance is then calculated for the separation of the dementia and no dementia groups in each of the PCA scatterplots.


The Figures folder contains 7 folders used in the final project write up.
- Figure1-Differentially_Expressed_Genes.png
  
This figure contains volcanon plots showing the upregulated, downregulated, and non-significant genes for either dementia among all donors with TBI or TBI among all donors with dementia stratified by APOE genotype as determined by limma analysis for the parietal cortex samples. It also contains the number of differentially expressed genes in each brain region for dementia from the TBI + dementia vs. TBI + no dementia or TBI from the TBI + dementia vs. No TBI + dementia comparisons by APOE genotype. The overlapping number in the middle of the venn diagrams represents the quantity of overlapping differentially expressed genes between dementia and TBI.

- Figure2-Chordplots_of_enriched_pathways.png
  
This figure shows pathways enriched from the differentially expressed genes for TBI and dementia in either the APOE4+ or APOE4- populations where the laetters represent the molecular function of the enriched pathway and the enriched genes are labeled on the right.

- Figure3-PCA_analysis.png
  
This figure contains plots where the X and Y axes are the top principal components identified by either the APOE4+ (left) or APOE4- (right) stratification and DEG analysis. The PCAs contain plots of either the APOE4+ donors with TBI (top) or APOE4- donors with TBI (bottom). All plots are colored by ground truth dementia diagnosis with pruple representing dementia and orange representing no dementia.

- Supplementary_Figure1.png
  
This figure shows the amount of variance explained by the number of principal components from PCA analysis of DEGs identified for each genotype by donor subpopulation.

- Supplementary_Figure2.png
  
This figure contains volcanon plots showing the upregulated, downregulated, and non-significant genes for either dementia among all donors with TBI or TBI among all donors with dementia stratified by APOE genotype as determined by limma analysis for all brain region samples. It also contains the number of differentially expressed genes in each brain region for dementia from the TBI + dementia vs. TBI + no dementia or TBI from the TBI + dementia vs. No TBI + dementia comparisons by APOE genotype. The overlapping number in the middle of the venn diagrams represents the quantity of overlapping differentially expressed genes between dementia and TBI.

- Supplementary_Figure3.png
  
This figure shows the significantly enriched gene ontology terms for biological processes in APOE4+ donor DEG analysis from the parietal cortex.
  
- Supplementary_Figure4.png
  
This figure shows the significantly enriched gene ontology terms for biological processes in APOE4- donor DEG analysis from the parietal cortex.


The Gene_ontology_results folder contains excel files that have information from the GO term enrichment and DAVID Functional Annotation.
- GO_biological_processes.xlsx
  
This file contains the significantly enriched gene ontology terms as determined by PANTHER for both the APOE4+ and APOE4- analysis identified differentially expressed genes with the FDR corrected p-value and -log10(FDR corrected p-value).
  
- apoe4_pathway.csv

This file contains the pathway annotations and associated genes from the DAVID Functional Annotation Tool for the differentially expressed genes identified by the APOE4+ analysis.
  
- non_apoe4_pathway.csv

This file contains the pathway annotations and associated genes from the DAVID Functional Annotation Tool for the differentially expressed genes identified by the APOE4- analysis.

**Installation:**

To run the code, python (Jupyter Notebook optional), R, and RStudio need to be downloaded. First all of the datasets in the data folder need to be downloaded. Then the data_preprocessing.py file should be executed. The packages needed to run this scipt include pandas and numpy. The output of this code should then be downloaded and stored with the other data files. An R project then needs to be created in the same computer pathway as the two .Rmd files and the datasets are stored in. The Differential_Gene_Expression_Analysis.Rmd should be run next using the outputs of the data_preprocessing.py code as well as the rows-genes.csv. To run this .Rmd one also needs to install the following packages:
- limma (one may have to install BioCManager to install limma if not already installed)
- ggplot2
- dplyr
- -VennDiagram
Once the packages are installed and if the datasets are in the file path with the R project and Differential_Gene_Expression_Analysis.Rmd file the code should run from there to generate the lists of differentially expressed genes, volcano plots, venn diagrams, and datasets to use in the PCA analysis. The two data files from the Differential_Gene_Expression_Analysis.Rmd can then be input into the pca_analysis.py file along with loading the following packages: pandas, numpy, scipy, stats from scipy, multiple tests from statsmodels.stats.multitest, StandardScaler from sklearn.preprocessing, PCA from sklearn.decomposition, seaborn, and matplotlib.pyplot. Finally, the apoe4_pathway.csv and non_apoe4_pathway.csv files can be downloaded and stored in the same file pathway as the R project to run the Gene_Ontology.Rmd code to generate the chord diagrams. The only package needed for this code to run is circlize.

**Citations:**

2016 Allen Institute for Brain Science. Aging, Dementia and TBI study. Available fromL http://aging.brain-map.org/
Miller JA et al. (2017) Neuropathological and transcriptomic characteristics of the aged brain, Elife doi: 10.7554/eLife.31126
Ritchie ME, et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research doi: 10.1093/nar/gkv007
