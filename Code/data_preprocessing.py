# -*- coding: utf-8 -*-
"""Data_Preprocessing.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1SM5x73GF17C390haTZ5QGMHMdyXwJPpA
"""

# Import necessary packages
import pandas as pd
import numpy as np

# Download data
# Insert your own file path where the variable file_path is
fpkm_table_normalized = pd.read_csv('/file_path/fpkm_table_normalized.csv')
donor_info = pd.read_csv('/file_path/donor_information.csv')
columns_samples = pd.read_csv('/file_path/columns-samples.csv')

# Format data table with gene_id as the columns and rnaseq_profile_id as the row index
fpkm_table_normalized.set_index('gene_id \ rnaseq_profile_id', inplace = True, drop = True)
fpkm_table_normalized.index.name = None

# Transpose the dataframe so rnaseq_profile_id is the rows and gene_id is the columns
fpkm_transpose = fpkm_table_normalized.T
# Make rnaseq_profile_id a column
fpkm_transpose.reset_index(inplace = True)
fpkm_transpose.rename(columns = {'index': 'rnaseq_profile_id'}, inplace = True)

# Merge the information about a samples rnaseq_profile_id with whether the donor
# that sample is from had a TBI or their dementia diagnosis by the donor_id
sample_diagnosis_info = pd.merge(columns_samples[['rnaseq_profile_id', 'donor_id',
                                                  'structure_acronym']],
                                 donor_info[['donor_id', 'apo_e4_allele',
                                             'ever_tbi_w_loc', 'act_demented']],
                                  on = 'donor_id')
sample_diagnosis_info.dropna(inplace = True)
# Remove the donor_id column
sample_diagnosis_info.drop(columns = ['donor_id'], inplace = True)
# Create a new variable with the group the patient is in by combining the
# information about ever having a TBI with their dementia diagnosis
sample_diagnosis_info['group'] = sample_diagnosis_info['apo_e4_allele'] + '_' + sample_diagnosis_info['ever_tbi_w_loc'] + '_' + sample_diagnosis_info['act_demented']
# Remove all columns except the rnaseq_profile_id and the group
sample_diagnosis_grouped = sample_diagnosis_info[['rnaseq_profile_id',
                                                  'structure_acronym', 'group']]
# Change all spaces in the values for group to underscores
sample_diagnosis_grouped['group'] = sample_diagnosis_grouped['group'].str.replace(" ", "_")
# Convert rnaseq_profile_id to a string to match the type of rnaseq_profile id
# in the fpkm_filtered_transpose dataframe
sample_diagnosis_grouped['rnaseq_profile_id'] = sample_diagnosis_grouped['rnaseq_profile_id'].astype(str)
# Merge the group information with the fpkm_filtered_transpose dataframe by the rnaseq_profile_id
samples_grouped = pd.merge(fpkm_transpose, sample_diagnosis_grouped, on = 'rnaseq_profile_id')
# Set the rnaseq_profile_id as the row index
samples_grouped.set_index('rnaseq_profile_id', inplace = True, drop = True)
samples_grouped.index.name = None

# Create separate dataframes subsetted by each brain region and drop the
# strucutre_acronym column
samples_grouped_HIP = samples_grouped[samples_grouped['structure_acronym'] == 'HIP']
samples_grouped_HIP.drop(columns = ['structure_acronym'], inplace = True)
samples_grouped_TCx = samples_grouped[samples_grouped['structure_acronym'] == 'TCx']
samples_grouped_TCx.drop(columns = ['structure_acronym'], inplace = True)
samples_grouped_PCx = samples_grouped[samples_grouped['structure_acronym'] == 'PCx']
samples_grouped_PCx.drop(columns = ['structure_acronym'], inplace = True)
samples_grouped_FWM = samples_grouped[samples_grouped['structure_acronym'] == 'FWM']
samples_grouped_FWM.drop(columns = ['structure_acronym'], inplace = True)

# Calculate the percentage of non-zero values in each column for each brain region
non_zero_percentage_HIP = (samples_grouped_HIP.astype(bool).sum() / len(samples_grouped_HIP)) * 100
selected_columns_HIP = non_zero_percentage_HIP[non_zero_percentage_HIP > 50]
column_names_HIP = selected_columns_HIP.index.tolist()

non_zero_percentage_TCx = (samples_grouped_TCx.astype(bool).sum() / len(samples_grouped_TCx)) * 100
selected_columns_TCx = non_zero_percentage_TCx[non_zero_percentage_TCx > 50]
column_names_TCx = selected_columns_TCx.index.tolist()

non_zero_percentage_PCx = (samples_grouped_PCx.astype(bool).sum() / len(samples_grouped_PCx)) * 100
selected_columns_PCx = non_zero_percentage_PCx[non_zero_percentage_PCx > 50]
column_names_PCx = selected_columns_PCx.index.tolist()

non_zero_percentage_FWM = (samples_grouped_FWM.astype(bool).sum() / len(samples_grouped_FWM)) * 100
selected_columns_FWM = non_zero_percentage_FWM[non_zero_percentage_FWM > 50]
column_names_FWM = selected_columns_FWM.index.tolist()

# Find the list of genes that have over 50% of the samples as non-zero values in
# all brain regions
common_genes = set(column_names_HIP) & set(column_names_TCx) & set(column_names_PCx) & set(column_names_FWM)

# Create new filtered dataframes containing only the commonly shared genes that
# have over 50% of samples as non-zero values in all brain regions for each of
# the brain regions
samples_grouped_HIP_filtered = samples_grouped_HIP[list(common_genes)]
samples_grouped_TCx_filtered = samples_grouped_TCx[list(common_genes)]
samples_grouped_PCx_filtered = samples_grouped_PCx[list(common_genes)]
samples_grouped_FWM_filtered = samples_grouped_FWM[list(common_genes)]
# Concatenate the individual brain region dataframes to create a combined brain
# region dataframe
samples_grouped_filtered = pd.concat([samples_grouped_HIP_filtered,
                                      samples_grouped_TCx_filtered,
                                      samples_grouped_PCx_filtered,
                                      samples_grouped_FWM_filtered])

# Transpose the filtered dataframes for the samples_grouped_filtered for each
# brain region and drop the group column to make a dataframe with rows of genes
# and columns of samples ids
fpkm_HIP_trans = samples_grouped_HIP_filtered.drop(columns = ['group'])
fpkm_HIP = fpkm_HIP_trans.T
fpkm_TCx_trans = samples_grouped_TCx_filtered.drop(columns = ['group'])
fpkm_TCx = fpkm_TCx_trans.T
fpkm_PCx_trans = samples_grouped_PCx_filtered.drop(columns = ['group'])
fpkm_PCx = fpkm_PCx_trans.T
fpkm_FWM_trans = samples_grouped_FWM_filtered.drop(columns = ['group'])
fpkm_FWM = fpkm_FWM_trans.T
fpkm_combined_trans = samples_grouped_filtered.drop(columns = ['group'])
fpkm_combined = fpkm_combined_trans.T

# Save processed dataframes for analysis in R
fpkm_HIP.to_csv('/content/drive/My Drive/440 project/Processed_Data/fpkm_HIP.csv')
fpkm_TCx.to_csv('/content/drive/My Drive/440 project/Processed_Data/fpkm_TCx.csv')
fpkm_PCx.to_csv('/content/drive/My Drive/440 project/Processed_Data/fpkm_PCx.csv')
fpkm_FWM.to_csv('/content/drive/My Drive/440 project/Processed_Data/fpkm_FWM.csv')
fpkm_combined.to_csv('/content/drive/My Drive/440 project/Processed_Data/fpkm_combined.csv')
samples_grouped_HIP_filtered.to_csv('/content/drive/My Drive/440 project/Processed_Data/samples_grouped_HIP.csv')
samples_grouped_TCx_filtered.to_csv('/content/drive/My Drive/440 project/Processed_Data/samples_grouped_TCx.csv')
samples_grouped_PCx_filtered.to_csv('/content/drive/My Drive/440 project/Processed_Data/samples_grouped_PCx.csv')
samples_grouped_FWM_filtered.to_csv('/content/drive/My Drive/440 project/Processed_Datasamples_grouped_FWM.csv')
samples_grouped_filtered.to_csv('/content/drive/My Drive/440 project/Processed_Data/samples_grouped_combined.csv')