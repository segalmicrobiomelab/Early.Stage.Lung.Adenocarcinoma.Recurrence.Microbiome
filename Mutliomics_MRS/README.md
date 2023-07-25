# Multiomics and Microbial Risk Score Analysis 

# Required data 
1. Microbe Abundance (table.filtered.qza):An abundance table with microbial species in rows and samples in columns.
2. Taxonomy Table (taxonomy.filtered.qza): Table of the taxonomic annotation
3. MetaData (Surgical.Cohort.Map.txt): A spreadsheet with samples in rows and metadata in columns
4. Tree File (rooted-tree_quality.filtered.qza)
5. RNAseq Raw Abundance (NovaSeq_Raw_Counts_Gene_Symbols):An abundance table with genes in rows and samples in columns

# COX_PH analysis 
open R 
Follow code Analysis1_Cox_PH.R

# MRS analysis 
Follow code MRS script.R

# Multiomics analysis 
Follow code Multi_omics_analysis.R
