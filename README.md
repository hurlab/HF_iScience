# HF_iScience
Scripts for "Single-cell RNA sequencing identifies hippocampal microglial dysregulation in diet-induced obesity"
## Description
Data were collected after mapping all fastq files to the mm10 genome with the cellranger version 4.0.0 software with default paramaters by using "cellranger count".

## Required packages:
Seurat v3, ggplot2 , richR, DEseq2, CellChat v1.0.0, velocyto.R 

## Script details
1. Read data, integration, clustering

run the run_main.r for clustering, cell annotation, differential expressed analysis and functional enrichment analysis. 

2. cell cell communication analysis

run the cell_cell.r for the cell cell communication analysis analysis and figures. 

3. RNA velocity

run the velocity.r for all RNA velocity analysis and figure generation. 
