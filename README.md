# Oral_Cancer_Psuedotime_Analysis

##Project Overview

The project focuses on processing and analyzing a single-cell RNA sequencing dataset to explore the cellular heterogeneity of oral cancer. The analysis includes filtering cells, identifying clusters, removing batch effects, and inferring cell trajectories using pseudotime analysis.

##Tools and Libraries
The following R packages were used in this project:
Seurat: For single-cell RNA-seq data processing, normalization, clustering, and visualization.
Monocle3: For trajectory inference and pseudotime analysis.
dplyr: For data manipulation.
ggplot2: For data visualization.
harmony: For batch effect correction.


##Data Description
The dataset used in this analysis is an RDS file taken from cellxgene census, which contains single-cell RNA sequencing data from oral cancer samples. The metadata includes cell type, development stage, and batch information.

Analysis Workflow
1. Data Preprocessing
Loaded the oral cancer dataset from an RDS file.
Created a Seurat object using the expression matrix to extract feature and expression counts.
Added mitochondrial content percentage to the metadata.

2. Quality Control
Visualized RNA features, counts, and mitochondrial content to identify optimal filtering cutoffs.
Applied filters to retain cells with a reasonable number of features and counts.

3. Clustering and Visualization
Identified and visualized distinct cell clusters.
Used UMAP for dimensionality reduction and visualization of batch effects.

4. Batch Effect Correction
Applied the Harmony algorithm to remove batch effects based on library preparation batch information.
Re-clustered the data after batch effect correction and visualized the updated clusters.

5. Trajectory Inference
Converted the Seurat object to a Cell Data Set (CDS) for Monocle3 analysis.
Assigned cluster and UMAP coordinates from Seurat to the CDS.
Learned the cell trajectory graph and set the root cells based on developmental stage.
Ordered cells along the pseudotime trajectory and visualized the pseudotime progression.

6. Differential Gene Expression Analysis
Identified differentially expressed genes along the pseudotime trajectory.
Visualized the expression of specific genes of interest in clusters and along pseudotime.

Key Insights
The analysis identified distinct cell clusters in the oral cancer dataset, including specific cell types such as oral mucosa squamous cells.
Batch effects were successfully corrected using the Harmony algorithm.
Pseudotime analysis revealed potential developmental trajectories of cells within the oral cancer samples.
Differential gene expression analysis highlighted genes of interest across different stages of pseudotime.

##How to Use This Repository
Clone the repository to your local machine.
Ensure that the required R packages are installed.
Run the analysis script to reproduce the results.

##Future Work
Further analysis of differentially expressed genes to identify potential biomarkers.
Integration of additional datasets to validate findings.
Functional analysis of identified gene clusters.

##Contact
For any questions or feedback, please open an issue or contact the repository maintainer.

