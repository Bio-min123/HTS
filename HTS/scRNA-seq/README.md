🧬 Single-cell RNA-seq Analysis Pipeline (Scanpy)

This repository contains a fully reproducible single-cell RNA-seq (scRNA-seq) analysis pipeline implemented in Python using Scanpy.
The workflow follows widely accepted best practices for preprocessing, quality control, clustering, and cell type annotation of 10x Genomics data.

📌 Features

Load 10x Genomics count matrices

Quality control (QC) and filtering

Normalization and log-transformation

Highly variable gene (HVG) selection

Dimensionality reduction (PCA, UMAP)

Graph-based clustering (Leiden)

Marker gene detection

Manual cell type annotation

Publication-quality visualizations

Export of processed data and marker tables

📂 Project Structure
.
├── data/
│   └── filtered_feature_bc_matrix/   # 10x Genomics output
├── results/
│   ├── figures/                       # All generated plots
│   ├── marker_genes.csv               # Differential expression results
│   └── adata_processed.h5ad           # Processed AnnData object
├── scanpy_pipeline.py                 # Main analysis script
└── README.md

🧪 Requirements

Install the required Python packages:

pip install scanpy pandas numpy matplotlib seaborn


Recommended Python version: Python ≥ 3.8

▶️ How to Run

From the project root directory:

python scanpy_pipeline.py


All figures will be automatically saved to:

results/figures/

🧬 Pipeline Overview
Step 1: Load Data

Reads 10x Genomics formatted data using sc.read_10x_mtx

Uses gene symbols as variable names

Ensures gene names are unique

Step 2: Quality Control (QC)

Identifies mitochondrial genes (MT-)

Computes per-cell QC metrics:

Number of genes

Total counts

Percent mitochondrial reads

Filters:

Cells with too few or too many genes

Cells with high mitochondrial content

Genes expressed in fewer than 3 cells

Generates violin plots for QC inspection

Step 3: Normalization

Library-size normalization to 10,000 counts per cell

Log-transformation (log1p)

Stores normalized data in adata.raw

Step 4: Highly Variable Genes (HVGs)

Selects genes based on:

Mean expression

Dispersion (variability)

Retains only HVGs for downstream analysis

Produces HVG diagnostic plots

Step 5: Scaling & PCA

Scales genes to unit variance (with value clipping)

Computes principal components (PCA)

Visualizes:

Variance explained by PCs

PCA colored by mitochondrial percentage

Step 6: Neighbors & UMAP

Constructs a k-nearest neighbor graph

Computes UMAP embedding for visualization

Produces low-dimensional representations of cells

Step 7: Clustering (Leiden)

Performs graph-based clustering using the Leiden algorithm

Visualizes clusters on UMAP

Step 8: Marker Gene Detection

Identifies differentially expressed genes per cluster

Uses the Wilcoxon rank-sum test

Saves results to:

results/marker_genes.csv

Step 9: Cell Type Annotation

Uses known marker genes for major immune cell types

Visualizes marker expression with dot plots

Manually maps clusters to biological cell types

Displays annotated cell types on UMAP

Step 10: Save Results

Saves the fully processed dataset as:

results/adata_processed.h5ad


This file can be reloaded later without rerunning the pipeline.

📊 Outputs
Output	Description
adata_processed.h5ad	Processed single-cell dataset
marker_genes.csv	Differential expression results
QC plots	Cell and gene filtering diagnostics
PCA plots	Dimensionality reduction
UMAP plots	Clusters and cell types
Dot plots	Marker gene expression
🧠 Notes

Cell type annotation is manual and should be adjusted based on marker expression.

Parameters (e.g., QC thresholds, clustering resolution) may need tuning for different datasets.

This pipeline is suitable for exploratory analysis and teaching and can be extended for:

Batch correction

Trajectory inference

Automated annotation

📖 References

Scanpy documentation: https://scanpy.readthedocs.io

10x Genomics scRNA-seq: https://www.10xgenomics.com

Wolf et al., Genome Biology, 2018