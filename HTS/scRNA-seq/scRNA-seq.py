import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# -----------------------------
# Settings
# -----------------------------
sc.settings.verbosity = 3
sc.settings.figdir = "results/figures"
sc.settings.set_figure_params(dpi=100, facecolor="white")

DATA_DIR = Path("data/filtered_feature_bc_matrix")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# -----------------------------
# Step 1: Load data
# -----------------------------
adata = sc.read_10x_mtx(
    DATA_DIR,
    var_names="gene_symbols",
    cache=True
)

adata.var_names_make_unique()

print(adata)

# -----------------------------
# Step 2: Quality Control (QC)
# -----------------------------
# Mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# QC plots
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True
)

# Filtering cells
adata = adata[
    (adata.obs.n_genes_by_counts > 200) &
    (adata.obs.n_genes_by_counts < 6000) &
    (adata.obs.pct_counts_mt < 10),
    :
]

# Filtering genes
sc.pp.filter_genes(adata, min_cells=3)

print("After QC:", adata)

# -----------------------------
# Step 3: Normalization
# -----------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save raw data
adata.raw = adata

# -----------------------------
# Step 4: Highly Variable Genes
# -----------------------------
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var.highly_variable]

# -----------------------------
# Step 5: Scaling & PCA
# -----------------------------
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")

sc.pl.pca_variance_ratio(adata, log=True)
sc.pl.pca(adata, color="pct_counts_mt")

# -----------------------------
# Step 6: Neighbors & UMAP
# -----------------------------
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)

sc.pl.umap(adata)

# -----------------------------
# Step 7: Clustering (Leiden)
# -----------------------------
sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata, color=["leiden"])

# -----------------------------
# Step 8: Marker Gene Detection
# -----------------------------
sc.tl.rank_genes_groups(
    adata,
    "leiden",
    method="wilcoxon"
)

sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# Save marker table
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv("results/marker_genes.csv", index=False)

# -----------------------------
# Step 9: Cell Type Annotation (manual)
# -----------------------------
marker_dict = {
    "T_cells": ["CD3D", "CD3E", "TRBC1"],
    "B_cells": ["MS4A1", "CD79A"],
    "NK_cells": ["NKG7", "GNLY"],
    "Monocytes": ["LYZ", "S100A8", "S100A9"]
}

sc.pl.dotplot(adata, marker_dict, groupby="leiden")

# Example mapping (adjust after inspection)
cluster_to_celltype = {
    "0": "T_cells",
    "1": "B_cells",
    "2": "Monocytes",
    "3": "NK_cells"
}

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_celltype)

sc.pl.umap(adata, color="cell_type")

# -----------------------------
# Step 10: Save results
# -----------------------------
adata.write("results/adata_processed.h5ad")
