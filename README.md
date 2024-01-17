# cellSight

cellSight is an R package for automated single-cell RNA sequencing (scRNA-seq) data analysis. It streamlines the process for both beginners and experienced bioinformaticians, requiring only data files.

* **Citation:** if you use the cellSight software, please cite: Ranojoy Chatterjee, Chiraag Gohel, and Ali Rahnavard (2023+). **cellSight: Deciphering underlying dynamics of cells and gene expression.** R package version 0.0.1. https://www.rahnavard.org/cellSight.
[![DOI](https://zenodo.org/badge/429576005.svg)](https://zenodo.org/doi/10.5281/zenodo.10041146)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Funtions](#functions)
- [Example](#example)
- [Contributing](#contributing)

## Installation
Install Bioconductor packages before installing cellSight on the R console:
```R
# Install Bioconductor Manager and required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DelayedMatrixStats", "glmGamPoi", "metap", "multtest"))

```

To install cellSight using devtools:

```R
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install cellSight from GitHub
devtools::install_github("omicsEye/cellSight")
```

## Usage

For cellSight, ensure your transcriptomics data is in the format:

Single dataset: <br>"path/to/data/outs/filtered_feature_bc_matrix/sample"

Multiple datasets:<br>"path/to/data/outs/filtered_feature_bc_matrix/sample1",<br>                         "path/to/data/outs/filtered_feature_bc_matrix/sample2",
<br>
                  "path/to/data/outs/filtered_feature_bc_matrix/sample3"

```R
# Load the cellSight package from the github repo
library(cellSight)

# Provide your data files (e.g., expression matrix )

# Run cellSight
cellSight("path/to/data", "path/to/output")


```

## Functions

cellSight comes with a wide range of features, including:

- **Quality Control:** Automated quality control (QC) checks for cells and genes.
```R
cellSight::qc_plot() 
#Takes the seurat object and the path to the output as parameter to this function

```
**The distribution of the number of feature, number of counts and the mitochondria percentage present in the data.**

<img src="plots/qcplot_violin_1.png" width="900" height="900">

**The distribution of the counts for each feature in the data**

<img src="plots/qcplot_scatter_1.png" width="900" height="900">

**The joined distribution of the counts for each feature in the data and the violin plot for the distribution of the number of feature, number of counts and the mitochondria percentage present in the data**

<img src="plots/qcplot_grid_joined_1.png" width="1000" height="900">

- **Filtering:** Data filtering to ensure accurate downstream analysis by removing doublets.

```R
cellSight::filtering() The parameters for this function are the seurat object and the output destination. It filters and the data within a given range. 

```
- **Clustering:** Cell clustering based on gene expression profiles.
```R
cellSight::pca_clustering() The fnction expects 2 parameters. 1) The seurat object 2) The output path

```
<p float="left">
  <img src="plots/integrated_snn_res(0.2).png" width="250" />
  <img src="plots/integrated_snn_res(0.6).png" width="250" />
  <img src="plots/integrated_snn_res(0.8).png" width="250" />
</p>


- **Differential Expression Analysis:** Identifying differentially expressed genes between clusters.
- **Visualization:** Visualize the results with various plots, including t-SNE, UMAP, and more.
- **Annotation:** Annotate cell clusters with known cell types or states.
- **Cell communication** Explore cell communication by inferring the possible Ligand-Receptor interactions.



## Contributing

We welcome contributions from the community. If you'd like to contribute to cellSight, please follow our [Contributing Guidelines](CONTRIBUTING.md).

---

Thank you for choosing cellSight for your scRNA-seq analysis. We hope this tool simplifies your research and helps you gain valuable insights from your single-cell data. If you have any questions, encounter issues, or want to contribute, please don't hesitate to get in touch with us.

For updates and news, follow us on Twitter [@RahnavardLab](https://twitter.com/RahnavardLab). Happy analyzing!
