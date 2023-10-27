# cellSight

cellSight is an R package that offers a fully automated pipeline for the analysis of single-cell RNA sequencing (scRNA-seq) data. This package is designed to streamline the analysis process, making it accessible and efficient for both beginners and experienced bioinformaticians. The best part? It requires only the data files, taking care of the rest.

* **Citation:** if you use the cellSight software, please cite: Ranojoy Chatterjee, Chiraag Gohel, and Ali Rahnavard (2023+). **cellSight: Deciphering underlying dynamics of cells and gene expression.** R package version 0.0.1. https://www.rahnavard.org/cellSight.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Example](#example)
- [Contributing](#contributing)

## Installation
Install Bioconductor packages before installing cellSight on the R console:
```R
# Install Bioconductor Manager and install DelayedMatrixStats
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DelayedMatrixStats")

#Install glmGamPoi
BiocManager::install("glmGamPoi")

```

To install cellSight, you'll need have the 'devtools' package installed. Run the following commands in your R console:

```R
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install cellSight from GitHub
devtools::install_github("omicsEye/cellSight")
```

## Usage

Using cellSight is straightforward. Load the package, provide your data files, and let it do the magic. Here's a basic example:

```R
# Load the cellSight package
library(cellSight)

# Provide your data files (e.g., expression matrix )

# Run cellSight
results <- cellSight(expression_matrix)

# Explore the results and visualize the output
```

## Features

cellSight comes with a wide range of features, including:

- **Quality Control:** Automated quality control (QC) checks for cells and genes.
- **Normalization:** Data normalization to ensure accurate downstream analysis.
- **Clustering:** Cell clustering based on gene expression profiles.
- **Differential Expression Analysis:** Identifying differentially expressed genes between clusters.
- **Visualization:** Visualize the results with various plots, including t-SNE, UMAP, and more.
- **Annotation:** Annotate cell clusters with known cell types or states.
- **Cell communication** Explore cell communication by infering the possible Ligand receptor interactions.



## Contributing

We welcome contributions from the community. If you'd like to contribute to cellSight, please follow our [Contributing Guidelines](CONTRIBUTING.md).

---

Thank you for choosing cellSight for your scRNA-seq analysis. We hope this tool simplifies your research and helps you gain valuable insights from your single-cell data. If you have any questions, encounter issues, or want to contribute, please don't hesitate to get in touch with us.

For updates and news, follow us on Twitter [@RahnavardLab](https://twitter.com/RahnavardLab). Happy analyzing!
