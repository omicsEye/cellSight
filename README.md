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

# To use the sample datafiles available with the package, please download the "datafiles" folder available in the cellSight replo and provide the path for the downloaded datafiles and your desired output file locations repectively:
Example - cellSight("C:/Users/Desktop/cellSight/datafiles/","C:/Users/Desktop/test")

```

## Functions
cellSight is an automated R package, which is an amulgation of several processes for performing numerous downstream analysis. When the main function "cellSight" is invoked, the function recursively calls the rest of the subroutines. 
The various subprogram's under cellSight performs various task:

- **Data Acquisition:** The function fetches all the various single cell datasets located in the input directory and reads the raw files and creates a list of seurat objects.
- Input parameters: The location of the datasets and the desired output directory.
- Output: It creates a "data" folder in the output directory and saves the R object as an RDS file.
```R
cellSight::data_directory("C:/Users/Desktop/cellSight/datafiles/","C:/Users/Desktop/test")
```
- **Quality Control:** Exploring QC Metrics with cellSight

cellSight provides a convenient platform that facilitates the exploration of quality control (QC) metrics and enables the filtering of cells based on user-defined criteria. In the field, several QC metrics are commonly employed, aiding in the identification and removal of problematic cells. Some of these metrics include:

1. **Number of Unique Genes Detected in Each Cell:**
   - Low-quality cells or empty droplets often exhibit very few genes.
   - Cell doublets or multiplets may show an unusually high gene count.

2. **Total Number of Molecules Detected within a Cell:**
   - Correlates strongly with the count of unique genes.
   - Aberrations in this metric may indicate low-quality cells.

3. **Percentage of Reads Mapping to the Mitochondrial Genome:**
   - Low-quality or dying cells may display extensive mitochondrial contamination.

To calculate mitochondrial QC metrics, the `PercentageFeatureSet()` function is employed. This function computes the percentage of counts originating from a specific set of features. In this case, the set of mitochondrial genes is defined as all genes starting with "MT-".


```R
cellSight::qc_plot() 
#Takes the seurat object and the path to the output as parameter to this function

```

<p float="left">
  <img src="plots/qcplot_violin_1.png" width="200" />
  <img src="plots/qcplot_violin_2.png" width="200" />
  <img src="plots/qcplot_violin_3.png" width="200" />
  <img src="plots/qcplot_violin_4.png" width="200" />
</p>

QC Metrics Visualized in the Violin Plot:

1. **Number of Detected Features (Genes) per Cell (nFeature_RNA)**:

    Purpose: This metric represents the number of unique genes detected in each cell.
    
    Interpretation: A violin plot for this metric helps identify cells with unusually low      or high gene detection. Outliers may indicate low-quality cells or potential cell          doublets.
    
2. **Total RNA Molecule Counts per Cell (nCount_RNA)**:

    Purpose: Indicates the total number of RNA molecules detected in each cell.
    
    Interpretation: Examining the distribution of RNA molecule counts helps identify cells     with unusually low or high total counts, which may be indicative of low-quality or         outlier cells.
    
3. **Percentage of Reads Mapping to Mitochondrial Genes (percent.mt)**:

    Purpose: Measures the percentage of RNA reads originating from mitochondrial genes.
    
    Interpretation: High percentages may suggest mitochondrial contamination, often            observed in low-quality or dying cells.
    
*Interpretation of the Violin Plot*:
The violin plot visually represents the distribution of each QC metric across all cells in the dataset. Each vertical violin corresponds to one metric, and the width of the violin at different points indicates the density of cells with specific metric values.

**nFeature_RNA Violin**: Check for outliers or asymmetry, which may indicate cells with aberrant gene detection.

**nCount_RNA Violin**: Assess the spread of total RNA molecule counts, identifying cells with extreme values.

**percent.mt Violin**: Examine the distribution of mitochondrial contamination percentages, identifying cells with potential issues.

<p float="left">
  <img src="plots/qcplot_scatter_1.png" width="200" />
  <img src="plots/qcplot_scatter_2.png" width="200" />
  <img src="plots/qcplot_scatter_3.png" width="200" />
  <img src="plots/qcplot_scatter_4.png" width="200" />
</p>

**Understanding the Plot**:

**X-Axis**: `nCount_RNA`

- **Definition:** Represents the total number of RNA molecules (UMIs) detected in each cell.

- **Interpretation:** Cells with higher values on the x-axis have a greater number of RNA molecules detected, indicating higher RNA content.

**Y-Axis**: `nFeature_RNA`

- **Definition:** Represents the number of unique genes detected in each cell.

- **Interpretation:** Cells with higher values on the y-axis have a greater diversity of genes expressed, indicating a higher complexity of transcriptional activity.

**Scatter Plot**:

- **Purpose:** The scatter plot visualizes the relationship between the two QC metrics (`nCount_RNA` and `nFeature_RNA`) for each cell in the dataset.

- **Interpretation:** 
  
  - Cells clustering towards the upper-right corner of the plot have both high RNA molecule counts (`nCount_RNA`) and a large number of detected unique genes (`nFeature_RNA`). These cells are typically considered of higher quality, exhibiting robust transcriptional activity.
  
  - Cells clustering towards the lower-left corner of the plot have lower RNA molecule counts and fewer detected unique genes, possibly indicating low-quality or empty cells.

**Conclusion**:
The scatter plot generated by the `FeatureScatter` function provides valuable insights into the distribution of two important QC metrics, `nCount_RNA` and `nFeature_RNA`, across the single-cell RNA-seq dataset. It enables researchers to identify and assess the quality of individual cells based on their RNA content and transcriptional complexity.

<p float="left">
  <img src="plots/qcplot_grid_joined_1.png" width="200" />
  <img src="plots/qcplot_grid_joined_2.png" width="200" />
  <img src="plots/qcplot_grid_joined_3.png" width="200" />
  <img src="plots/qcplot_grid_joined_4.png" width="200" />
</p>
The above plot shows the joined distribution of the counts for each feature in the data and the violin plot for the distribution of the number of feature, number of counts and the mitochondria percentage present in the data.

- **Filtering:** Data filtering to ensure accurate downstream analysis by removing doublets.

```R
cellSight::filtering() The parameters for this function are the seurat object and the output destination. It filters and the data within a given range. 

```

The filtering function performs data subsetting in cellSight, filtering cells based on specified quality control (QC) criteria. Let's break down the conditions used for subsetting:

- `nFeature_RNA > min_rna`: Selects cells with a number of detected unique genes (`nFeature_RNA`) greater than a specified minimum threshold (`min_rna`).
- `nFeature_RNA < nfeature_rna`: Selects cells with a number of detected unique genes (`nFeature_RNA`) less than a specified maximum threshold (`nfeature_rna`).
- `nCount_RNA < ncount_rna`: Selects cells with a total number of RNA molecules (UMIs) (`nCount_RNA`) less than a specified maximum threshold (`ncount_rna`).
- `percent.mt < mit.percent`: Selects cells with a percentage of reads mapping to the mitochondrial genome (`percent.mt`) less than a specified threshold (`mit.percent`).

**Interpretation**:
This code is used to create a subset of cells from the original dataset based on specific QC metrics. The conditions specified in the subset are designed to filter out cells that may exhibit characteristics associated with low quality or potential contamination.

- **Filtered Cells:** The resulting subset includes cells that meet all the specified conditions, indicating they have a desirable range of unique gene counts, total RNA molecules, and mitochondrial contamination levels.

- **Downstream Analysis:** This subset of cells can then be used for downstream analyses, such as clustering, dimensionality reduction, or differential expression analysis, with increased confidence in the quality of the selected cells.

**Conclusion**:
Subsetting based on QC criteria is a crucial step in the analysis of single-cell RNA-seq data. This process allows researchers to focus on a subset of high-quality cells for further investigation and analysis, improving the accuracy and reliability of downstream results.

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
