# CytoExploreR-Interactive-visualization
CytoExploreR: A Comprehensive Open-Source R Package for Interactive Cytometry Data Analysis

```markdown
# Interactive Analysis of Cytometry Data in R

**Author:** Pritam Kumar Panda  
**Institution:** German Cancer Research Center (DKFZ)  
**Country:** Germany  

This project is dedicated to providing a comprehensive guide for the interactive analysis of cytometry data using R. It leverages a variety of R packages to preprocess, analyze, and visualize cytometry data effectively.

## Prerequisites

Before you can run the scripts, make sure you have R installed on your system. You will also need to install some dependencies.

### Set Your Working Directory

```r
setwd("/Users/pritam/facs_data_analysis_R/FACS/CytoExploreR")
```

### Initialize renv (Optional)

```r
# renv::init()
```

### Increase Memory and Use Cache for Larger Files

```r
knitr::opts_chunk$set(cache = TRUE, warning = FALSE, message = FALSE, cache.lazy = FALSE)
```

### Install Required Libraries

```r
# Install these packages if not already installed
install.packages("tidyverse")
install.packages("knitr")
install.packages("ggplot2")
install.packages("remotes")
install.packages("BiocManager")
install.packages("devtools")
install.packages("kableExtra")

# Bioconductor packages
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("openCyto")
BiocManager::install("flowAI")
BiocManager::install("ggcyto")
BiocManager::install("CytoML")

# CytoExploreR
devtools::install_github("DillonHammill/CytoExploreR")
devtools::install_github("DillonHammill/CytoExploreRData")
```

### Load Required Packages

```r
library(flowCore)
library(CytoML)
library(flowAI)
library(flowWorkspace)
library(ggcyto)
library(tidyverse)
library(openCyto)
library(knitr)
library(kableExtra)
library(dplyr)
library(CytoExploreR)
library(CytoExploreRData)
```

## Data Preparation

### Save FCS Files

- **Compensation FCS Files**

  ```r
  cyto_save(Compensation, save_as = "Compensation-Samples")
  ```

- **Activation FCS Files**

  ```r
  cyto_save(Activation, save_as = "Activation-Samples")
  ```

## Analysis Workflow

1. **Compensation of Fluorescent Spillover**

   - Prepare and compute the spillover matrix.
   - Interactively edit spillover matrices.

2. **Sample Analysis**

   - Load, annotate, and compensate samples.
   - Perform gating strategies for cell populations.

3. **Visualizations**

   - Visualize compensation, gating trees, and gating schemes.


## Acknowledgments

- Thanks to the creators of the R packages used in this project. https://github.com/DillonHammill/CytoExploreR/
```

This README.md template provides a structured way to present your project on GitHub, including how to set up the environment, install dependencies, prepare data, and proceed through the analysis workflow. Adjust the paths and commands according to your specific project needs.
