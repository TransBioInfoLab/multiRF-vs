## An Integrative Multi-Omics Random Forest Framework for Robust Biomarker Discovery

Wei Zhang, Hanchen Huang, Lily Wang, Brian D. Lehmann, Steven X. Chen

This repo contains code associated with the manuscript, our method enhanced the detection of key variables in multi-omics integration.

### Abstract

High-throughput technologies now produce a wide array of omics data, from genomic and transcriptomic profiles to epigenomic and proteomic measurements. Integrating these diverse data types can yield deeper insights into the biological mechanisms driving complex traits and diseases. Yet, extracting key shared biomarkers from multiple data layers remains a major challenge. We present a multivariate random forest (MRF)–based framework enhanced by a novel inverse minimal depth (IMD) metric for integrative variable selection. By assigning response variables to tree nodes and employing IMD to rank predictors, our approach efficiently identifies essential features across different omics types, even when confronted with high-dimensionality and noise. Through extensive simulations and analyses of multi-omics datasets from The Cancer Genome Atlas, we demonstrate that our method outperforms established integrative techniques in uncovering biologically meaningful biomarkers and pathways. Our findings show that selected biomarkers not only correlate with known regulatory and signaling networks but can also stratify patient subgroups with distinct clinical outcomes. The method’s scalable, interpretable, and user-friendly implementation ensures broad applicability to a range of research questions. This MRF-based framework advances robust biomarker discovery and integrative multi-omics analyses, accelerating the translation of complex molecular data into tangible biological and clinical insights. 

### Code Directory

**Method**

- `code/function/method`

  - `mrf3_init.R`: The initial multi-omics MRF model fitting function. Runs before variable selection 

  - `mrf3_vs.R`: Main function for MRF-based multi-omics variable selection
  - `mrf3-util.R`: Utility functions supporting initial model execution and variable selection.

**Simulation**

- `code/simulation`: Runs simulation of the MRF-based methods and other benchmark methods
- `code/function/utility`: Simulation utilities (data generation, simulation wrappers)

**Real data analysis**

- `data_prepare.Rmd`: TCGA data preprocessing
- `$DATA$.R`: Performs model fitting for real data analysis on each TCGA data
- `real_data_analysis.Rmd`: Contains figures generated for the manuscript

### Instruction

To conduct MRF multi-omics variable selection using the method described in the manuscript, 


### For Reproducible Research

Install all the R packages from the `load_package.R` file and source all the scripts in the `code/function` folder.

```R
─ Session info ───────────────────────────────────────────────
 setting  value
 version  R version 4.4.2 (2024-10-31)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2024-11-20
 rstudio  2024.04.2+764.pro1 Chocolate Cosmos (server)
 pandoc   3.1.11 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
```
