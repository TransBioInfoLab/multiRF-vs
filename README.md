# multiRF-variable-selection
## Enhanced Multi-omics Integration via Multivariate Random Forests for Robust Biomarker Discovery

Wei Zhang, Hanchen Huang, Lily Wang, Brian D. Lehmann, Steven X. Chen

This repo contains code associated with the manuscript, our method enhanced the detection of key variables in multi-omics integration.

### Abstract

The rapid advancement of high-throughput omics technologies has led to an explosion of multi-dimensional biological data. Integrating these diverse data types holds the potential to uncover novel insights into complex biological systems. However, identifying shared and relevant variables across different omics datasets remains a significant challenge. Here, we propose a novel approach using multivariate random forests to perform variable selection across multiple omics datasets. Our method, which leverages both maximal splitting response variable (MSRV) and inverse minimal depth (IMD), effectively identifies key variables by reducing dimensionality while maintaining predictive power. We demonstrate its superiority over traditional methods like sparse partial least squares (sPLS) and canonical correlation analysis (CCA) through comprehensive simulations and real-world multi-omics data. Our results show that the proposed method can robustly detect biologically relevant features, even in high-dimensional, noisy datasets. This tool offers a scalable, interpretable, and accurate framework for multi-omics data integration, with wide applications in biomarker discovery and disease prediction.

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

- `data_prepare.Rmd`: Data preprocessing
- `$DATA$.R`: Performs model fitting for real data analysis on each TCGA data
- `real_data_analysis.Rmd`: Contains figures generated for the manuscript

### For Reproducible Research

Install all the R packages from the `load_package.R` file and source all the scripts in the `code/function` folder.

```
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
