# An Integrative Multi-Omics Random Forest Framework for Robust Biomarker Discovery

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

- `data_prepare.Rmd`: Real data download and preprocessing of TCGA data
- `BRCA_COAD.Rmd`: Comprehensive Analysis of Individual Cancer Data: Breast Cancer and Colorectal Cancer
- `PAN.R`: TCGA PAN Cancer Clustering Analysis
- `ADNI.Rmd`: Integrative Analysis Enhances Prediction of Dementia Progression in the ADNI Cohort

### Instruction

For instructions on performing MRF multi-omics variable selection as described in the manuscript, a detailed vignette can be found here: [[Link](http://rpubs.com/noblegasss/multiRF-vs-vignette)]

### Data Accessibility

- TCGA data
  - RNAseq datasets: from R package TCGAbiolinks
  - Others: from UCSC Xena: [https://xena.ucsc.edu/]([https://xena.ucsc.edu/)

- ADNI data: [adni.loni.usc.edu](adni.loni.usc.edu) 

| Dataset   | Data Type                         | Link                                                         |
| --------- | --------------------------------- | ------------------------------------------------------------ |
| TCGA-BRCA | Gene Expression                   | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | miRNA                             | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FmiRNA_HiSeq_gene&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | DNA Methylation                   | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHumanMethylation450&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | Survival Data                     | [View Data](https://xenabrowser.net/datapages/?dataset=survival%2FBRCA_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | Phenotype                         | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
| TCGA-COAD | Gene Expression                   | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | miRNA                             | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FmiRNA_HiSeq_gene&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | DNA Methylation                   | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FHumanMethylation450&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | Survival Data                     | [View Data](https://xenabrowser.net/datapages/?dataset=survival%2FCOAD_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | Phenotype                         | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
| TCGA-PAN  | ATAC-Seq                          | [View Data](https://xenabrowser.net/datapages/?dataset=TCGA_ATAC_peak_Log2Counts_dedup_sample&host=https%3A%2F%2Fatacseq.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
|           | RNA-Seq                           | [View Data Download Code](https://github.com/TransBioInfoLab/multiRF-vs/blob/main/code/real_data/data_prepare.Rmd) |
| ADNI      | Gene Expression & DNA Methylation | [ADNI study website](adni.loni.usc.edu)                      |


### For Reproducible Research

Install all the R packages from the `load_package.R` file and source all the scripts in the `code/function` folder.

License: GPL-3.0

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
