# ChIP-seq Analysis: RUNX1 Binding in MCF-7 Breast Cancer Cells
[![Nextflow](https://img.shields.io/badge/nextflow-%23EAD69C.svg?style=for-the-badge&logo=nextflow&logoColor=000000)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.bbagrm.2016.100-00A3E0)](https://doi.org/10.1016/j.bbagrm.2016.100)

A comprehensive re-analysis of ChIP-seq data from Barutcu et al. (2016) investigating RUNX1's role as a transcriptional regulator in breast cancer cells, implemented as a reproducible Nextflow pipeline.

## Project Overview
This project implements a reproducible Nextflow pipeline for ChIP-seq analysis, focusing on RUNX1 transcription factor binding in MCF-7 breast cancer cells. The analysis includes quality control, peak calling, motif discovery, and integration with gene expression data to identify direct RUNX1 targets.

**Key Findings:**
- Identified 6,462 high-confidence RUNX1 binding sites
- Discovered RUNX1's role as a master regulator coordinating multiple cancer hallmarks
- Validated direct transcriptional regulation of key breast cancer genes
- Revealed extensive involvement in cell cycle, DNA damage response, and metabolic pathways
### Prerequisites
- Nextflow >= 22.10.0
- Singularity/Apptainer or Docker
- Java 8 or higher

### Installation
```bash
git clone https://github.com/yourusername/chipseq-runx1-analysis.git
cd chipseq-runx1-analysis
```

### Running the Pipeline
```bash
nextflow run main.nf -profile singularity --samplesheet samplesheet.csv
```
### Project Strucutre

``` text
chipseq-runx1-analysis/
├── main.nf                          # Main Nextflow pipeline
├── modules/                         # Custom Nextflow modules
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bowtie2_align/
│   ├── homer_findpeaks/
│   └── ...
├── results/                         # Analysis outputs
│   ├── MultiQC_Report.pdf          # Comprehensive QC report
│   ├── profile_plot.png            # Signal coverage across genes
│   ├── known_homer_motif.png       # Motif enrichment results
│   ├── figure_2F_accurate.png      # Integration with expression data
│   └── ...
├── config/                         # Configuration files
└── docs/                          # Documentation
```

## Citation
Barutcu, A. R., Hong, D., Lajoie, B. R., McCord, R. P., van Wijnen, A. J., Lian, J. B., Stein, J. L., Dekker, J., Imbalzano, A. N., & Stein, G. S. (2016). RUNX1 contributes to higher-order chromatin organization and gene regulation in breast cancer cells. Biochimica et biophysica acta, 1859(11), 1389–1397. https://doi.org/10.1016/j.bbagrm.2016.08.003
