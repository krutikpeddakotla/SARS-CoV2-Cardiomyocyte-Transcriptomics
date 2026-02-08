# Differential Gene Expression in SARS-CoV-2 Infected Adult Human Cardiomyocytes

This repository contains a bioinformatics pipeline for analyzing RNA-Seq transcriptomic data from adult human cardiomyocytes infected with SARS-CoV-2. The analysis identifies key genes that are up-regulated or down-regulated in response to viral infection, providing insights into how COVID-19 affects heart tissue at a cellular level.

## Project Overview
The study compares **Mock-infected** (control) vs. **SARS-CoV-2 infected** cardiomyocytes (3 replicates each). The pipeline performs:
* **Preprocessing:** Data cleaning and removal of zero-count genes.
* **Normalization:** Conversion of raw counts to Counts Per Million (CPM).
* **Statistical Analysis:** Log2 transformation, Fold Change calculation, and Welch’s T-Test ($p < 0.05$).
* **Visualization:** Generation of Volcano Plots and Heatmaps for biomarker identification.

## Data Source
The raw data is derived from the **GSE151879** dataset available on the Gene Expression Omnibus (GEO).
* **Tissue:** Adult Human Cardiomyocytes.
* **Comparison:** Mock vs. SARS-CoV-2 Infection (24h post-infection).
## Installation & Requirements
To run the analysis, ensure you have Python 3.x installed along with the following dependencies:
```bash
pip install pandas numpy seaborn matplotlib scipy
