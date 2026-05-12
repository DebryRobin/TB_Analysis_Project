# TB Airway Immune Signature Analysis (GSE328391)

This project performs a complete RNA-seq analysis comparing human tuberculosis household contacts who progressed to active disease versus those who remained healthy.

## Prerequisites
- **R** (version 4.0 or higher)
- **Data**: A `counts.csv` file containing the gene expression matrix (place in `/data`).

## Pipeline Steps
1. **Metadata Integration**: Automatically downloads clinical characteristics from NCBI GEO using `GEOquery`.
2. **Normalization**: Uses `DESeq2` to normalize raw counts and estimate data dispersion.
3. **Statistical Testing**: Identifies differentially expressed genes (DEGs) with an FDR < 0.05.
4. **Annotation**: Maps Ensembl gene IDs to official Gene Symbols.
5. **Visualization**: Generates Volcano plots and PCA clustering for analysis.

## Usage
Open RStudio, set your working directory to this folder, and run:
`source("main_analysis.R")`