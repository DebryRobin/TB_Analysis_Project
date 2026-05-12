# ==============================================================================
# PROJECT  : TB Airway Immune Signatures Analysis (GSE328391)
# MODULES  : Metadata Fetching | DESeq2 | Advanced Plotting 
# AUTHOR   : [Your Name]
# ==============================================================================

# ── 0. SETUP & DEPENDENCIES ──────────────────────────────────────────────────
# This block automatically installs and loads everything needed.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "DESeq2", "GEOquery", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", 
  "enrichplot", "ggplot2", "ggrepel", "pheatmap", "dplyr", "tidyr", "tibble"
)

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# Create directory structure properly
out_dir <- "results_output"
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

cat("✓ Environment and directories ready\n")

# ── 1. ROBUST METADATA LOADING (GSE328391) ───────────────────────────────────
# Fetches patient information directly from NCBI GEO
load_metadata <- function(gse_id) {
  message("--> Fetching metadata from GEO...")
  gse <- getGEO(gse_id, destdir = "data/", getGPL = FALSE)
  metadata <- pData(gse[[1]])
  
  # FIX: Automatic column detection to avoid "Column doesn't exist" errors
  # We look for the column describing the patient group (Progressor vs Control)
  status_col <- grep("disease|status|group|type|contact", colnames(metadata), ignore.case = TRUE, value = TRUE)[1]
  
  if (is.na(status_col)) {
    stop("Could not find a 'Status' column in metadata. Check colnames(metadata).")
  }
  
  message("--> Detected status column: ", status_col)
  
  clean_meta <- metadata %>%
    select(geo_accession, title, all_of(status_col)) %>%
    rename(SampleID = geo_accession, Status = !!sym(status_col)) %>%
    mutate(Status = as.factor(Status))
  
  rownames(clean_meta) <- clean_meta$SampleID
  return(clean_meta)
}

# ── 2. CORE ANALYSIS PIPELINE ────────────────────────────────────────────────
# Uses your existing counts.csv file
run_analysis <- function(counts_path, metadata) {
  if (!file.exists(counts_path)) {
    stop("Missing counts file at: ", counts_path, ". Please place your CSV in the /data folder.")
  }
  
  message("--> Loading counts and running DESeq2...")
  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  
  # Sync Samples: Match CSV columns with GEO metadata rows
  common <- intersect(colnames(counts), rownames(metadata))
  if(length(common) == 0) stop("No matching sample IDs found between CSV and Metadata!")
  
  counts <- counts[, common]
  metadata <- metadata[common, ]
  
  # Create DESeq Object
  dds <- DESeqDataSetFromMatrix(countData = round(counts), 
                                colData = metadata, 
                                design = ~ Status)
  
  # Noise filtering: Keep genes with sufficient reads
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run Statistics
  dds <- DESeq(dds)
  return(dds)
}

# ── 3. MAIN EXECUTION ────────────────────────────────────────────────────────

# Step A: Load patient data
meta <- load_metadata("GSE328391")

# Step B: Run the math (Assumes your file is in data/counts.csv)
counts_file <- "data/counts.csv"

if (file.exists(counts_file)) {
  dds <- run_analysis(counts_file, meta)
  res <- results(dds, alpha = 0.05)
  vst_data <- vst(dds, blind = FALSE)
  
  # Step C: Annotation (Ensembl IDs -> Gene Symbols)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    mutate(
      direction = case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "UP",
        padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
        TRUE                               ~ "NS"
      )
    )
  
  gene_anno <- AnnotationDbi::select(org.Hs.eg.db, keys = res_df$ensembl_id, 
                                     columns = c("SYMBOL"), keytype = "ENSEMBL") %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  res_ann <- left_join(res_df, gene_anno, by = c("ensembl_id" = "ENSEMBL"))
  
  # ── 4. VISUALIZATION ───────────────────────────────────────────────────────
  
  # FIGURE 1: Volcano Plot
  fig1 <- ggplot(res_ann, aes(x=log2FoldChange, y=-log10(padj), color=direction)) +
    geom_point(alpha=0.5, size=1.5) +
    scale_color_manual(values=c("UP"="#E41A1C", "DOWN"="#377EB8", "NS"="grey70")) +
    theme_bw() +
    labs(title="TB Analysis: Volcano Plot", subtitle="GSE328391 | Progressors vs Controls")
  
  ggsave(file.path(out_dir, "figures/Figure1_Volcano.png"), fig1, width=8, height=6)
  
  # FIGURE 2: PCA (Clustering)
  fig2 <- plotPCA(vst_data, intgroup="Status") + 
    theme_bw() + 
    labs(title="Sample Clustering (PCA)")
  
  ggsave(file.path(out_dir, "figures/Figure2_PCA.png"), fig2, width=8, height=6)
  
  # Save Final Table
  write.csv(res_ann, file.path(out_dir, "tables/Final_Differential_Expression.csv"))
  
  cat("\n Pipeline complete. Check the '", out_dir, "' folder for results.\n")
  
} else {
  cat("\n Waiting for data. Place your 'counts.csv' in the 'data/' folder to start.\n")
}