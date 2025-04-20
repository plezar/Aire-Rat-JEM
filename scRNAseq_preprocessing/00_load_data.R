# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

#==== COMMAND-LINE ARGS ====#
args = commandArgs(trailingOnly = T)
data_dir <- args[1]
doublet_dir <- args[2]
out_dir <- args[3]
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
}

add_meta_seurat <- function(obj) {
  # Add number of genes per UMI for each cell to metadata
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  obj$pct_counts_in_top_20_features <- sapply(Cells(obj), function(cell){sum(sort(obj@assays$RNA@counts[,cell], decreasing = TRUE)[1:20])/sum(obj@assays$RNA@counts[,cell])})
  
  # Compute percent mito ratio
  obj$mitoRatio <- PercentageFeatureSet(object = obj, pattern = "^Mt-")
  obj$mitoRatio <- obj@meta.data$mitoRatio / 100
  
  obj$rpRatio <- PercentageFeatureSet(object = obj, pattern = "^Rp")
  obj$rpRatio <- obj$rpRatio / 100
  
  load("/Users/plezar/Documents/molpat/bioinf/single_cell/cycle_genes.rda")
  # Score cells for cell cycle
  obj <- CellCycleScoring(obj, 
                          g2m.features = g2m_genes,
                          s.features = s_genes)
  
  # Create metadata dataframe
  metadata <- obj@meta.data
  
  # Add cell IDs to metadata
  metadata$cells <- rownames(metadata)
  
  # Create sample column
  metadata$sample <- metadata$orig.ident
  
  # Rename columns
  metadata <- metadata %>%
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  
  # Add metadata back to Seurat object
  obj@meta.data <- metadata
  return(obj)
}

sample_ids <- list.files(data_dir)

# Create each individual Seurat object for every sample
for (file in sample_ids){
  doublets <- readLines(paste0(doublet_dir, file, ".txt"))
  seurat_data <- Read10X(data.dir = paste0(data_dir, file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   #min.features = 100, 
                                   project = file)
  seurat_obj <- add_meta_seurat(seurat_obj)
  assign(file, seurat_obj[,grepl(0, doublets)])
  #assign(file, seurat_obj)
}

# Create a merged Seurat object
merged_seurat_raw <- merge(x = T_HE1,
                       y = c(T_HE3, T_HE7, T_KO2, T_KO4, T_KO6, T_KO8),
                       add.cell.id = sample_ids,
                       project = "rat_tecs")

raw_seurat <- list(T_HE1, T_HE3, T_HE7, T_KO2, T_KO4, T_KO6, T_KO8)
names(raw_seurat) <- sample_ids

saveRDS(raw_seurat, paste0(out_dir, "seurat_raw.rds"))
