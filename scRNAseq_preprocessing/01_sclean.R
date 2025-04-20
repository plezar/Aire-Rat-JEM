library(Seurat)
library(scran)
library(tidyverse)
library(scater)

#==== COMMAND-LINE ARGS ====#
args = commandArgs(trailingOnly = T)
data_dir <- args[1]
out_dir <- args[2]
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
}

cluster_seurat <- function(in_data, n_pcs, cluster_res = c(0.6, 0.8, 1.0, 1.2, 1.4)) {
  
  Idents(in_data) <- "RNA"
  
  # run sctransform
  seurat_sct <- SCTransform(in_data,  verbose = FALSE)
  
  # These are now standard steps in the Seurat workflow for visualization and clustering
  seurat_sct <- RunPCA(seurat_sct, verbose = FALSE)
  seurat_sct <- RunUMAP(seurat_sct,
                        reduction = "pca",
                        dims = 1:n_pcs,
                        verbose = FALSE)
  
  # Determine the K-nearest neighbor graph
  seurat_sct <- FindNeighbors(object = seurat_sct, 
                              dims = 1:n_pcs,
                              verbose = FALSE)
  
  # Determine the clusters for various resolutions                                
  seurat_sct <- FindClusters(object = seurat_sct,
                             resolution = cluster_res,
                             verbose = FALSE)
  
  return(seurat_sct)
}

score_cells <- function(meta) {
  # Scores each cell based on "default filtering" described in PipeComp paper
  v.1 <- isOutlier(meta$mitoRatio, batch = meta$quick_clusters, nmads = 2.5, type = "higher")
  #v.2 <- meta$mitoRatio < 0.08
  
  #sc.1 <- as.integer(v.1 | v.2)
  sc.1 <- as.integer(v.1)
  
  v.1 <- isOutlier(log10(meta$nUMI), batch = meta$quick_clusters, nmads = 2.5, type = "higher")
  v.2 <- isOutlier(log10(meta$nUMI), batch = meta$quick_clusters, nmads = 5, type = "lower")
  
  sc.2 <- as.integer(v.1 | v.2)
  
  v.1 <- isOutlier(log10(meta$nGene), batch = meta$quick_clusters, nmads = 2.5, type = "higher")
  v.2 <- isOutlier(log10(meta$nGene), batch = meta$quick_clusters, nmads = 5, type = "lower")
  
  sc.3 <- as.integer(v.1 | v.2)
  
  v.1 <- isOutlier(meta$pct_counts_in_top_20_features, batch = meta$quick_clusters, nmads = 5, type = "higher")
  v.2 <- isOutlier(meta$pct_counts_in_top_20_features, batch = meta$quick_clusters, nmads = 5, type = "lower")
  
  sc.4 <- as.integer(v.1 | v.2)
  
  v.1 <- isOutlier(meta$log10GenesPerUMI, batch = meta$quick_clusters, nmads = 2.5, type = "higher")
  v.2 <- isOutlier(meta$log10GenesPerUMI, batch = meta$quick_clusters, nmads = 5, type = "lower")
  
  sc.5 <- as.integer(v.1 | v.2)
  
  filt.score <- rowSums(cbind(sc.1, sc.2, sc.3, sc.4, sc.5))
  
  return(filt.score)
}


sclean <- function(seurat_obj) {
  cell_counts_vec <- dim(seurat_obj)[2]
  
  # Main assumption is that most of the cells are living: removing cells that are definitely dead
  seurat_obj <- subset(x = seurat_obj, subset = mitoRatio < 0.2)
  

  seurat_obj$quick_clusters <- quickCluster(as.SingleCellExperiment(seurat_obj))
  
  mt_per_cluster <- seurat_obj@meta.data %>%
    group_by(quick_clusters) %>%
    summarise(med_mito = median(mitoRatio))
  
  seurat_obj$filt.score <- score_cells(seurat_obj@meta.data)
  seurat_obj <- seurat_obj[, seurat_obj$filt.score <= 1]
  cell_counts_vec <- c(cell_counts_vec, dim(seurat_obj)[2])
  

  return(seurat_obj)
}

plot_qc_metrics <- function(seurat_obj) {
  p1 <- plot_n_gene(seurat_obj@meta.data)
  p2 <- plot_n_umi(seurat_obj@meta.data)
  p3 <- plot_perc_mt(seurat_obj@meta.data) + geom_vline(xintercept = 0.2)
  p4 <- plot_novely(seurat_obj@meta.data)
  p5 <- plot_gene_vs_umi(seurat_obj@meta.data)
  
  p1 + p2 + p3 + p4 + p5
}

raw_seurat <- readRDS(data_dir)

clean_seurat <- lapply(raw_seurat, sclean)



filtering_stats <- data.frame(Before = t(sapply(raw_seurat, dim))[,2],
                              After = t(sapply(clean_seurat, dim))[,2])

filtering_stats$Perc_removed <- (filtering_stats$Before-filtering_stats$After)/filtering_stats$Before*100

for (i in 1:7) {
  ggsave(plot_qc_metrics(raw_seurat[[i]]), path = paste0(out_dir, "/figures/qc_metrics/unfiltered"), filename=paste0(i, ".png"), width=11, height=7, dpi=700)
  
  ggsave(plot_qc_metrics(clean_seurat[[i]]), path = paste0(out_dir, "/figures/qc_metrics/filtered"), filename=paste0(i, ".png"), width=11, height=7, dpi=700)
}



saveRDS(clean_seurat, paste0(out_dir, "/scleaned_seurat_objects.rds"))
write.table(filtering_stats, file= paste0(out_dir, "/filtering_statistics.tsv"), quote=F, row.names=F, sep="\t")


