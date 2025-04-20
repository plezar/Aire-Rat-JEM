library(Seurat)

#==== COMMAND-LINE ARGS ====#
args = commandArgs(trailingOnly = T)
data_dir <- args[1]
out_dir <- args[2]
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
}

integrate_seurat <- function(obj_list) {
  
  obj_list <- SplitObject(obj_list, split.by = "sample")
  options(future.globals.maxSize = 4000 * 1024^2)
  
  for (i in 1:length(obj_list)) {
    obj_list[[i]] <- SCTransform(obj_list[[i]])
  }
  
  
  integ_features <- SelectIntegrationFeatures(object.list = obj_list, 
                                              nfeatures = 3000) 
  
  # Prepare the SCT list object for integration
  obj_list <- PrepSCTIntegration(object.list = obj_list, 
                                 anchor.features = integ_features)
  
  # Find best buddies - can take a while to run
  integ_anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                          normalization.method = "SCT", 
                                          anchor.features = integ_features)
  
  # Integrate across conditions
  obj <- IntegrateData(anchorset = integ_anchors, 
                       normalization.method = "SCT")
  return(obj)
}

clean_seurat <- readRDS(data_dir)
clean_seurat <- integrate_seurat(clean_seurat)
#clean_seurat$genotype <- str_split(clean_seurat$seq_folder, "_", simplify = TRUE)[,1]
saveRDS(clean_seurat, paste0(out_dir, "/integrated_seurat_object.rds"))
