generate_pseudobulk_samples <- function(seurat_obj, filter=F) {
  
  ### generate_pseudobulk_samples
  ### 
  ### This function generates pseudobulk expression matrices from a Seurat object.
  ### It converts the Seurat object to a SingleCellExperiment format, aggregates
  ### raw counts by sample and cell type, optionally filters lowly expressed genes,
  ### normalizes the data using TMM, computes log-transformed CPM values, and then
  ### averages expression values by cell type to produce a final pseudobulk matrix.
  ###
  ### Args:
  ###   seurat_obj: A Seurat object with "sample", "cell.type", and "genotype" metadata.
  ###   filter: Logical; if TRUE, filter out lowly expressed genes using edgeR::filterByExpr.
  ###
  ### Returns:
  ###   A gene-by-cell type expression matrix (log2 CPM), averaged across samples.
  
  merged_sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
  merged_sce@assays@data$logcounts <- NULL
  summed <- aggregateAcrossCells(merged_sce, 
                                 id=colData(merged_sce)[,c("sample", "cell.type")],
                                 use.assay.type = 1)
  
  y <- DGEList(counts(summed), samples=colData(summed))
  
  if (filter) {
    keep <- filterByExpr(y, group=summed$genotype)
    y <- y[keep,]
  }
  
  y <- calcNormFactors(y)
  
  mat <- cpm(y, log = T, prior.count = 1)
  colnames(mat) <- 1:ncol(mat)
  
  exprs_mat <- mat %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(!Gene, names_to = "Sample", values_to = "Expression") %>%
    left_join(
      data.frame(Sample = as.character(1:nrow(colData(summed))), Cell.type = colData(summed)$cell.type, Sample_name = colData(summed)$sample),
      by="Sample"
    ) %>%
    dplyr::select(-Sample, Sample = Sample_name) %>%
    group_by(Gene, Cell.type) %>%
    summarise(Expression = mean(Expression)) %>%
    pivot_wider(names_from = "Cell.type", values_from = "Expression") %>%
    column_to_rownames("Gene")
  return(exprs_mat)
}

plot_smoothers_fun <- function(sce_obj, gene_name) {
  
  ### plot_smoothers_fun
  ###
  ### This function generates a smooth expression trend plot for a specified gene 
  ### using a tradeSeq-fitted SingleCellExperiment object.
  ### It visualizes gene expression along inferred trajectories with customized colors.
  ###
  ### Args:
  ###   sce_obj: A SingleCellExperiment object with fitted GAM curves (from tradeSeq).
  ###   gene_name: Character string; the name of the gene to plot.
  ###
  ### Returns:
  ###   A ggplot object showing smoothed expression curves for the specified gene.
  
  curvesCols <- rev(met.brewer("Nizami", 2))
  
  plotSmoothers(sce_obj, assays(sce_obj)$counts, gene_name, curvesCols = curvesCols, border = FALSE) +
    scale_color_manual(values = curvesCols) + ggtitle(gene_name) + theme_classic(base_size = 14)
}

dot_plot <- function (object, features, assay = NULL, cols = c("lightgrey", 
                                                   "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
          idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
          scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  ### DotPlot function: see Seurat documentation
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), 
                                        idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- FetchData(object = object, vars = split.by)[cells, 
                                                          split.by]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop(paste0("Need to specify at least ", length(x = unique(x = splits)), 
                    " colors using the cols parameter"))
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = log1p(data.use))
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- unlist(x = lapply(X = data.plot$id, FUN = function(x) sub(paste0(".*_(", 
                                                                                   paste(sort(unique(x = splits), decreasing = TRUE), 
                                                                                         collapse = "|"), ")$"), "\\1", x)))
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                           units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(data.plot)
}

run_custom_GSEA <- function(de_dta, TERM2GENE) {
  
  ### run_custom_GSEA
  ###
  ### This function performs gene set enrichment analysis (GSEA) using a ranked list 
  ### of log2 fold changes from a differential expression dataset.
  ### It removes genes with missing or duplicated symbols, ranks genes by Log2FC, 
  ### and runs the GSEA algorithm using a custom TERM2GENE annotation.
  ###
  ### Args:
  ###   de_dta: A data frame containing at least 'SYMBOL' and 'Log2FC' columns.
  ###   TERM2GENE: A data frame with two columns: 'term' and 'gene', specifying gene sets.
  ###
  ### Returns:
  ###   A GSEA result object from the clusterProfiler::GSEA function.
  
  foldchanges <- de_dta$Log2FC
  names(foldchanges) <- de_dta$SYMBOL
  foldchanges <- sort(foldchanges, decreasing = TRUE)
  foldchanges <- foldchanges[!is.na(names(foldchanges))]
  dups <- names(foldchanges)[duplicated(names(foldchanges))]
  foldchanges <- foldchanges[!names(foldchanges)%in%dups]
  
  GSEA_res <- GSEA(
    geneList = foldchanges,
    pvalueCutoff = 1,
    verbose = T,
    TERM2GENE = TERM2GENE,
    nPermSimple = 10000,
    eps = 0
  )
  
  return(GSEA_res)
}

first_upper <- function(x) {
  
  ### first_upper
  ###
  ### This function converts a character string to sentence case by capitalizing 
  ### the first letter and converting the rest to lowercase.
  ###
  ### Args:
  ###   x: A character string.
  ###
  ### Returns:
  ###   A character string with the first letter capitalized and the rest in lowercase.
  
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

plot_tonic_ISGs <- function(de.results) {
  
  ### plot_tonic_ISGs
  ###
  ### This function creates a customized volcano plot to visualize differential expression 
  ### results for interferon-stimulated genes (ISGs), highlighting tonic IFN-sensitive genes.
  ### Genes are colored by group ("Tonic-sensitive", "Tonic-insensitive", or "Unknown"),
  ### and significant tonic-sensitive genes are labeled.
  ###
  ### Args:
  ###   de.results: A data frame containing differential expression results with columns 
  ###               'logFC', 'FDR', 'Gene.name', and 'Group'.
  ###
  ### Returns:
  ###   A ggplot object showing the volcano plot of ISG expression changes.
  
  de.results <- de.results %>%
    mutate(logP = -log10(FDR)) %>%
    mutate(Label = ifelse(Group=="Tonic-sensitive" & FDR<0.05, Gene.name, NA))
  
  col_vec <- c("#E41A1C",
               "#d3d3d3",
               "#377EB8")
  names(col_vec) <- c("Tonic-sensitive", "Unknown", "Tonic-insensitive")
  
  
  lab_df = de.results %>%
    filter(!is.na(Label))
  
  p <- ggplot(de.results %>% arrange(desc(Group)),
              aes(x=logFC, y=logP, color=factor(Group))) +
    geom_point(size = 0.4) +
    coord_cartesian(xlim=c(-5.5, 3), ylim = c(0, 7.5)) +
    geom_vline(xintercept = log2(2), lty = 2) +
    geom_vline(xintercept = -log2(2), lty = 2) +
    geom_hline(yintercept = -log10(0.05), lty = 2) +
    geom_text_repel(data = lab_df,
                    aes(x = logFC, y = logP, label = Label),
                    min.segment.length = .5,
                    seed = 3,
                    box.padding = .5,
                    show.legend = FALSE,
                    max.overlaps =10,
                    size=5
    ) +
    scale_color_manual(values = col_vec[names(col_vec) %in% unique(de.results$Group)]) +
    theme_cowplot(16) +
    theme(legend.position = "none") +
    labs(x=expression(log[2]*"FC"),
         y = expression(-log[10]*"FDR"))
  
  return(p)
}

plot_scatter <- function(Cell_type) {
  
  ### plot_scatter
  ###
  ### This function creates a log-log scatter plot comparing average gene expression 
  ### between control and KO conditions for a given cell type. Differentially expressed 
  ### genes are color-coded as upregulated, downregulated, or not changed. ISGs that 
  ### are significantly differentially expressed are labeled.
  ###
  ### Args:
  ###   Cell_type: A character string specifying the cell type to subset and analyze.
  ###
  ### Returns:
  ###   A ggplot object showing a scatter plot of gene expression in control vs. KO cells.
  
  de.results_sub <- de.results %>% filter(Cell.type == Cell_type)
  
  ave_exprs <- AverageExpression(
    seurat_obj[,seurat_obj$cell.type == Cell_type],
    assays = "RNA",
    features = de.results_sub$Gene.name,
    slot = "data",
    group.by = "genotype"
  )$RNA
  
  # Extract normalized counts for only the significant genes
  sig_norm <- ave_exprs %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene.name") %>%
    left_join(de.results_sub %>%
                dplyr::select(Gene.name, FDR, logFC), by="Gene.name") %>%
    mutate(DE = if_else(FDR < .05 & logFC > 0, "Upregulated", "Not changed"),
           DE = if_else(FDR < .05 & logFC < 0, "Downregulated", DE)) %>%
    left_join(rat_isgs, by=c("Gene.name"="Rat_symbol")) %>%
    mutate(Label = ifelse(ISG==1 & FDR<0.05, Gene.name, ""))
  
  
  col_vec <- c(as.vector(met.brewer("Nizami", 5)[4]),
               "#d3d3d3",
               as.vector(met.brewer("Nizami", 5)[2]))
  names(col_vec) <- c("Downregulated", "Not changed", "Upregulated")
  
  xy_max <- max(c(sig_norm$HE, sig_norm$KO)) + 40
  xy_min <- min(c(sig_norm$HE, sig_norm$KO)) - .01
  
  p <- ggplot(sig_norm %>% arrange(desc(FDR)), aes(x=HE, y=KO, color=DE, label=Label)) +
    geom_point(size=0.5) +
    scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(xlim=c(xy_min, xy_max), ylim=c(xy_min, xy_max)) +
    geom_text_repel(min.segment.length = 1, seed = 3,
                    box.padding = 1,
                    show.legend = FALSE,
                    max.overlaps =500,
                    force = 3,
                    size=4
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_manual(values = col_vec[names(col_vec) %in% unique(sig_norm$DE)]) +
    theme_cowplot(16) +
    labs(x="Average expr., Control",
         y="Average expr., KO") +
    ggtitle(Cell_type) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p)
}

plot_volcano_color_ISG <- function(cell_type, save_path) {
  
  ### plot_volcano_color_ISG
  ###
  ### This function generates a volcano plot of differential expression results for a 
  ### specific cell type, highlighting interferon-stimulated genes (ISGs) using a binary 
  ### color scheme. Downregulated ISGs that are significantly differentially expressed 
  ### are labeled.
  ###
  ### Args:
  ###   cell_type: A character string specifying the cell type to plot.
  ###   save_path: A placeholder argument (not used within the function as written).
  ###
  ### Returns:
  ###   A ggplot object showing a volcano plot with ISGs highlighted and labeled.
  
  de.results_sub <- de.results %>%
    filter(Cell.type == cell_type) %>%
    mutate(logP = -log10(FDR)) %>%
    mutate(Label = ifelse(ISG==1 & FDR<0.05, Gene.name, "")) %>%
    mutate(DE = if_else(FDR < .05 & logFC > 0, "Upregulated", "Not changed"),
           DE = if_else(FDR < .05 & logFC < 0, "Downregulated", DE))
  
  col_vec <- c("#d3d3d3",
               "black")
  
  names(col_vec) <- c("0", "1")
  
  
  lab_df = de.results_sub %>%
    filter(ISG == 1, DE == "Downregulated")
  
  
  p <- ggplot(de.results_sub %>%
                arrange((ISG)),
              aes(x=logFC, y=logP, color=factor(ISG))) +
    geom_point(size = 0.4) +
    coord_cartesian(xlim=c(-5, 5), ylim = c(0, 8)) +
    geom_vline(xintercept = log2(2), lty = 2, col="orange") +
    geom_vline(xintercept = -log2(2), lty = 2, col="orange") +
    geom_hline(yintercept = -log10(0.05), lty = 2, col="orange") +
    geom_text_repel(data = lab_df,
                    aes(x = logFC, y = logP, label = Label),
                    min.segment.length = .5,
                    seed = 3,
                    box.padding = .5,
                    show.legend = FALSE,
                    max.overlaps =10,
                    size=5
    ) +
    scale_color_manual(values = col_vec[names(col_vec) %in% unique(de.results_sub$ISG)]) +
    theme_cowplot(16) +
    theme(legend.position = "none") +
    labs(x=expression(log[2]*"FC"),
         y = expression(-log[10]*"FDR")) +
    ggtitle(cell_type)
  
  return(p)
}