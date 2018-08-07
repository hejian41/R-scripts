########## SC3 ############
library(SC3)
library(SingleCellExperiment)
library(scater)
sce <- SingleCellExperiment(assay = list(counts = as.matrix(seu@raw.data[, rownames(seu@meta.data)]),
                                         logcounts = log2(as.matrix(seu@raw.data[, rownames(seu@meta.data)]) + 1)))

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
## RUN pca, tsne, diffusionmap, MDS
sce <- runPCA(sce)
plotPCA(sce)
sce <- runTSNE(sce)
plotTSNE(sce)
sce <- runDiffusionMap(sce)
plotDiffusionMap(sce)
sce <- runMDS(sce)
plotMDS(sce)
######
sce <- sc3_estimate_k(sce)
sce@metadata$sc3$k_estimation
sce <- sc3(sce, ks = 5:7, biology = TRUE)
cluster <- as.data.frame(sce@colData[, "sc3_6_clusters", drop = F])

plotPCA(sce, colour_by = "sc3_6_clusters")
plotTSNE(sce, colour_by = "sc3_5_clusters")

###ggplot showing
pc.frame <- as.data.frame(sce@reducedDims$PCA)
ggdata <- data.frame(row.names = rownames(sce@colData), cluster = sce@colData$sc3_6_clusters, stage = seu@meta.data$stage,
                     Dim1 = pc.frame$PC1, Dim2 = pc.frame$PC2)
ggplot(ggdata, aes(x = Dim1, y = Dim2, color = stage))+
  geom_point(size = 3) + scale_color_manual(values = pal_dolphin)

sc3_plot_consensus(sce, k = 2)
sc3_plot_consensus(sce, k = 4, show_pdata = "sc3_4_clusters")
sc3_plot_silhouette(sce, k = 4)
sc3_plot_expression(sce, k = 2, show_pdata = "sc3_2_clusters")
sc3_plot_cluster_stability(sce, k = 2)
sc3_plot_de_genes(sce, k = 2, show_pdata = "sc3_2_clusters", p.val = 0.01)
sc3_plot_markers(sce, k = 2, auroc = 0.85, p.val = 0.01,
                 show_pdata = "sc3_2_clusters")
cluster_sc3 <- cluster[rownames(seu@meta.data) ,, drop = F]
seu@meta.data$cluster <- cluster_sc3$sc3_4_clusters
