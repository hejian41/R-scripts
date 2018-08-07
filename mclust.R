###### mclust #######
library(mclust)
library(factoextra)
M <- Mclust(seu@dr$pca@cell.embeddings[, 1:10])
# M <- Mclust(t(seu@data[seu@var.genes, ]))
fviz_mclust_bic(M)
Mdr <- MclustDR(M)
plot(Mdr, what = "boundaries")

ph1_colors <- list(stage = setNames(myrainbow[1:length(unique(seu@meta.data$stage))], sort(unique(seu@meta.data$stage))),
                   mclust = setNames(pal_dolphin[1:length(unique(seu@meta.data$mclust))], sort(unique(seu@meta.data$mclust))),
                   location = setNames(bertie.color[1:length(unique(seu@meta.data$location))], sort(unique(seu@meta.data$location))))
seu@meta.data$mclust <- factor(M$classification)
DimPlot(seu, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, 
        pt.size = 2, cols.use = pal_dolphin, group.by = "mclust")
