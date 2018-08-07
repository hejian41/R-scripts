##### monocle #######
library(monocle)
pd <- new("AnnotatedDataFrame", data = cbind(seu@meta.data,t(seu@data)))
# HSMM <- newCellDataSet(as.matrix(seu@raw.data[, rownames(seu@meta.data)]), phenoData = pd, expressionFamily = VGAM::negbinomial.size() )
HSMM <- newCellDataSet(as.matrix(seu@data[, rownames(seu@meta.data)]), phenoData = pd, expressionFamily = VGAM::negbinomial.size() )
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table, 
                         # mean_expression >= 0.5 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id
ordering_genes <- setdiff(ordering_genes, cellcycle_all)
HSMM<- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM, reverse = FALSE)

plot_cell_trajectory(HSMM, color_by = "CD44", cell_size = 2, do.return = T) + scale_color_gradientn(colours = colorRampPalette(c("blue","green","yellow","red"))(100))
plot_cell_trajectory(HSMM, color_by = "State", cell_size = 2, do.return = T)
plot_cell_trajectory(HSMM, color_by = "cluster", cell_size = 2, do.return = T) + scale_color_manual(values = ph_colors$cluster)
