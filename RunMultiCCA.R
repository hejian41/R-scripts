options(stringsAsFactors = F)
setwd("")
rm(list = ls())

load("D:/data/human HSC analysis//human HSC/bzj20180427/cellcyclegenes.RData")
load("D:/data/human HSC analysis//human HSC/bzj20180427/labcolors.RData")

#####
# data_umi <- cbind(data_raw0330, data_raw0411, data_raw170228, data_raw170324, data_raw20180625, data_raw160707, data_raw170421)
# data_annot <- rbind(data_annot0330[, c("stage", "location")], data_annot0411[, c("stage", "location")], 
#                     data_annot170228[, c("stage", "location")], data_annot170324[, c("stage", "location")],
#                     data_annot20180625[, c("stage", "location")], data_annot160707[, c("stage", "location")],
#                     data_annnot170421[, c("stage", "location")])
# save(data_umi, data_annot, file = "")

load(file = "CS10_16_Rawdata_Annotation.RData")
load(file = "CS10_16_annot.RData")

data_annot <- subset(CS10_16_annot, Cluster == "EC_IE" & location %in% c("Body", "DA"))
data_umi <- data_umi[, rownames(data_annot)]
table(data_annot$stage, data_annot$location)

ERCC <- base::grep("ERCC-|RGC-|RPS|RPL", rownames(data_umi), value = F)
data_umi <- data_umi[-ERCC,]
####
#####
library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)


for (i in unique(data_annot$stage)) {
  seu_name <- paste("seu", i, sep="_")
  assign(seu_name, CreateSeuratObject(data_umi[, rownames(subset(data_annot, stage == i))], 
                                      meta.data = subset(data_annot, stage == i)))
}
  
seu.list <- list(seu_CS10, seu_CS11, seu_CS12, seu_CS13, seu_CS14, seu_CS15, seu_CS16)
seu.list <- lapply(seu.list, function(x){
  x = NormalizeData(x)
  cc.genes <- readLines(con = "D:/data/bzj2018human 小胚数据/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
  cc.genes <- toupper(cc.genes)
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]
  x <- CellCycleScoring(x, g2m.genes = g2m.genes, s.genes = s.genes, set.ident = T)
  x <- ScaleData(x, vars.to.regress = c("S.Score", "G2M.Score")) 
  x = FindVariableGenes(x)
})

hvg.union <- c()
for (i in 1:length(seu.list)) {
  hvg.union = c(hvg.union, head(rownames(seu.list[[i]]@hvg.info), 1000))
  hvg.union = unique(hvg.union)
}


seu.cca <- RunMultiCCA(object = seu.list, genes.use = hvg.union, num.ccs = 15)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = seu.cca, reduction.use = "cca", group.by = "stage", pt.size = 2, 
              do.return = TRUE)
p2 <- VlnPlot(object = seu.cca, features.plot = "CC1", group.by = "stage", do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = seu.cca, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
# choose CCs for downstream analysis and then ‘align them’
# here explore the CC dimensions as we have previously demonstrated for PCA.
DimHeatmap(object = seu.cca, reduction.type = "cca",  dim.use = 1:9, 
           do.balanced = TRUE)
DimHeatmap(object = seu.cca, reduction.type = "cca",  dim.use = 10:15, 
           do.balanced = TRUE)

# Before we align the subspaces, we first search for cells whose expression profile cannot \
# be well-explained by low-dimensional CCA, compared to low-dimensional PCA.
# seu.cca <- CalcVarExpRatio(object = seu.cca, reduction.type = "pca", grouping.var = "stage", dims.use = 1:9)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
seu.all.save <- seu.cca
# seu.cca <- SubsetData(object = seu.cca, subset.name = "var.ratio.pca", accept.low = 0.1)

# Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
seu.cca <- AlignSubspace(object = seu.cca, reduction.type = "cca", grouping.var = "stage", 
                      dims.align = 1:15)
# Visualize the aligned CCA and perform integrated analysis
p1 <- VlnPlot(object = seu.cca, features.plot = "ACC1", group.by = "stage", 
              do.return = TRUE)
p2 <- VlnPlot(object = seu.cca, features.plot = "ACC2", group.by = "stage", 
              do.return = TRUE)
plot_grid(p1, p2)

gene_matrix <- as.data.frame(as.matrix(seu.cca@data))
save(gene_matrix, file = "gene_matrix0-0.RData")

# Now we can run a single integrated analysis on all cells!
seu.cca <- RunPCA(seu.cca)
seu.cca <- ProjectPCA(seu.cca)
DimElbowPlot(seu.cca)
PCHeatmap(seu.cca, pc.use = 1:9)
PCHeatmap(seu.cca, pc.use = 10:18)
DimPlot(seu.cca,  pt.size = 2, dim.1 = 1, dim.2 = 3, reduction.use = "pca", cols.use = pal_dolphin, 
        do.return = T, group.by = "stage")
DimPlot(seu.cca,  pt.size = 2, dim.1 = 1, dim.2 = 3, reduction.use = "pca", cols.use = pal_dolphin, 
        do.return = T, group.by = "location")
seu.cca <- RunTSNE(seu.cca, reduction.use = "pca", dims.use = 1:10)
DimPlot(seu.cca,  pt.size = 2, dim.1 = 1, dim.2 = 2, reduction.use = "tsne", cols.use = pal_dolphin, 
        do.return = T, group.by = "stage")

myFeaturePlot(seu.cca, features.plot = c("CD34", "PECAM1", "CDH5", "KDR"), dr = "tsne", ncol = 4)
myFeaturePlot(seu.cca, features.plot = c("EFNB2","DLL4", "GJA5", "EPHB4","NR2F2", "NRP2","PROX1", "LYVE1"), dr = "pca", ncol = 4)
myFeaturePlot(seu.cca, features.plot = c("RUNX1", "SOX17", "CD44", "NOTCH1","GATA1", "PTPRC","SPN", "GFI1", "ACE", "TAL1"), dr = "pca", ncol = 4)
