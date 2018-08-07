library(Seurat)
library(dplyr)
library(Matrix)
library(enrichR)

CST.data <- Read10X(data.dir = "/data2/hej/R/RPE/hg19")
#?鿴ϸ???ĸ???
ncol(as.matrix(CST.data))
#????ʣ???ڴ?
dense.size <- object.size(as.matrix(CST.data))
sparse.size <- object.size(CST.data)
dense.size/sparse.size

#????Seurat
CST <- CreateSeuratObject(raw.data = CST.data, min.cells = 3, min.genes = 2000, 
                           project = "10X_CST")
#??׼Ԥ????????
#?ʿ?QC????????��??????percent.mito values
mito.genes <- grep(pattern = "^MT-", x = rownames(x = CST@data), value = TRUE)
percent.mito <- Matrix::colSums(CST@raw.data[mito.genes, ])/Matrix::colSums(CST@raw.data)
CST <- AddMetaData(object = CST, metadata = percent.mito, col.name = "percent.mito")
#չʾ????????UMI?????Լ?percent.mito
VlnPlot(object = CST, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave("RPE/results/nGene_nUMI_percent_mito.pdf", 
       height = 4, width = 7,units = "in")

#ɸѡϸ????һ??ʹ??geneplot?鿴ngenes/nUMI/percent.mito֮??????????(from seurat website)
par(mfrow = c(1, 2))
GenePlot(object = CST, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = CST, gene1 = "nUMI", gene2 = "nGene")
CST <- FilterCells(object = CST, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(2000, -Inf), high.thresholds = c(6000, 0.15))
#ɸѡϸ?????????ֱ?????ͼ??????ɸѡ
nGene_nUMI <- CST@meta.data
# boxplot of nUMI
ggplot(nGene_nUMI, aes(x = orig.ident, y = nUMI)) + 
  geom_jitter() + 
  geom_boxplot(fill = "lightblue", color = "lightblue", alpha = 0.5) + 
  geom_violin(alpha = 0, color = "blue")
ggsave("RPE/results/boxplot_nUMI.pdf", 
       width = 6, height = 6, units = "in")

# nGene vs nUMI
ggplot(nGene_nUMI, aes(x = nUMI, y = nGene)) + geom_point() + 
  geom_density_2d(color = 'red')
ggsave("RPE/results/nUMI_vs_nGene.smalldot.pdf",width=6,height=6,units="in")

# nGene distribution
pdf("RPE/results/dist_of_nGene.pdf", 
    width = 7, height = 4)
hist(nGene_nUMI$nGene)
dev.off()
nGene_thres <- unname(quantile(nGene_nUMI[,"nGene"],probs=c(0.95)))

# nUMI distribution
pdf("RPE/results/dist_of_nUMI.pdf", width=7,height=4)
hist(nGene_nUMI$nUMI)
dev.off()
nUMI_thres <- unname(quantile(nGene_nUMI[,"nUMI"],probs=c(0.95)))

# percent.mito distribution
pdf("RPE/results/dist_of_mito.pdf", width = 7,
    height = 4)
hist(nGene_nUMI$percent.mito)
dev.off()
# mito_thres <- unname(quantile(nGene_nUMI[,"percent.mito"],probs=c(0.95)))

# ɸѡϸ??
CST <- FilterCells(object = CST, 
                    subset.names = c("nUMI", "percent.mito"), 
                    low.thresholds = c(10000, -Inf), 
                    high.thresholds = c(25000, 0.12))
#nhne <- FilterCells(object = nhne, 
#                    subset.names = c("nGene","nUMI", "percent.mito"), 
#                    low.thresholds = c(-Inf, 10000, -Inf), 
#                    high.thresholds = c(nGene_thres, nUMI_thres, mito_thres))

# nhne <- FilterCells(object = nhne, 
#                     subset.names = c("nGene","nUMI", "percent.mito"), 
#                     low.thresholds = c(-Inf, 10000, -Inf), 
#                     high.thresholds = c(6000, 45000, 0.10))
#?鿴һ??ɸѡ????ϸ????��
CST@meta.data %>% nrow()
cells_filtered <- colnames(x = CST@data)

#???ݱ?׼??
CST <- NormalizeData(object = CST, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
# scale.data: center the mean expression to 0
CST <- ScaleData(object = CST)

# ???Ҳ???????find high variable genes??ע??????ֵ?????ã?????ʹ?ò???????????��??1000-2000֮?䣩
CST <- FindVariableGenes(object = CST, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 4, 
                          y.cutoff = 0.4, 
                          do.plot = T)
length(CST@var.genes)
#????marker?????Ƿ??ڲ?????????(??ѡ)
#Mitif RPE65 Bestrophin Rlbp1
c("MITF", "BEST1", "POU5F1", "RLBP1", "RPE65") %in% CST@var.genes
#????mito.gene?Ƿ??ڲ?????????
length(CST@var.genes[! CST@var.genes %in% mito.genes])
# ȷ???߱??????򣨳?ȥ??��????????
hvg_seurat <- CST@var.genes[! CST@var.genes %in% mito.genes]

# ??????????ϸ???????еķ???
rowVar_expMat <- apply(as.matrix(CST@data), 1, var)
hist(rowVar_expMat, nclass = 20)
# ??????????5000??????
VarTop5000 <- head(names(rowVar_expMat[order(rowVar_expMat, 
                                             decreasing=TRUE)]), n = 5000)
VarTop5000 <- VarTop5000[!VarTop5000 %in% mito.genes]

#PCA????????
CST <- RunPCA(CST, pc.genes = hvg_seurat,
               pcs.compute = 40, do.print = F)
pdf("40PCA_genes.pdf", height = 40, width = 7)
PCHeatmap(object = CST, pc.use = 1:40, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()
# ʹ??Jackstraw???????????ɷ?: Significant PCs 
CST <- JackStraw(object = CST, num.pc = 20,
                  num.replicate = 100, do.print = FALSE)
pdf(file = "JackStraw_plot.pdf", height = 30, width = 7)
JackStrawPlot(object = CST, PCs = 1:20)
dev.off()
# ElbowPlot ??ʾȡ???ɷֵĸ???
pdf(file = "PCElbowPlot.pdf", height = 4, width = 9)
PCElbowPlot(object = CST, num.pc = 20)
dev.off()

# ѡ?????ɷ?: Cluster with seurat (time-consuming)????��ȷ?????ɷָ????????Ѱ취
NumCluster <- rep(0,40)
for (k in 2:40) {
CST <- FindClusters(object = CST, reduction.type = "pca", dims.use = 1:k, 
                    print.output = 0)
NumCluster[k] = length(levels(CST@ident))
  }
ggplot(data.frame(PC_used = 1:40, NumCluster),aes(x = PC_used,y = NumCluster)) + 
   geom_point() + scale_y_continuous(breaks = seq(0,14,2))
ggsave("NumCluster_seurat_vs_PCused.pdf",width = 8, height = 4, units = "in")

# ʹ??tSNE with seurat????
CST <- FindClusters(object = CST, reduction.type = "pca", dims.use = 1:20, 
                     print.output = FALSE, force.recalc = TRUE,, ave.SNN = TRUE),
                    temp.file.location = '/data2/hej/R/RPE/results')CST <- RunTSNE(object = CST, dims.use = 1:20, do.fast = TRUE)

comp <- data.frame(CST@dr$tsne@cell.embeddings,
                   cluster = CST@ident,
                   nUMI = nGene_nUMI[cells_filtered, "nUMI"],
                   nGene = nGene_nUMI[cells_filtered, "nGene"],
                   pct_mito = nGene_nUMI[cells_filtered, "percent.mito"])

TSNEPlot(object = CST, 
         colors.use = c("darkorange","red","navy",
                        "darkgreen","olivedrab1","cyan", "sienna",
                        "pink"))
#length(comp$cluster[comp$cluster == 5])

CST <- RunPCA(object = CST)
PCAPlot(object = CST)

# Feature plot

FeaturePlot(object = CST, 
            features.plot = c('MITF', 'BEST1', 'RLBP1', 'RPE65', 
                              'POU5F1', 'NANOG','MYBL2' , 'nGene', 'nUMI'), 
            cols.use = c("grey", "red"), 
            reduction.use = "tsne")


# Find markers
CST.markers <- FindAllMarkers(object = CST, only.pos = TRUE, min.pct = 0.2)
 write.table(CST.markers, "cluster.markers.seurat.hvg.PC20.csv",
            sep = ",", row.names = F, col.names = T)

# maker GO-BP
 dbs <- c("GO_Biological_Process_2017")
 for (i in unique(as.character(CST.markers$cluster))) {
   enriched <- enrichr(CST.markers[CST.markers$cluster==i,"gene"],dbs)
   x <- enriched[["GO_Biological_Process_2017"]]
   write.table(x,
               paste0("cluster_",i,".enrichr.csv"),
               sep=",", row.names=FALSE
   )
 }
# enrich plot
cluster77_enrich <- read.csv(file = '/data2/hej/R/RPE/cluster_7.enrichr.csv',nrows = 15)
library(reshape2)
tmp7 <- colsplit(cluster7_enrich$Overlap, '/', c('count', 'total'))
cluster7_enrich$counts <- tmp7count
library(graphics)
cluster7_enrich$Term <- reorder(cluster7_enrich$Term, rev(cluster7_enrich$P.value), median)
 ggplot(data = cluster7_enrich, aes(x = Term, y = counts, fill = P.value)) + 
   geom_bar(stat = 'identity', position = 'stack') + coord_flip() +
   scale_fill_gradientn(colours = colorRampPalette(c("red","blue"))(100)) +
   theme(axis.text.y = element_text(size = 8)) 
 
 
# preview
CST.markers %>% group_by(cluster) %>% top_n(2, avg_logFC) %>% ungroup() -> tmp

FeaturePlot(object = CST, features.plot = c('DCT', 'TFPI2'), 
            cols.use = c("grey", "olivedrab1"), reduction.use = "tsne")
ggsave("/data2/hej/R/RPE/results/top2_marker_genes.pdf", width=12, height=12,units="in")

FeaturePlot(object = nhne,
            features.plot = c("Snai1", "Cxcl12",
                              "Igf1","Prrx1",
                              "Sox2", "Nes", "Pax6","Notch1",
                              "Prrx2","Snai2","Plat",
                              "Olig2",
                              "Dcx","Rbfox3","Th","Chat","Nrg1",
                              "Gdnf",
                              "Bmp4","Hand1","Gata6",
                              "Tgfb2","Tgfb1","Kitl","Gata4",
                              "Ascl1","Sox10","Erbb3",
                              "Epcam","Gata3","Bmp7","Wnt4",
                              "Pou5f1","Dppa3","Esrrb","Kit",
                              "Esam","Cdh5","Kdr"),
            cols.use = c("grey", "red"),
            reduction.use = "tsne")
ggsave("celltype_spec_marker_genes.pdf", width=8, height=20,units="in")


# SCF corrleated genes
BarPlot <- function(object,genes.plot,cells.plot,clusters.plot) {
  if (! is.null(clusters.plot)) {
    x = object@ident
    cells.plot = names(x[x %in% clusters.plot])
  }
  data = reshape2::melt(mtx)
  colnames(data) = c("gene","cell","expression")
  if (! is.null(clusters.plot)) {
    x = object@ident
    data$cluster = x[as.character(data$cell)]
    ggplot(data,aes(x=cell,y=expression,fill=cluster)) + 
      geom_bar(stat="identity") + 
      facet_grid(gene~cluster) + 
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            strip.text.y=element_text(angle=0))
  }
  else {
    ggplot(data,aes(x=cell,y=expression)) + 
      geom_bar(stat="identity") + 
      facet_grid(gene~.) + 
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            strip.text.y=element_text(angle=0))
  }
}

BarPlot(nCSTgenes.plot=c("Wt1","Kitl"),clusters.plot=c("5","6"))
ggsave("Wt1_SCF.pdf",width=12,height=4,units="in")


save(nhne, nGene_nUMI, cells_filtered, hvg_seurat, mito.genes, file="1537.RData")











