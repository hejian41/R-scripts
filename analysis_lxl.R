library(Seurat)
library(dplyr)
library(Matrix)
library(enrichR)


nhne.data <- Read10X(data.dir="filtered_gene_bc_matrices/mm10")
ncol(as.matrix(nhne.data))

nhne <- CreateSeuratObject(raw.data = nhne.data, 
                           min.cells = 3, min.genes = 2000, 
                           project = "10X_nhne")

# evaluate threshold for percent.mito
mito.genes <- grep(pattern = "^mt-", x = rownames(x = nhne@data), 
                   value = TRUE)
percent.mito <- Matrix::colSums(nhne@raw.data[mito.genes, ]) / 
                Matrix::colSums(nhne@raw.data)

nhne <- AddMetaData(object = nhne, metadata = percent.mito, 
                    col.name = "percent.mito")

VlnPlot(object = nhne, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)
ggsave("./filter/nGene_nUMI_percent_mito.pdf", height=4, width=7,units="in")


# seurat
# evaluate the nUMI/nGene threshold
# VlnPlot(object = nhne, features.plot = c("nGene", "nUMI"), nCol = 2)
# GenePlot(object = nhne, gene1 = "nUMI", gene2 = "nGene")

nGene_nUMI <- nhne@meta.data
# boxplot of nUMI
ggplot(nGene_nUMI,aes(x=orig.ident,y=nUMI)) + 
  geom_jitter() + 
  geom_boxplot(aes(fill="orange",color="orange",alpha=0.5)) + 
  geom_violin(alpha=0,color="green")
ggsave("./filter/boxplot_nUMI.pdf",width=6,height=6,units="in")

# nGene vs nUMI
ggplot(nGene_nUMI, aes(x=nUMI, y=nGene)) + geom_point() + geom_density_2d()
ggsave("./filter/nUMI_vs_nGene.smalldot.pdf",width=6,height=6,units="in")

# nGene distribution
pdf("./filter/dist_of_nGene.pdf", width=7,height=4)
hist(nGene_nUMI$nGene)
dev.off()
nGene_thres <- unname(quantile(nGene_nUMI[,"nGene"],probs=c(0.95)))

# nUMI distribution
pdf("./filter/dist_of_nUMI.pdf", width=7,height=4)
hist(nGene_nUMI$nUMI)
dev.off()
nUMI_thres <- unname(quantile(nGene_nUMI[,"nUMI"],probs=c(0.95)))

# percent.mito distribution
pdf("./filter/dist_of_mito.pdf", width=7,height=4)
hist(nGene_nUMI$percent.mito)
dev.off()
# mito_thres <- unname(quantile(nGene_nUMI[,"percent.mito"],probs=c(0.95)))

# filter cells
nhne <- FilterCells(object = nhne, 
                    subset.names = c("nUMI", "percent.mito"), 
                    low.thresholds = c(10000, -Inf), 
                    high.thresholds = c(50000, 0.10))
#nhne <- FilterCells(object = nhne, 
#                    subset.names = c("nGene","nUMI", "percent.mito"), 
#                    low.thresholds = c(-Inf, 10000, -Inf), 
#                    high.thresholds = c(nGene_thres, nUMI_thres, mito_thres))

# nhne <- FilterCells(object = nhne, 
#                     subset.names = c("nGene","nUMI", "percent.mito"), 
#                     low.thresholds = c(-Inf, 10000, -Inf), 
#                     high.thresholds = c(6000, 45000, 0.10))
nhne@meta.data %>% nrow()
# nGene            nUMI         percent.mito
# low    high    low    high    low    high    nFilteredCells
# -      -       10000  50000   -Inf   0.10    1537
# -      6238    10000  49907   -Inf   0.047   1483
# -      6000    10000  45000   -Inf   0.10    1470
cells_filtered <- colnames(x = nhne@data)

# a global-scaling normalization method “LogNormalize” that normalizes the 
# gene expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), 
# and log-transforms the result.
nhne <- NormalizeData(object = nhne, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# nhne@scale.data: center the mean expression to 0
nhne <- ScaleData(object=nhne)

# find high variable genes
# FindVariableGenes(object = nhne, 
#                  x.low.cutoff = 0.0125, x.high.cutoff = 3, 
#                  y.cutoff = 0.5, 
#                  do.plot = TRUE)

nhne <- FindVariableGenes(object = nhne, 
                          x.low.cutoff = 0.01, x.high.cutoff = 3.5, 
                          y.cutoff = 0.5, 
                          do.plot = FALSE)
length(x = nhne@var.genes)
# parameters result
# nVar.Genes    x.low.cutoff    x.high.cutoff    y.cutoff
# 1806/1791     0.01            3.5              0.5
# 1610          0.0125          3                0.5

# check if some marker genes were included
c("Pou5f1", "Nanog", "Dppa3", "Wt1", "Cdh5", "Kitl") %in% nhne@var.genes
# check var.genes without mitochondrial genes
length(nhne@var.genes[! nhne@var.genes %in% mito.genes])
# high variable genes
hvg_seurat <- nhne@var.genes[! nhne@var.genes %in% mito.genes]

# variations accross cells for each genes
# rowVar_expMat <- apply(as.matrix(nhne@data), 1, var)
# hist(rowVar_expMat, nclass=20)
# VarTop5000 <- head(
#   names(rowVar_expMat[order(rowVar_expMat, decreasing=TRUE)]), 
#   n = 5000)
# VarTop5000 <- VarTop5000[!VarTop5000 %in% mito.genes]


# PCA
nhne <- RunPCA(nhne, pc.genes = hvg_seurat,
               pcs.compute=40, do.print = FALSE)

pdf("PCA_genes.pdf", height=40, width=7)
PCHeatmap(object = nhne, pc.use = 1:40, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()

# PCA: Significant PCs 
JackStraw
nhne <- JackStraw(object = nhne, num.pc = 40,
                 num.replicate = 100, do.print = FALSE)
pdf(file = "JackStraw_plot.pdf",height=30, width=7)
JackStrawPlot(object = nhne, PCs = 1:40)
dev.off()

# ElbowPlot
pdf(file = "PCElbowPlot.pdf",height=4,width=9)
PCElbowPlot(object = nhne, num.pc = 40)
dev.off()

# dim selection: Cluster with seurat (time-consuming)
# NumCluster <- rep(0,40)
#for (k in 2:40) {
#nhne <- FindClusters(object = nhne, reduction.type = "pca", dims.use = 1:k, 
#                    print.output = 0)
  
#NumCluster[k] = length(levels(nhne@ident))
#}
#ggplot(data.frame(PC_used=1:40,NumCluster),aes(x=PC_used,y=NumCluster)) + 
#  geom_point() + scale_y_continuous(breaks = seq(0,14,2))
#ggsave("NumCluster_seurat_vs_PCused.pdf",width=8,height=4,units="in")



# tSNE with seurat
nhne <- FindClusters(object = nhne, reduction.type = "pca", dims.use = 1:20, 
                     print.output = FALSE, force.recalc = TRUE, save.SNN = TRUE)
nhne <- RunTSNE(object = nhne, dims.use = 1:20, do.fast = TRUE)

comp <- data.frame(nhne@dr$tsne@cell.embeddings,
                   cluster=nhne@ident,
                   nUMI=nGene_nUMI[cells_filtered,"nUMI"],
                   nGene=nGene_nUMI[cells_filtered,"nGene"],
                   pct_mito=nGene_nUMI[cells_filtered,"percent.mito"])

TSNEPlot(object = nhne, 
         colors.use = c("darkorange","red","navy",
                        "darkgreen","olivedrab1","cyan",
                        "sienna","grey50","magenta4", 
                        "pink","steelblue1"))

# Feature plot
FeaturePlot(object = nhne, 
            features.plot = c("Dppa3","Nanog","Cdh5","Esam",
                              "Kit","Kitl","nGene","nUMI"), 
            cols.use = c("grey", "red"), 
            reduction.use = "tsne")


# Find markers
nhne.markers <- FindAllMarkers(object = nhne, only.pos = TRUE, 
                               min.pct = 0.4,
                               min.diff.pct = 0.3)
write.table(nhne.markers,"cluster.markers.seurat.hvg.20PC.csv",
            sep=",",row.names=F,col.names=T)

# maker GO-BP
dbs <- c("GO_Biological_Process_2017")
for (i in unique(as.character(nhne.markers$cluster))) {
  enriched <- enrichr(nhne.markers[nhne.markers$cluster==i,"gene"],dbs)
  x <- enriched[["GO_Biological_Process_2017"]]
  write.table(x,
              paste0("cluster_dim20/cluster_",i,".enrichr.csv"),
              sep=",", row.names=FALSE
  )
}


# preview
nhne.markers %>% group_by(cluster) %>% top_n(2, avg_logFC) %>% ungroup() -> tmp

FeaturePlot(object = nhne,
  features.plot = tmp$gene,
  cols.use = c("grey", "red"),
  reduction.use = "tsne")
ggsave("top2_marker_genes.pdf", width=8, height=12,units="in")

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

BarPlot(nhne,genes.plot=c("Wt1","Kitl"),clusters.plot=c("5","6"))
ggsave("Wt1_SCF.pdf",width=12,height=4,units="in")


save(nhne, nGene_nUMI, cells_filtered, hvg_seurat, mito.genes, file="1537.RData")


