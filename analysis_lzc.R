options(stringsAsFactors = F)
setwd("/data2/hej/R/bzj20180328/rawdata/")
################################
# Prepare initial data
################################
# umi <- read.table(file = "refGene.UMI.count.xls", header = T, sep = "\t",row.names = 1)
# stat <- read.table(file = "qc_align_gene.stat.xls", header = T, sep = "\t", row.names = 1)
# tmp <- reshape2::colsplit(colnames(umi),"_sc", c("Library","Barcode"))
# unique(tmp$Library)
# 
# libs <- unique(tmp$Library)
# select_cells <- c(
#   paste0(rep(libs[1], each = 13), "_sc",
#          rep(1:13, times = 1)),
#   paste0(rep(libs[-1], each = 48), "_sc",
#          rep(setdiff(1:50,c(23,34)), times = 6))
# )
# 
# data_umi1 <- umi[,select_cells]
# data_stat1 <- stat[select_cells,]
# data_annot1 <- reshape2::colsplit(colnames(data_umi1), "_", c("Type",2:4))[,1,drop = F]
# rownames(data_annot1) <- colnames(data_umi1)
# data_annot1$Type <- gsub("[0-9]", replacement = "1", data_annot1$Type)
data_annot1$Type[1:13] <- paste0("BH","1")
save(data_umi1, data_stat1, data_annot1, file = "data1.Rdata")

load("/data2/hej/R/bzj20180328/rawdata/data1.Rdata")
library(Seurat)
library(Rtsne)
library(dplyr)
library(pheatmap)

# seurat

pbmc <- CreateSeuratObject(raw.data = data_umi[setdiff(rownames(data_umi), ""),], meta.data = data_annot,
                           min.cells = 0, min.genes = 0, is.expr = 0, 
                           normalization.method = NULL, scale.factor = 1e6, do.scale = F, do.center = F, 
                           names.field = 1)

pdf(file = "nUMI-nGene.pdf", width = 8, height = 8)
plot(pbmc@meta.data$nUMI, pbmc@meta.data$nGene, col = pal_dolphin[c(1:7)][as.factor(pbmc@meta.data$Type)], 
     pch = 19, cex=1.6, xlim = c(0,1e6), ylim = c(0, 9000), 
     xlab = "nUMI - the number of transcripts", ylab = "nGene - the number of genes")
abline(v = 2.5e4, col = "red"); abline(h = 1000, col = "red")
text(x = 7e5, y = 2000, labels = "nGene > 1000")
text(x = 6e4, y = 8000, labels = "nUMI > 25000", adj = 0)
legend(x = "right", inset = 0, legend = unique(pbmc@meta.data$Type), col = pal_dolphin[c(1:7)][as.factor(unique(pbmc@meta.data$Type))],
       pch = 19, cex= 0.8)
dev.off()
# origin cells: BH BN1 BN3 FH1 FH2 FN1 FN2 
#               13  48  48  48  48  48  48   total  301 
pbmc <- FilterCells(pbmc, subset.names = "nUMI", low.thresholds = 25000)
pbmc <- FilterCells(pbmc, subset.names = "nGene", low.thresholds = 1000)
# filtered cells : BH BN1 BN3 FH1 FH2 FN1 FN2 
#                  9  44  42  44  25  16  36 total  216

data_annot1 <- data_annot[colnames(pbmc@data), , drop = F]
tmp <- substr(data_annot1$Type, 1, 2)
# tmp <- grep("[a-zA-Z]", data_annot1$Type, value = T)
data_annot_union <- data.frame(row.names = rownames(data_annot1), type = tmp)
pbmc <- AddMetaData(pbmc, metadata = data_annot_union)

pbmc <- SetAllIdent(pbmc, id = 'Type')
pdf("nGene_nUMI_Vlnplot.pdf", width = 10, height = 8)
VlnPlot(pbmc, features.plot = c("nGene","nUMI"), cols.use = pal_dolphin[c(1:7)][as.factor(unique(pbmc@meta.data$Type))], x.lab.rot = T)
dev.off()


# points(pbmc@meta.data$nUMI, pbmc@meta.data$nGene, col = pal_dolphin[as.factor(pbmc@meta.data$Library)], 
#        pch = 19, cex=1.6, lwd = 1.5)
length(pbmc@cell.names)/ncol(data_umi)
pbmc@data <- apply(pbmc@data, 2, function(x) {log2(1e5*x/sum(x) + 1)})
pbmc <- ScaleData(pbmc, do.scale = T, do.center = T)
pbmc <- SetAllIdent(pbmc,id = "Library")
pbmc <- SetAllIdent(pbmc,id = "type")

pdf("featureplot.pdf", width = 10, height = 8)
VlnPlot(pbmc, features.plot = c("PTPRC", "CD34", "CD38", "THY1", "VNN2",  "ACE",
                                "PROCR"), cols.use =  pal_dolphin[c(1:4)], x.lab.rot = T, nCol = 4)
dev.off()

ph_colors <- sapply(colnames(pbmc@meta.data), function(x){
  cats <- sort(unique(pbmc@meta.data[[x]]))
  if(length(cats) <= length(pal_dolphin))
    setNames(pal_dolphin[1:length(cats)], cats)
})

pdf("pheatmap_genes.pdf", width = 8, height = 10)
pheatmap::pheatmap(pbmc@data[c("CD34", "CD9", "KIT", "TFRC", "GYPA", "CD36", "CD38", "PTPRC", 
                               "EPOR", "HBA1", "ACE", "THY1"),],
                   # which(pbmc@meta.data$Cell_type %in% "IVS.1")], 
                   cluster_rows = F, cluster_cols = F, annotation_col = pbmc@meta.data[,"type", drop=F],
                   annotation_colors = ph_colors, gaps_col = cumsum(table(pbmc@meta.data$type)),
                   show_colnames = F, color = gplots::colorpanel(n = 11, "darkblue","white","red"), border_color = NA)
dev.off()


pbmc <- FindVariableGenes(pbmc, x.low.cutoff = 1, x.high.cutoff = Inf, y.cutoff = 1)
length(pbmc@var.genes)
pbmc <- RunPCA(object = pbmc, pc.genes = setdiff(rownames(pbmc@data), cellcycle))
pbmc <- RunPCA(object = pbmc, pc.genes = setdiff(pbmc@var.genes, cellcycle))
pbmc <- ProjectPCA(pbmc,do.print = F)
DimElbowPlot(pbmc, reduction.type = "pca", dims.plot = 20, xlab = "",
             ylab = "", title = "")
pbmc <- FindClusters(pbmc, dims.use = 1:10, k.param = 10, 
                      plot.SNN = T, resolution = 0.8)
PCAPlot(pbmc, do.hover = F, group.by = "type", do.identify = F, 
        cols.use = pal_dolphin[c(1:4)][as.factor(unique(pbmc@meta.data$type))], pt.size =3, do.label = F)

PCAPlot(pbmc, do.hover = F, group.by = "Type", do.identify = F, 
        cols.use = pal_dolphin, pt.size =3, do.label = F)
PCAPlot(pbmc, do.hover = F, group.by = "res.0.8", do.identify = F, 
        cols.use = pal_dolphin, pt.size =3, do.label = F)
pbmc <- RunTSNE(pbmc, perplexity = 5, genes.use = setdiff(pbmc@var.genes,""),dims.use = 1:5)
pbmc <- RunTSNE(pbmc, perplexity = 5, genes.use = setdiff(pbmc@var.genes, cellcycle))

TSNEPlot(pbmc, pt.size = 3, group.by = "type", colors.use = pal_dolphin[c(1:4)][as.factor(unique(pbmc@meta.data$type))], do.label = F)
TSNEPlot(pbmc, pt.size = 3, group.by = "Library", colors.use = pal_dolphin, do.hover = F)

library(Rtsne)
d_euclidean <- stats::dist(t(pbmc@data[setdiff(pbmc@var.genes,cellcycle),]))
d_euclidean <- pairwise.topgene.corr(pbmc@data[setdiff(rownames(pbmc@data), cellcycle),], ngene = 1000)
d_euclidean <- pairwise.topgene.equal(pbmc@data[setdiff(rownames(pbmc@data), cellcycle),], ngene = 1000)

set.seed(10)
tsne_out <- Rtsne(d_euclidean, is_distance=TRUE, pca = F,perplexity=5, verbose = TRUE,max_iter=1500, theta = 0.0) 
#tsne_out <- Rtsne(as.dist((1-cor(pbmc@data))/2), is_distance=TRUE, perplexity=30, verbose = TRUE,max_iter=2500) 
#plot(tsne_out$Y, pch=16, main='tSNE')
pbmc@dr$tsne@cell.embeddings[,c(1,2)] <- tsne_out$Y[,c(1,2)]
pbmc <- FindClusters(pbmc, reduction.type = "tsne", dims.use = 1:2, force.recalc = T, resolution = 0.8, plot.SNN = T)
# pbmc <- FindClusters(pbmc, reduction.type = "pca", dims.use = 1:10, force.recalc = T, resolution = 0.8, plot.SNN = T, 
#                      k.param = 10)
# pbmc <- FindClusters(pbmc, genes.use = setdiff(pbmc@var.genes, cellcycle), force.recalc = T, resolution = 0.8, plot.SNN = T, 
#                      k.param = 10)
TSNEPlot(pbmc, pt.size = 3, colors.use = pal_dolphin[-c(1:4)], do.label = F)
PCAPlot(pbmc, do.hover = F,  do.identify = F, 
        cols.use = pal_dolphin[-c(1:4)], pt.size =3, do.label = F)

pbmc <- SetAllIdent(pbmc, id = "type")
markers <- FindAllMarkers(pbmc, logfc.threshold = log(1.5), only.pos = T, test.use = "wilcox")
pbmc.markers <- subset(markers, p_val <= 0.001 & avg_logFC > 0)

library(dplyr)
library(pheatmap)
pbmc.markers %>% group_by(cluster) %>% top_n(20, -p_val) -> top10
DoHeatmap(pbmc, genes.use = top10$gene,#col.low = "darkblue", col.mid = "white",col.high = "red",
          use.scaled = T, slim.col.label = TRUE, remove.key = F, cex.col = 0.6)
# pheatmap(pbmc@data[c(top10$gene), order(pbmc@meta.data[,"orig.ident"])],
#          annotation_colors = ph_colors, fontsize_row = 5,
#          cluster_rows = F, cluster_cols = F, annotation_col = pbmc@meta.data[,c("orig.ident"), drop=F],
#          show_colnames = F, color = colorRampPalette(colors = c("darkblue","white","red"))(11), 
#          border_color = NA, gaps_col = cumsum(table(pbmc@meta.data$orig.ident)))
pheatmap(pbmc@data[c(top10$gene), order(pbmc@meta.data[,"type"])],
         annotation_colors = ph_colors, fontsize_row = 5,
         cluster_rows = F, cluster_cols = F, annotation_col = pbmc@meta.data[,c("type"), drop=F],
         show_colnames = F, color = gplots::colorpanel(n = 11, "darkblue","white","red"), 
         border_color = NA, gaps_col = cumsum(table(pbmc@meta.data$type)))
pheatmap(pbmc@data[c("CD34","CD38","THY1","PTPRC","ITGA6","MME","FLT3","SELL","IL3RA","TFRC","GYPA",
                     "PECAM1","CDH5","EPCAM","LEPR","LEPROT","PDGFRA","PDGFRB",
                     "CD19","CD33","ITGAM","ITGA2B","KIT","CD36","EPOR",
                     "PROCR","VNN2","CD9","GFI1B"), order(pbmc@meta.data[,"type"])],
         annotation_colors = ph_colors, fontsize_row = 10,
         cluster_rows = F, cluster_cols = F, annotation_col = pbmc@meta.data[,c("type"), drop=F],
         show_colnames = F, color = gplots::colorpanel(n = 11, "darkblue","white","red"), 
         border_color = NA, gaps_col = cumsum(table(pbmc@meta.data$type)))

pbmc <- AddCCscore(pbmc, add.dr = T)
pbmc <- AddAVscore(pbmc, gene.max = 10, add.dr = T)

myGenePlot(pbmc, gene1 = "FABP4", gene2 = "NPR3", cols.use = pal_dolphin[-c(2,3)], group.by = "Type")
VlnPlot(pbmc, features.plot = c("Fabp4","Npr3"), group.by = "Type",
        cols.use =  pal_dolphin[-c(2,3)], x.lab.rot = T, nCol = 2)
myGenePlot(pbmc, gene1 = "Fabp4", gene2 = "Npr3", cols.use = pal_dolphin[-c(1:4)], group.by = "Cluster")
VlnPlot(pbmc, features.plot = c("Fabp4","Npr3"), group.by = "Cluster",
        cols.use =  pal_dolphin[-c(1:4)], x.lab.rot = T, nCol = 2)
myFeaturePlot(pbmc, features.plot = c("Fabp4","Npr3"), dr = "avscore")
myFeaturePlot(pbmc, features.plot = c("Fabp4","Npr3"), dr = "ccscore")
myFeaturePlot(pbmc, features.plot = c("Fabp4","Npr3"), dr = "tsne")
myFeaturePlot(pbmc, features.plot = c("Fabp4","Npr3"), dr = "pca")

myGenePlot(pbmc, "G1S.Score", "G2M.Score", group.by = "Type", cols.use = pal_dolphin[-c(2,3)],
           xlim = c(0,4), ylim = c(0,6), main = "Cell Cycle Score")
segments(2,0, 2,2)
segments(0,2, 4,2)
segments(2,2, 4,6)
text(2,0.3,labels = c("Quiescent"), adj = 1.1)
text(4,0.3,labels = c("G1"), adj = 1.1)
text(4,4,labels = c("S"), adj = 1.1)
text(1,5,labels = c("G2/M"), adj = 1.1)

pbmc@meta.data$Cluster <- pbmc@meta.data$res.0.8

myCCPlot(pbmc, group.by = "type", th.g1s = 2, th.g2m = 2, cols.use = pal_dolphin[1:4])
myCCPlot(pbmc, group.by = "Cluster", th.g1s = 2, th.g2m = 2, cols.use = pal_dolphin[-c(1:4)])

myGenePlot(pbmc, "Artery.Score", "Venous.Score", group.by = "Type", cols.use = pal_dolphin[-c(2,3)],
           xlim = c(0,10), ylim = c(0,10), main = "Arteriovenous score")


select.genes <- subset(pbmc.markers, cluster %in% "VW")$gene
select.genes <- names(sort(pbmc@dr$pca@gene.loadings.full[,2], decreasing = F)[1:300])
write.table(select.genes, file = "clipboard", sep = ",",quote = F, row.names = F, col.names = F)

go <- myGOEnrich(pbmc.markers)



myGOEnrich <- function(pbmc.markers, clusters = NULL, width = 600, height = 350){
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(stringr)
  if(is.null(clusters)){
    clusters <- sort(unique(pbmc.markers$cluster))
  }
  for(i in clusters){
    gogenes <- subset(pbmc.markers, cluster == i)$gene
    
    engo <- enrichGO(gene         = gogenes,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "SYMBOL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
    png(filename = paste0("GO_BP results of ",i,".png"), width = width, height = height)
    print(dotplot(engo, title = paste0("GO_BP results of ",i))+ scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
    dev.off()
  }
}

pairwise.topgene.corr <- function(dataMat, ngene = 500, method = "spearman"){
  order_gene <- apply(dataMat, 2, function(x){order(x, decreasing = T)[1:ngene]})
  resMat <- matrix(NA, nrow = ncol(dataMat), ncol = ncol(dataMat), dimnames = list(colnames(dataMat), colnames(dataMat)))
  for(i in 1:ncol(dataMat)){
    for(j in i:ncol(dataMat)){
      if(i == j){
        resMat[i,j] <- 0
      }else{
        data.x <- dataMat[,i,drop = T]
        data.y <- dataMat[,j,drop = T]
        sel_gene <- unique(c(order_gene[,i], order_gene[,j]))
        resMat[i,j] <- (1-cor(data.x[sel_gene], data.y[sel_gene], method = method))/2
        resMat[j,i] <- resMat[i,j]
      }
    }
  }
  return(resMat)
}

pairwise.topgene.equal <- function(dataMat, ngene = 500){
  order_gene <- apply(dataMat, 2, function(x){order(x, decreasing = T)[1:ngene]})
  resMat <- matrix(NA, nrow = ncol(dataMat), ncol = ncol(dataMat), dimnames = list(colnames(dataMat), colnames(dataMat)))
  for(i in 1:ncol(dataMat)){
    for(j in i:ncol(dataMat)){
      if(i == j){
        resMat[i,j] <- 0
      }else{
        equal <- length(intersect(order_gene[,i], order_gene[,j]))
        resMat[i,j] <- 1 - equal/(2*ngene-equal)
        resMat[j,i] <- resMat[i,j]
      }
    }
  }
  return(resMat)
}


focusTSNEPlot <- function(object, group.by = 'orig.ident', focus = NULL, each = F, colors.use = NULL,...){
  cats <- sort(unique(object@meta.data[[group.by]]))
  focus = intersect(focus, cats)
  if(each){
    sapply(focus, function(x){
      TSNEPlot(object, group.by =group.by, ...)
    })
  } else{
    if(is.null(names(colors.use)))
      grey_ind <- which(!(cats %in% focus))
    else
      grey_ind <- which(!(names(colors.use) %in% focus))
    colors.use[grey_ind] <- "grey"
    top <- object@meta.data[[group.by]] %in% focus
    toprow <- c(rownames(object@meta.data)[which(!top)],rownames(object@meta.data)[which(top)])
    object@dr$tsne@cell.embeddings <- object@dr$tsne@cell.embeddings[toprow,]
    TSNEPlot(object, group.by = group.by, colors.use = colors.use, ...)
  }
}

AddCCscore <- function(pbmc, th.g1s = 2, th.g2m = 2, org = c("mouse","human"), gene.max = NULL, add.dr = F, use.scale = F){
  geneset1 <- toupper(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"))
  geneset2 <- toupper(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))
  org = org[1]
  if(org %in% c("human","hsa","hg19")){
    geneset1 <- toupper(geneset1)
    geneset2 <- toupper(geneset2)
  }
  use.data <- if(use.scale) pbmc@scale.data else pbmc@data
  use.data <- use.data[c(geneset1,geneset2),]
  use.data <- if(is.null(gene.max)) use.data else t(apply(use.data, 1, function(x) gene.max*x/max(x)))
  
  pbmc@meta.data$G1S.Score <- colMeans(use.data[geneset1, rownames(pbmc@meta.data)])
  pbmc@meta.data$G2M.Score <- colMeans(use.data[geneset2, rownames(pbmc@meta.data)])
  
  pbmc@meta.data$cc.phase[pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m] <- "Quiescent"
  pbmc@meta.data$cc.phase[!(pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m) & (pbmc@meta.data$G1S.Score < pbmc@meta.data$G2M.Score)] <- "G2/M"
  pbmc@meta.data$cc.phase[!(pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m) & (pbmc@meta.data$G1S.Score > pbmc@meta.data$G2M.Score) & pbmc@meta.data$G2M.Score < th.g2m] <- "G1"
  pbmc@meta.data$cc.phase[!(pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m) & (pbmc@meta.data$G1S.Score > pbmc@meta.data$G2M.Score) & pbmc@meta.data$G2M.Score > th.g2m] <- "S"
  pbmc@meta.data$cc.phase <- factor(pbmc@meta.data$cc.phase, levels = c("Quiescent","G1","S","G2/M"))
  
  
  if(add.dr){
    pbmc@dr$ccscore <- pbmc@dr$pca
    pbmc@dr$ccscore@cell.embeddings <- as.matrix(pbmc@meta.data[,c("G1S.Score","G2M.Score")])
  }
  pbmc
}

myGenePlot <- function(pbmc, gene1, gene2, use.scaled = F, group.by = NULL, cols.use = colors(),...){
  group.colors <- if(is.null(group.by)) pbmc@ident else pbmc@meta.data[[group.by]]
  plot(FetchData(pbmc, vars.all = gene1, use.scaled = use.scaled)[,1], FetchData(pbmc, vars.all = gene2,use.scaled = use.scaled)[,1],
       pch = 19, cex=1.6, col = cols.use[as.factor(group.colors)], xlab = gene1, ylab = gene2, ...)
}
myFeaturePlot <- function(pbmc, features.plot, nrow = NULL, ncol = NULL, dr = c("tsne","pca","ccscore","avscore"), cc.args = list(th.g1s = 2, th.g2m = 2),...){
  require(ggplot2)
  require(gridExtra)
  dr <- dr[1]
  ggData <- as.data.frame(cbind(pbmc@dr[[dr]]@cell.embeddings,FetchData(pbmc, features.plot)))
  colnames(ggData) <- c(colnames(pbmc@dr[[dr]]@cell.embeddings),gsub("-",".",features.plot))
  # print(feature.tmp)
  # ggData[,feature.tmp] <- t(pbmc@data[feature.tmp,])
  if(dr == 'tsne') {
    xx <- "tSNE_1"
    yy <- "tSNE_2"
  }
  if(dr == 'pca'){
    xx <- "PC1"
    yy <- "PC2"
  }
  if(dr == 'ccscore'){
    xx <- "G1S.Score"
    yy <- "G2M.Score"
  }
  if(dr == 'avscore'){
    xx <- "Artery.Score"
    yy <- "Venous.Score"
  }
  ggl <- lapply(features.plot, function(feature){
    p <- ggplot(ggData) + geom_point(mapping = aes_string(x = xx, y = yy, color = gsub("-",".",feature)), size = 2) + 
      scale_color_gradientn(colours = c("grey","yellow","red")) + 
      xlab(label = xx) + ylab(label = yy) +
      theme(legend.title = element_blank()) + ggtitle(feature) 
    if(dr == "ccscore"){
      ccx <- ceiling(max(pbmc@meta.data$G1S.Score))
      ccy <- ceiling(max(pbmc@meta.data$G2M.Score))
      p <- p + geom_linerange(mapping = aes_(x = th.g1s, ymin = 0, ymax = th.g1s)) + 
        geom_segment(mapping = aes_(x = 0, y = th.g2m, xend = ccx, yend = th.g2m)) + 
        geom_segment(mapping = aes_(x = th.g1s, y = th.g2m, xend = ccx, yend = ccy)) +
        geom_text(mapping = aes_(th.g1s, quote(0.3), label = quote("Quiescent"), hjust = 1.1)) +
        geom_text(mapping = aes_(ccx, quote(0.3), label = quote("G1"), hjust = 1.1)) +
        geom_text(mapping = aes_(ccx, ccy-2, label = quote("S"), hjust = 1.1)) +
        geom_text(mapping = aes_(th.g1s/2, ccy-1, label = quote("G2M"), hjust = 1.1))
    }
    p
  })
  grid.arrange(grobs = ggl, nrow= nrow,ncol = ncol)
}

myCCPlot <- function(pbmc, group.by = NULL, th.g1s = 2, th.g2m = 2, text.size = 3, cols.use = palette(),...){
  group.colors <- sort(unique(if(is.null(group.by)) pbmc@ident else pbmc@meta.data[[group.by]]))
  require(dplyr)
  require(ggplot2)
  ccData <- pbmc@meta.data 
  ccx <- ceiling(max(ccData$G1S.Score))
  ccy <- ceiling(max(ccData$G2M.Score))
  p1 <- ggplot(ccData) + geom_point(mapping = aes_string("G1S.Score", "G2M.Score", color = group.by), size = 3) + 
    xlab(label = "G1/S phase score") + ylab("G2/M phase score") + ggtitle("Cell Cycle Analysis") +
    scale_color_manual(values = cols.use[as.factor(group.colors)]) +
    #scale_x_continuous(expand = c(0, 0), limits = c(0, ccx)) + 
    #scale_y_continuous(expand = c(0, 0), limits = c(0, ccy)) +
    geom_linerange(mapping = aes_(x = th.g1s, ymin = 0, ymax = th.g1s)) + 
    geom_segment(mapping = aes_(x = 0, y = th.g2m, xend = ccx, yend = th.g2m)) + 
    geom_segment(mapping = aes_(x = th.g1s, y = th.g2m, xend = ccx, yend = ccy)) +
    geom_text(mapping = aes_(th.g1s, quote(0.3), label = quote("Quiescent"), hjust = 1.1)) +
    geom_text(mapping = aes_(ccx, quote(0.3), label = quote("G1"), hjust = 1.1)) +
    geom_text(mapping = aes_(ccx, ccy-2, label = quote("S"), hjust = 1.1)) +
    geom_text(mapping = aes_(th.g1s/2, ccy-1, label = quote("G2M"), hjust = 1.1)) +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  
  ccData <- pbmc@meta.data %>% count_(vars = c(group.by, "cc.phase")) %>% group_by_(group.by) %>% arrange(desc(cc.phase)) %>%
    mutate(pct = n/sum(n),
           ypos = cumsum(pct) - 0.5*pct)
  p2 <- ggplot(ccData) + 
    geom_bar(mapping = aes_string(x = group.by, y = "pct", fill = "cc.phase"), stat = "identity", width = 0.8) + 
    geom_text(mapping = aes_string(x = group.by, y= "ypos",label = quote(paste0(ccData$n,", ",sprintf("%1.1f", 100*ccData$pct),"%"))), size = text.size) +
    ggtitle(label = "Cell cycle Distribution of Cells") + theme(plot.title = element_text(hjust = 0.5))
  #  pheatmap::pheatmap(expr, cluster_cols = F, cluster_rows = F, annotation_col = annot, border_color = NA)
  require(gridExtra)
  grid.arrange(grobs = list(p1,p2), ncol = 2)
}
