######## monocle ########
setwd("D:/data/bzj2018human 小胚数据/monocle/")
options(stringsAsFactors = F)
rm(list = ls())

load("D:/data/bzj2018human 小胚数据/filtered_data_Matrix.RData")
load("D:/data/human HSC analysis//human HSC/bzj20180427/cellcyclegenes.RData")
load("D:/data/human HSC analysis//human HSC/bzj20180427/labcolors.RData")

#####
library(monocle)
library(scater, quietly = T)
library(knitr)

####### gene annotation from ensembl gene names ######
library(org.Hs.eg.db)
eg2symbol=toTable(org.Hs.egSYMBOL)
eg2ensembl=toTable(org.Hs.egENSEMBL)
# egid=eg2ensembl[ match(rownames(molecules),eg2ensembl$gene_id),'gene_id']
# symbol=eg2symbol[match( egid ,eg2symbol$gene_id),'symbol']
# gene_annotation = data.frame(ensembl=rownames(molecules),
#                              gene_short_name=symbol,
#                              egid=egid)

data_annot_meta$clusters <- data_annot_meta$res.1.2
pd <- new("AnnotatedDataFrame", data = data_annot_meta)
HSMM <- newCellDataSet(as.matrix(data_filtered), phenoData = pd)
#RNA count
# RC_matrix <- relative2abs(HSMM) 
# HSMM <- newCellDataSet(RC_matrix, phenoData = pd, expressionFamily = VGAM::negbinomial.size() )
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

####### Cell QC ########
# pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
# HSMM <- HSMM[, pData(HSMM)$Total_mRNAs < 1e6]
# upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
# lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
# table(pData(HSMM)$stage)

####### unsupervised cluster ########
disp_table <- dispersionTable(HSMM)
unsup_clustering_gene <- subset(disp_table, mean_expression >= 0.1)
 #### 14201 genes were used to cluster 
HSMM <- setOrderingFilter(HSMM, unsup_clustering_gene$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = "log"

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6, reduction_method = "tSNE", verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 8)
p1 <- plot_cell_clusters(HSMM, 1, 2, color_by = "stage", cell_size = 1, do.return = T)
p2 <- plot_cell_clusters(HSMM, 1, 2, color_by = "location", cell_size = 1 , do.return =T)
gridExtra::grid.arrange(p1, p2, nrow = 1)
plot_cell_clusters(HSMM, 1, 2, color_by = "Cluster", cell_size = 1)
plot_cell_clusters(HSMM, 1, 2, color_by = "clusters", cell_size = 1)
plot_cell_clusters(HSMM, 1, 2, color_by = "State", cell_size = 1)

##### semisupvised cluster ####
### DEG ###
markers <- differentialGeneTest(HSMM[unsup_clustering_gene$gene_id, ], fullModelFormulaStr = "~Cluster",cores = 1)

####### Pseudotime analysis #######
### 1. ordering based on genes that differ between clusters ###
if (F) {
  differ_test_res <- differentialGeneTest(HSMM[unsup_clustering_gene$gene_id, ], 
                                           fullModelFormulaStr = "~Cluster",cores = 1)
  ordering_gene <- row.names(subset(differ_test_gene, qval < 0.01))
}

### 2. Selecting genes with high dispersion across cells ###
disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table, 
                         # mean_expression >= 0.5 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id
ordering_genes <- setdiff(ordering_genes, cellcycle)
HSMM<- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM, reverse = FALSE)

t1 <- plot_cell_trajectory(HSMM, color_by = "stage", cell_size = 2, do.return = T)
t2 <- plot_cell_trajectory(HSMM, color_by = "location", cell_size = 2, do.return = T)
t3 <- plot_cell_trajectory(HSMM, color_by = "clusters", cell_size = 2, do.return = T)
t4 <- plot_cell_trajectory(HSMM, color_by = "State", cell_size = 2, do.return = T)
gridExtra::grid.arrange(t1, t2, t3, t4, nrow =2)
plot_cell_trajectory(HSMM, color_by = "Pseudotime", cell_size = 2, do.return = T)

Hemato_markers <- c("CD34","CD38","THY1","PTPRC","ITGA6","MME","FLT3","SELL","IL3RA","TFRC","GYPA", "PECAM1",
                    "CDH5","EPCAM","LEPR","LEPROT","PDGFRA","PDGFRB", "CD19","CD33","ITGAM","ITGA2B",
                    "KIT","CD36","EPOR", "PROCR","VNN2","CD9","GFI1B")
EC_markers <- c("RUNX1", "GATA1", "MYB", "GATA3", "BMP4", "LMO2", "ITGA2B", "THY1", "KIT", "PTPRC", "CD38",
                "CDH5", "PECAM1", "TEK", "TAL1", "FLI1")
plot_genes_in_pseudotime(HSMM[EC_markers[1:8], ], color_by = "clusters", ncol = 2, cell_size = 1)
plot_genes_in_pseudotime(HSMM[EC_markers[9:16], ], color_by = "clusters", ncol = 2, cell_size = 1)
plot_genes_in_pseudotime(HSMM[Hemato_markers[1:10], ], color_by = "clusters", ncol = 2, cell_size = 1)
plot_genes_in_pseudotime(HSMM[Hemato_markers[11:20], ], color_by = "clusters", ncol = 2, cell_size = 1)
plot_genes_in_pseudotime(HSMM[Hemato_markers[21:29], ], color_by = "clusters", ncol = 2, cell_size = 1)








library(pheatmap)
ph_colors <- list(stage = setNames(myrainbow[1:length(unique(HSMM@phenoData$stage))], sort(unique(HSMM@phenoData$stage))), 
                  location = setNames(pal_dolphin[1:length(unique(HSMM@phenoData$location))], sort(unique(HSMM@phenoData$location))),
                  clusters = setNames(bertie.color[1:length(unique(HSMM@phenoData$clusters))], sort(unique(HSMM@phenoData$clusters))),
                  State = setNames(mycolors[1:length(unique(HSMM@phenoData$State))], sort(unique(HSMM@phenoData$State))))
heatmap_Data <- apply(data_filtered, 2, function(x) {log2(1e5*x/sum(x) + 1)})
data_annot_heatmap <- HSMM@phenoData@data[, c("clusters", "State")]
pheatmap(heatmap_Data[ordering_genes, order(data_annot_heatmap[,"clusters"])],
         annotation_colors = ph_colors, fontsize_row = 7, show_rownames = F,
         cluster_rows = T, cluster_cols = F, annotation_col = data_annot_heatmap[,c("clusters", "State")],
         show_colnames = F, color =colorRampPalette(colors = c("#FF00FF","#000000","#FFFF00"))(100), 
         border_color = NA)




