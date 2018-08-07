##### MAGIC ###
library(Rmagic)
library(ggplot2)
library(readr)
library(viridis)
library(phateR)
data_EC <- t(data_umi)
keep_cols <- colSums(data_EC > 0) > 10
data_EC <- data_EC[, keep_cols]

ggplot() + geom_histogram(aes(x = rowSums(data_EC)), bins = 50) +
  geom_vline(xintercept = 100, color='red')

# keep_rows <- rowSums(data_EC) < 6e5 & rowSums(data_EC) > 1e4
# data_EC <- data_EC[keep_rows,]

data_EC <- library.size.normalize(data_EC)
data_EC <- sqrt(data_EC)

data_magic <- Rmagic::magic(data_EC, genes = "all_genes", 
                            seed = 123, k = 10, alpha = 15, t = 4)

ggplot(as.data.frame(cbind(data_EC, data_annot[rownames(data_EC),]))) +
  geom_point(aes(x = RUNX1, y = CD44, colour = stage), size = 2) +
  scale_colour_manual(values = pal_dolphin)



# data_magic <- magic(data, genes = c("RUNX1", "SOX17", "CD44"), seed = 123, t = 4)
ggplot(data_magic) +
  geom_point(aes(x = RUNX1, y = CD44, colour = CDH5), size = 1) +
  scale_colour_viridis(option="B")
ggplot(as.data.frame(cbind(data_magic$result, data_annot[rownames(data_EC),]))) +
  geom_point(aes(x = RUNX1, y = CD44, colour = stage), size = 2) +
  scale_colour_manual(values = pal_dolphin)


data_magic_pca <- magic(data_EC, genes="pca_only", t = 4, init = data_magic, seed = 1234)


data_pca <- cbind(data_magic_pca$result, data_annot[rownames(data_EC),])
data_pca <- cbind(data_pca, data_magic$result)
ggplot(data_pca) +
  geom_point(aes(x = PC1, y = PC2, color = stage), size = 2) +
  scale_color_manual(values = pal_dolphin) 
