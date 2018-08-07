######## SOM analysis ###
library(scrat)
# Apply a variance cutoff 
df <- seu@data[which(apply(seu@data,1,var)>1), ] 
df <- df[which(rowSums(df > 0) > 2 & rowSums(df > 5) < ncol(df)), ] 

# Set up the environment 
env <- scrat.new(list(dataset.name = "HEC_cluster", dim.1stLvlSom = "auto", dim.2ndLvlSom = "auto", training.extension = 1, 
                      rotate.SOM.portraits = 0, flip.SOM.portraits = F, database.dataset = "auto", 
                      database.id.type = "auto", geneset.analysis = T, geneset.analysis.exact = F, 
                      spot.coresize.modules = 3, spot.threshold.modules = 0.95, spot.coresize.groupmap = 5, 
                      spot.threshold.groupmap = 0.75, 
                      pseudotime.estimation = list(n.waypoints = 20,  n.iterations = 20,  k = 30,  I = 5, initiator.sample = 1 ), 
                      feature.centralization = T,  sample.quantile.normalization = T,  pairwise.comparison.list = list() ) )                     


# Load input data into environment 
env$indata <- df                       
# Define sample groups 
env$group.labels <- seu@meta.data$cluster
#env$group.labels <- "auto"                       
# Define sample colors (optional) 
seu@meta.data$color <- ph_colors$cluster[match(seu@meta.data$cluster, names(ph_colors$cluster))] 
env$group.colors <-  seu@meta.data$color
#env$group.colors<- c("col1","col2",...)[match(env$group.labels, unique(env$group.labels))]                       
# Define colorschemes (optional) 
env$color.palette.portraits <- colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2")) 
# env$color.palette.portraits <- colorRampPalette(c("white", "grey", "grey8")) # 
env$color.palette.heatmaps <- colorRampPalette(c("darkblue","white","red"))                         

# Pipeline execution 
scrat.run(env)  
