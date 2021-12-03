library(vegan)  
library(ComplexHeatmap)

#load_data
my.data.taxa <- read.csv("omega_genus_20210408.csv", header=T, check.names = FALSE,row.names="sampleID" )

TSS.divide = function(x){x/sum(x)}

data.TSS <- t(apply(my.data.taxa[,18:145], 1, TSS.divide))

### Heatmap (#Figure 1)
data.dist.g.taxa <- vegdist(data.TSS, method = "euclidean")
col.clus.taxa <- hclust(data.dist.g.taxa, method = "ward.D2")
dend.taxa<-as.dendrogram(col.clus.taxa)
plot(dend.taxa)

metadata <- my.data.taxa[,1:2]
metadata <- data.frame(my.data.taxa[, c("stage", "studyarm")])

ha <- HeatmapAnnotation(df = metadata, col = list(stage =c("week 0" = "#d7191c", "week 12" = "#fdae61"), 
                                                  studyarm = c("Biotene" = "#abd9e9", "Listerine" = "#2c7bb6")))

col.heatmap <- colorRampPalette(c("gray95", "purple", "midnightblue"), space = "rgb")(99)

#Figure 1
Heatmap(t(data.TSS[,1:25]), name = "Relative abundance", col = col.heatmap,top_annotation = ha, 
        column_dend_height = unit(3, "cm"), column_dend_reorder = FALSE, show_row_names = TRUE, 
        show_column_names = FALSE, cluster_columns=dend.taxa, cluster_rows=FALSE, show_row_dend=FALSE)
