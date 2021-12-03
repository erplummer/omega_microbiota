library(ggplot2)
library(mixOmics)
library(vegan)

#load_data
my.data.taxa <- read.csv("omega_asv_20210408.csv", header=T, check.names = FALSE, row.names="sampleID")

############################################################
##Function to remove ASVs with 0 read
############################################################
remove.zero.cols <- function(x){
  if(is.numeric(x)){
    sum(x) > 0
  } else {
    TRUE
  }
}

#####################################################################################
##Low count removal - from http://mixomics.org/mixmc/mixmc-pre-processing/
#####################################################################################
low.count.removal = function(
  data, 
  percent=0.01 
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

###### Colour selection
omega.cols <- c("#e66101","#fdb863", "#5e3c99", "#b2abd2")


######################################################################################
##Pre-processing of data - from http://mixomics.org/mixmc/mixmc-pre-processing/
######################################################################################
dataOffset <- my.data.taxa[,20:1105] + 1

result.filter.main <- low.count.removal(dataOffset, percent=0.01)
data.filter.main <- result.filter.main$data.filter

#####################################################################
##PCA (Principal Component Analysis) - All samples included (Fig 2a)
#####################################################################
pca.result <- pca(as.matrix(data.filter.main), logratio = 'CLR', ncomp=10)

plot(pca.result, ncomp = 10, explained.var = TRUE)

plotIndiv(pca.result, group = my.data.taxa$Group4, comp=c(1,2), legend = TRUE, 
          legend.title="Sample category", style="ggplot2", ellipse = TRUE, 
          ellipse.level = 0.95,  centroid = FALSE, title = "All samples", 
          ind.names = F, col = omega.cols, size.xlabel = rel(1.5), size.ylabel = rel(1.5), 
          size.legend.title = rel(2.5), size.legend = rel(2.5), point.lwd = 1.5)

#############################################################
##Looking at W0 vs W12 samples in Listerine - (Fig 2b)
#############################################################

my.data.Listerine<-my.data.taxa[my.data.taxa$studyarm== "Listerine",]
my.data.Listerine<-my.data.Listerine[, sapply(my.data.Listerine,  remove.zero.cols)]

dataOffset.L <- my.data.Listerine[,20:966] + 1

result.filter.L <- low.count.removal(dataOffset.L, percent=0.01)
data.filter.L <- result.filter.L$data.filter

pca.result.Listerine <- pca(as.matrix(data.filter.L), logratio = 'CLR', ncomp=10)

plot(pca.result.Listerine, ncomp = 10, explained.var = TRUE)

pListerine<-plotIndiv(pca.result.Listerine, group = my.data.Listerine$stage, comp=c(1,2), 
          legend = FALSE, style="ggplot2", ellipse = TRUE, 
          ellipse.level = 0.95,  centroid = FALSE, title = 'Listerine', 
          ind.names = F,col = omega.cols[3:4], pch=c(3,4), size.xlabel = rel(1.5), 
          size.ylabel = rel(1.5), point.lwd = 1.5, )


#############################################################
##Looking at W0 vs W12 samples in Biotene - (Fig 2c)
#############################################################

my.data.Biotene<-my.data.taxa[my.data.taxa$studyarm== "Biotene",]
my.data.Biotene<-my.data.Biotene[, sapply(my.data.Biotene,  remove.zero.cols)]

dataOffset.B <- my.data.Biotene[,20:982] + 1

result.filter.B <- low.count.removal(dataOffset.B, percent=0.01)
data.filter.B <- result.filter.B$data.filter

pca.result.Biotene<- pca(as.matrix(data.filter.B), logratio = 'CLR', ncomp=10)

plot(pca.result.Biotene, ncomp = 10, explained.var = TRUE)

pBiotene<-plotIndiv(pca.result.Biotene, group = my.data.Biotene$stage, comp=c(1,2), 
          legend = FALSE, style="ggplot2", ellipse = TRUE, 
          ellipse.level = 0.95,  centroid = FALSE, title = 'Biot\u00e8ne', 
          ind.names = F,col = omega.cols[1:2], pch=c(1,2), size.xlabel = rel(1.5), 
          size.ylabel = rel(1.5), point.lwd = 1.5)


########################################################################################
##Looking at W0 samples - any difference between NG infection vs no infection? (Supp Fig 4a)
########################################################################################

my.data.NG<-my.data.taxa[(my.data.taxa$stage== "W0")&(my.data.taxa$ThroatNGw0== "No"|my.data.taxa$ThroatNGw0== "Yes"),]

my.data.NG<-my.data.NG[, sapply(my.data.NG,  remove.zero.cols)]

dataOffset.NG <- my.data.NG[,20:1018] + 1

result.filter.NG <- low.count.removal(dataOffset.NG, percent=0.01)
data.filter.NG <- result.filter.NG$data.filter

col_NG=c("#d01c8b", "#4dac26")

pca.result.NG <- pca(as.matrix(data.filter.NG), logratio = 'CLR', ncomp=10)

plot(pca.result.NG, ncomp = 10, explained.var = TRUE)

plotIndiv(pca.result.NG, group = my.data.NG$ThroatNGw0, comp=c(1,2), 
          legend = TRUE, legend.title="N.gonorrhoea detected by NAAT", 
          style="ggplot2", ellipse = TRUE, ellipse.level = 0.95, 
          title = 'Oropharyngeal gonorrhoea detected by NAAT', 
          ind.names = F, col = col_NG, 
          size.xlabel = rel(1.8), size.ylabel = rel(1.8),  
          size.legend.title = rel(1.8), size.legend = rel(1.8), point.lwd = 1.5, legend.position = "bottom")

#################################################################################
##Looking at D0 samples - any difference between smokers vs non smokers? (Supp Fig 6a)
#################################################################################

my.data.smoker<-my.data.taxa[(my.data.taxa$stage== "W0")&(my.data.taxa$Currentsmoker== "no"|my.data.taxa$Currentsmoker== "yes"),]

my.data.smoker<-my.data.smoker[, sapply(my.data.smoker,  remove.zero.cols)]
dim(my.data.smoker)

dataOffset.smk <- my.data.smoker[,20:1007] + 1

result.filter.smk <- low.count.removal(dataOffset.smk, percent=0.01)
data.filter.smk <- result.filter.smk$data.filter

pca.result.smk <- pca(as.matrix(data.filter.smk), logratio = 'CLR', ncomp=10)

plot(pca.result.smk, ncomp = 10, explained.var = TRUE)

col_smk=c("#a6611a", "#018571")

plotIndiv(pca.result.smk, group = my.data.smoker$Currentsmoker, comp=c(1,2), 
          legend = TRUE, legend.title="Current smoker", style="ggplot2", 
          ellipse = TRUE, title = 'Smoking', ind.names = F, col = col_smk, 
          size.xlabel = rel(1.8), size.ylabel = rel(1.8),  
          size.legend.title = rel(1.8), size.legend = rel(1.8), point.lwd = 1.5, legend.position = "bottom")


######################## Adonis analysis ########################
TSS.divide = function(x){x/sum(x)}

##Listerine week0 vs week12
my.data.List.TSS<- t(apply(as.matrix(data.filter.L), 1, TSS.divide))
adon_ng<-adonis(my.data.List.TSS ~ my.data.Listerine$stage, method="bray",perm=999)
print(adon_ng)

##Biotene week0 vs week12
my.data.Biot.TSS<- t(apply(as.matrix(data.filter.B), 1, TSS.divide))
adon_ng<-adonis(my.data.Biot.TSS ~ my.data.Biotene$stage, method="bray",perm=999)
print(adon_ng)

##NG at baseline
my.data.NG.TSS<- t(apply(as.matrix(data.filter.NG), 1, TSS.divide))
adon_ng<-adonis(my.data.NG.TSS ~ my.data.NG$ThroatNGw0, method="bray",perm=999)
print(adon_ng)

##Smoking
my.data.smoker.TSS<- t(apply(as.matrix(data.filter.smk), 1, TSS.divide))
adon_smk<-adonis(my.data.smoker.TSS ~ my.data.smoker$Currentsmoker, method="bray",perm=999)
print(adon_smk)