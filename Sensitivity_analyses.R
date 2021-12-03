library(ggplot2)
library(mixOmics)
library(vegan)

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

############################################################
##Function for total sum scale (TSS)
############################################################
TSS.divide = function(x){x/sum(x)}

############################################################
##Low count removal - from http://mixomics.org/mixmc/mixmc-pre-processing/
############################################################
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}


#############################################################
## PERMANOVA W0 vs W12 samples in Listerine
#############################################################

my.data.Listerine<-my.data.taxa[my.data.taxa$studyarm== "Listerine"&my.data.taxa$dailyMouthwash== "No",]

my.data.Listerine<-my.data.Listerine[, sapply(my.data.Listerine,  remove.zero.cols)]
head(my.data.Listerine[,20:916])

dataOffset.L <- my.data.Listerine[,20:916] + 1

result.filter.L <- low.count.removal(dataOffset.L, percent=0.01)
data.filter.L <- result.filter.L$data.filter

my.data.Listerine.TSS <- t(apply(data.filter.L, 1, TSS.divide))

adon_List<-adonis(my.data.Listerine.TSS ~ my.data.Listerine$stage, method="bray",perm=999)
print(adon_List)

#############################################################
##PERMANOVA W0 vs W12 samples in Biotene
#############################################################

my.data.Biotene<-my.data.taxa[my.data.taxa$studyarm== "Biotene"&my.data.taxa$dailyMouthwash== "No",]

my.data.Biotene<-my.data.Biotene[, sapply(my.data.Biotene,  remove.zero.cols)]
dim(my.data.Biotene)
head(my.data.Biotene[,20:916])

dataOffset.B <- my.data.Biotene[,20:916] + 1

result.filter.B <- low.count.removal(dataOffset.B, percent=0.01)
data.filter.B <- result.filter.B$data.filter
my.data.Biotene.TSS <- t(apply(data.filter.B, 1, TSS.divide))

adon_Biot<-adonis(my.data.Biotene.TSS ~ my.data.Biotene$stage, method="bray",perm=999)
print(adon_Biot)


#############################################################
##Looking at mouthwash use at baseline
#############################################################

my.data.moutwashUse<-my.data.taxa[(my.data.taxa$stage== "W0")&(my.data.taxa$dailyMouthwash== "No"|my.data.taxa$dailyMouthwash== "Yes"),]
dim(my.data.moutwashUse)
head(my.data.moutwashUse[,20:1105])

my.data.moutwashUse<-my.data.moutwashUse[, sapply(my.data.moutwashUse,  remove.zero.cols)]
dim(my.data.moutwashUse)

dataOffset.mwashuse <- my.data.moutwashUse[,20:1014] + 1

result.filter.mwashuse <- low.count.removal(dataOffset.mwashuse, percent=0.01)
data.filter.mwashuse <- result.filter.mwashuse$data.filter

my.data.mwashuse.TSS <- t(apply(data.filter.mwashuse, 1, TSS.divide))

pca.result.mwashUse <- pca(my.data.mwashuse.TSS, logratio = 'CLR', ncomp=10)

plot(pca.result.mwashUse, ncomp = 10, explained.var = TRUE)

##Supp Figure 3A
plotIndiv(pca.result.mwashUse, group = my.data.moutwashUse$dailyMouthwash, comp=c(1,2), 
          legend = TRUE, legend.title="Daily mouthwash use", style="ggplot2", ellipse = TRUE, 
          ellipse.level = 0.95,  centroid = FALSE, title = 'Daily  mouthwash use at week 0', 
          ind.names = F, col = c("#2c7bb6","#d7191c"), pch=c(2,1), size.xlabel = rel(1.8), 
          size.ylabel = rel(1.8),  size.legend.title = rel(1.8), size.legend = rel(1.8), 
          point.lwd = 1.5, legend.position = "bottom")

adon_Mouthwash_bin<-adonis(my.data.mwashuse.TSS ~ my.data.moutwashUse$dailyMouthwash, method="bray",perm=999)
print(adon_Mouthwash_bin)

##Supp Figure 3B
plotIndiv(pca.result.mwashUse, group = my.data.moutwashUse$mouthwash_freq, comp=c(1,2),
          legend = TRUE, legend.title="Mouthwash use frequency", style="ggplot2", ellipse = TRUE,
          ellipse.level = 0.95,  centroid = FALSE, title = 'Mouthwash use at week 0',
          ind.names = F, col = c("#d7191c","#abd9e9","#2c7bb6","#fdae61"), size.xlabel = rel(1.8), 
          size.ylabel = rel(1.8),  size.legend.title = rel(1.8), size.legend = rel(1.8), 
          point.lwd = 1.5, legend.position = "bottom")

adon_Mouthwash_4grp<-adonis(my.data.mwashuse.TSS ~ my.data.moutwashUse$mouthwash_freq, method="bray",perm=999)
print(adon_Mouthwash_4grp)


#########################################################################################
##Looking at W0 vs W12 samples in Listerine among men reporting no mouthwash use at baseline
#########################################################################################

my.data.Listerine.sens2<-my.data.taxa[my.data.taxa$studyarm== "Listerine"&my.data.taxa$mouthwash_freq== "never",]

my.data.Listerine.sens2<-my.data.Listerine.sens2[, sapply(my.data.Listerine.sens2,  remove.zero.cols)]
dim(my.data.Listerine.sens2)

dataOffset.L.sens2 <- my.data.Listerine.sens2[,20:740] + 1

result.filter.L.sens2 <- low.count.removal(dataOffset.L.sens2, percent=0.01)
data.filter.L.sens2 <- result.filter.L.sens2$data.filter

my.data.Listerine.sens2.TSS <- t(apply(data.filter.L.sens2, 1, TSS.divide))

adon.List.sens2<-adonis(my.data.Listerine.sens2.TSS ~ my.data.Listerine.sens2$stage, method="bray",perm=999)
print(adon.List.sens2)


######################################################################################
##Looking at W0 vs W12 samples in Biotene among men reporting no mouthwash use at baseline
######################################################################################

my.data.Biotene.sens2<-my.data.taxa[my.data.taxa$studyarm== "Biotene"&my.data.taxa$mouthwash_freq== "never",]

my.data.Biotene.sens2<-my.data.Biotene.sens2[, sapply(my.data.Biotene.sens2,  remove.zero.cols)]
dim(my.data.Biotene.sens2)

dataOffset.B.sens2 <- my.data.Biotene.sens2[,20:671] + 1

result.filter.B.sens2 <- low.count.removal(dataOffset.B.sens2, percent=0.01)
data.filter.B.sens2 <- result.filter.B.sens2$data.filter

my.data.Biotene.sens2.TSS <- t(apply(data.filter.B.sens2, 1, TSS.divide))

adon.Biot.sens2<-adonis(my.data.Biotene.sens2.TSS ~ my.data.Biotene.sens2$stage, method="bray",perm=999)
print(adon.Biot.sens2)