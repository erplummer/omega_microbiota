library(ALDEx2)
library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr)


#Alpha diversity
my.data.asv <- read.csv("omega_asv_20210408.csv", header=T, check.names = FALSE, row.names="sampleID")

shannon.diversity<-diversity(my.data.asv[,20:1105], index="shannon")

shannon.diversity<-data.frame(shannondiversity)

alpha.boxplot <- data.frame(c(my.data.asv[, c("stage", "studyarm")], shannondiversity))

pAlpha <- ggplot(alpha.boxplot, aes(x=stage, y=shannondiversity)) + geom_boxplot(aes(colour=stage), lwd=0.8)  + scale_colour_manual(name="Study stage", values=c("#d7191c", "#fdae61")) + facet_wrap(~studyarm)
pAlpha <- pAlpha + theme_bw()
pAlpha <- pAlpha + xlab(" ") + ylab("Shannon diversity index")
pAlpha <- pAlpha + theme(axis.text.y=element_text(hjust=1, size=12))
pAlpha <- pAlpha + theme(axis.title.y=element_text(size=14))
pAlpha <- pAlpha + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
pAlpha <- pAlpha + theme(strip.text = element_text(size=14))
pAlpha <- pAlpha + theme(legend.title=element_text(size=14), legend.text=element_text(size=14))
pAlpha


##Aldex2 analysis
my.data.taxa <- read.csv("omega_genus_20210408.csv", header=T, check.names = FALSE, row.names="sampleID")

dim(my.data.taxa)
head(my.data.taxa[,18:22])

omega.cols <- c("#e66101","#fdb863", "#5e3c99", "#b2abd2")


remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

#############################
#######   W0 vs W12   ########
#############################

###Biotene

my.data.Biotene<-my.data.taxa[(my.data.taxa$stage== "week 0"|my.data.taxa$stage== "week 12")& my.data.taxa$studyarm == "Biotene",]
my.data.Biotene<-my.data.Biotene[order(my.data.Biotene$stage, my.data.Biotene$StudyID),]

Biotene.sub <- t(my.data.Biotene[,18:145])

Biotene.sub.rare.removed <- remove_rare(table=Biotene.sub, cutoff_pro=0.1)

Y.Biotene  <- my.data.Biotene$stage

aldex.Biotene <- aldex.clr(Biotene.sub.rare.removed, Y.Biotene, mc.samples=128, denom="all", verbose=F)

aldex.Biotene.tt <- aldex.ttest(aldex.Biotene, paired.test=TRUE, verbose=FALSE)

aldex.Biotene.effect <- aldex.effect(aldex.Biotene, verbose=FALSE, CI=TRUE, include.sample.summary = TRUE)

aldex.Biotene.all <- data.frame(aldex.Biotene.tt, aldex.Biotene.effect)

par(mfrow=c(1,2))

aldex.plot(aldex.Biotene.all, type="MA", test="wilcox")
aldex.plot(aldex.Biotene.all, type="MW", test="wilcox")

write.csv(as.matrix(aldex.Biotene.all), "Biotene_week12_genus_20210505.csv")


###Listerine

my.data.Listerine<-my.data.taxa[(my.data.taxa$stage== "week 0"|my.data.taxa$stage== "week 12")& my.data.taxa$studyarm == "Listerine",]
my.data.Listerine<-my.data.Listerine[order(my.data.Listerine$stage, my.data.Listerine$StudyID),]

Listerine.sub <- t(my.data.Listerine[,18:145])

Listerine.sub.rare.removed <- remove_rare(table=Listerine.sub, cutoff_pro=0.1)

Y.Listerine  <- my.data.Listerine$stage

aldex.Listerine <- aldex.clr(Listerine.sub.rare.removed, Y.Listerine, mc.samples=128, denom="all", verbose=F)

aldex.Listerine.tt <- aldex.ttest(aldex.Listerine, paired.test=TRUE, verbose=FALSE)

aldex.Listerine.effect <- aldex.effect(aldex.Listerine, verbose=FALSE, CI=TRUE, include.sample.summary = TRUE)

aldex.Listerine.all <- data.frame(aldex.Listerine.tt, aldex.Listerine.effect)

par(mfrow=c(1,2))

aldex.plot(aldex.Listerine.all, type="MA", test="wilcox")
aldex.plot(aldex.Listerine.all, type="MW", test="wilcox")

write.csv(as.matrix(aldex.Listerine.all), "Listerine_week12_genus_20210505.csv")

#############################
#####   Smoking Status  #####
#############################

my.data.smoker<-my.data.taxa[(my.data.taxa$stage== "week 0")&(my.data.taxa$Currentsmoker== "no"|my.data.taxa$Currentsmoker== "yes"),]
my.data.smoker<-my.data.smoker[order(my.data.smoker$Currentsmoker, my.data.smoker$StudyID),]

smoke.sub <- t(my.data.smoker[,18:145])

smoke.sub.rare.removed <- remove_rare(table=smoke.sub, cutoff_pro=0.1)

Y.Smoker  <- my.data.smoker$Currentsmoker

aldex.smoke <- aldex.clr(smoke.sub.rare.removed, Y.Smoker, mc.samples=128, denom="all", verbose=F)

aldex.smoke.tt <- aldex.ttest(aldex.smoke, paired.test=FALSE, verbose=FALSE)

aldex.smoke.effect <- aldex.effect(aldex.smoke, verbose=FALSE, CI=TRUE, include.sample.summary = TRUE)

aldex.smoke.all <- data.frame(aldex.smoke.tt, aldex.smoke.effect)

par(mfrow=c(1,2))

aldex.plot(aldex.smoke.all, type="MA", test="wilcox")
aldex.plot(aldex.smoke.all, type="MW", test="wilcox")

write.csv(as.matrix(aldex.smoke.all), "Smoking_genus_20210505.csv")

#############################
### NG at baseline  ###
#############################

my.data.NG <-my.data.taxa[(my.data.taxa$stage== "week 0")&(my.data.taxa$ThroatNGw0== "No"|my.data.taxa$ThroatNGw0== "Yes"),]
my.data.NG<-my.data.NG[order(my.data.NG$ThroatNGw0, my.data.NG$StudyID),]

NG.sub <- t(my.data.NG[,18:145])

NG.sub.rare.removed <- remove_rare(table=NG.sub, cutoff_pro=0.1)

Y.NG  <- my.data.NG$ThroatNGw0

aldex.NG  <- aldex.clr(NG.sub.rare.removed, Y.NG, mc.samples=128, denom="all", verbose=F)

aldex.NG.tt <- aldex.ttest(aldex.NG, paired.test=FALSE, verbose=FALSE)

aldex.NG.effect <- aldex.effect(aldex.NG, verbose=FALSE, CI=TRUE, include.sample.summary = TRUE)

aldex.NG.all <- data.frame(aldex.NG.tt, aldex.NG.effect)

par(mfrow=c(1,2))

aldex.plot(aldex.NG.all, type="MA", test="wilcox")
aldex.plot(aldex.NG.all, type="MW", test="wilcox")

write.csv(as.matrix(aldex.NG.all), "ThroatNG_genus_20210505.csv")


####Generate boxplots ####

#Read in data that has been re-formatted for drawing boxplots
Biotene.boxplot <- read.csv("Biotene boxplots 20210505.csv", header=T, check.names = FALSE,row.names="SampleID" )
Listerine.boxplot <- read.csv("Listerine boxplots 20210505.csv", header=T, check.names = FALSE,row.names="SampleID" )
NG.boxplot <- read.csv("Oral NG boxplots 20210505.csv", header=T, check.names = FALSE,row.names="SampleID" )
smk.boxplot <- read.csv("Smoking boxplots 20210505.csv", header=T, check.names = FALSE,row.names="SampleID" )


##Biotene and Listerine boxplots
Biotene.boxplot.m <- melt(Biotene.boxplot, id.var = "Stage")
Listerine.boxplot.m <- melt(Listerine.boxplot, id.var = "Stage")

pBiotene <- ggplot(data=Biotene.boxplot.m, aes(x=variable, y=value))  +   geom_boxplot(aes(colour=Stage), lwd=0.8)  + scale_colour_manual(name="Study stage", values=c("#d7191c", "#fdae61"))
pBiotene <- pBiotene + xlab(" ") + ylab("CLR transformed abundance") + ggtitle("Biot\u00e9ne")
pBiotene <- pBiotene + theme_bw()
pBiotene <- pBiotene + theme(plot.title = element_text(hjust = 0.5 , size = 16))
pBiotene <- pBiotene + theme(axis.text.x=element_text(angle=45,hjust=1, size=14), axis.title =element_text(size =14))
pBiotene <- pBiotene + theme(axis.text.y=element_text(hjust=1, size=12))
pBiotene <- pBiotene + annotate("text", x = 1, y = 12.5, hjust=0, size=6, label = "*")
pBiotene <- pBiotene + annotate("text", x = 2, y = 12.5, hjust=0, size=6, label = "*")
pBiotene <- pBiotene + theme(legend.position = "none")
pBiotene

pListerine <- ggplot(data=Listerine.boxplot.m, aes(x=variable, y=value)) + geom_boxplot(aes(colour=Stage), lwd=0.8)  + scale_colour_manual(name="Study stage", values=c("#d7191c", "#fdae61"))
pListerine <- pListerine + xlab(" ") + ylab("CLR transformed abundance") + ggtitle("Listerine")
pListerine <- pListerine + theme_bw()
pListerine <- pListerine + theme(plot.title = element_text(hjust = 0.5, size = 16))
pListerine <- pListerine + theme(axis.text.x=element_text(angle=45,hjust=1, size=14), axis.title =element_text(size =14))
pListerine <- pListerine + theme(axis.text.y=element_text(hjust=1, size=12))
pListerine <- pListerine + theme(legend.position = "none")
pListerine


##Generate a panel with alpha diversity plot and the aldex2 plots (Figure 3)
ggarrange(pAlpha, ggarrange(pBiotene, pListerine, nrow = 2, labels = c("B", "C")), widths =  c(1.0, 2.0), ncol=2, labels = c("A"), common.legend = TRUE, legend = "bottom") 


########### Oral NG and smoking box plots ######### 
NG.boxplot.m <- melt(NG.boxplot, id.var = "OralNG")
smk.boxplot.m <- melt(smk.boxplot, id.var = "Smoker")


## Supplementary Figure 4b
pNG <- ggplot(data=NG.boxplot.m, aes(x=variable, y=value))  + geom_boxplot(aes(colour=OralNG), lwd = 0.8)  + scale_colour_manual(name="Oropharyngeal gonorrhoea \ndetected by NAAT", values=c("#d01c8b", "#4dac26"))
pNG <- pNG + xlab(" ") + ylab("CLR transformed abundance")
pNG <- pNG + theme_bw()
pNG <- pNG + theme(plot.title = element_text(hjust = 0.5, size = 20))
pNG <- pNG + theme(axis.text.x=element_text(angle=45,hjust=1, size=20), axis.title =element_text(size =20))
pNG <- pNG + theme(axis.text.y=element_text(hjust=1, size=20))
pNG <- pNG + theme(legend.position = "none")
pNG


## Supplementary Figure 6b
psmk <- ggplot(data=smk.boxplot.m, aes(x=variable, y=value)) + geom_boxplot(aes(colour=Smoker), lwd = 0.8)  + scale_colour_manual(name="Current smoker", values=c("#a6611a", "#018571"))
psmk <- psmk + xlab(" ") + ylab("CLR transformed abundance")
psmk <- psmk + theme_bw()
psmk <- psmk + theme(plot.title = element_text(hjust = 0.5, size = 20))
psmk <- psmk + theme(axis.text.x=element_text(angle=45,hjust=1, size=20), axis.title =element_text(size =20))
psmk <- psmk + theme(axis.text.y=element_text(hjust=1, size=20))
psmk <- psmk + theme(legend.position = "none")
psmk
