##wordclouds for Figure 4, Ivan Reyna-Llorens

library(reshape)
library(ggplot2)
library(plyr)
library(reshape2)
library(wordcloud)
library(tm)

setwd("E:/Documents/DNAse-seq/Motifs")

#Read Ranking files for motif families

zmwl<- read.delim("Data_Figure_4_Zm_WL_known_Families.txt", head=T, sep="")
sbwl<- read.delim("Data_Figure_4_Sb_WL_known_Families.txt", head=T, sep="")
siwl<- read.delim("Data_Figure_4_Si_WL_known_Families.txt", head=T, sep="")
bdwl<- read.delim("Data_Figure_4_Bd_WL_known_Families.txt", head=T, sep="")
atwl<- read.delim("Data_Figure_4_At_WL_Families.txt", head=T, sep="")
sbwlde<- read.delim("Data_Figure_4_Sb_WL_DAP_DE_families.txt", head=T, sep="")
sbbsde<- read.delim("Data_Figure_4_Sb_BS_DAP_DE_families.txt", head=T, sep="")

head(zmwl)
head(sbwl)
head(siwl)
head(bdwl)
head(atwl)
head(sbwlde)
head(sbbsde)

# Worclouds
#maize
pdf("Zmwl_wc.pdf")
wordcloud(rep(zmwl$Var1, zmwl$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()
#sorghum
pdf("Sbwl_wc.pdf")
wordcloud(rep(sbwl$Var1, sbwl$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()
#setaria
pdf("Siwl_wc.pdf")
wordcloud(rep(siwl$Var1, siwl$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()

#brachy
pdf("Bdwl_wc.pdf")
wordcloud(rep(bdwl$Var1, bdwl$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()

#arabidopsis
pdf("Atwl_wc.pdf")
wordcloud(rep(atwl$Var1, atwl$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()
#sorghum DE WL
pdf("Sbwlde_wc.pdf")
wordcloud(rep(sbwlde$Var1, sbwlde$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()

#sorghum DE BS
pdf("Sbbsde_wc.pdf")
wordcloud(rep(sbbsde$Var1, sbbsde$Freq),rot.per = 0,
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()

