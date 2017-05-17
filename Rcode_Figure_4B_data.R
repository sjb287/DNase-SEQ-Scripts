setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/")

#Read in DGFs annotated with motifs Sb_WL
Sb_WL_IRL<-read.table("Sb_WL_DGF_fdr_0.01_fimo_hits_all_precise.bed")
Sb_WL_IRL<-Sb_WL_IRL[!(duplicated(Sb_WL_IRL)),]
#Subset DGFs to only keep those annotated with motif 372
x<-(subset(Sb_WL_IRL,Sb_WL_IRL$V4=="motif_/s_372"))
x<-x[,c(1:3)]
write.table(x,"Data_Figure_4B.bed",quote=F,row.name=F,col.names=F,sep="\t")