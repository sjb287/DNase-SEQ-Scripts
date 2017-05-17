setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP")
#Code required to make Figure 4D and F
#Figure 4D
#Read in files
Bd<-read.table("Data_Figure_4_Bd_WL_known.txt",head=T)
Sb<-read.table("Data_Figure_4_Sb_WL_known.txt",head=T)
Si<-read.table("Data_Figure_4_Si_WL_known.txt",head=T)
Zm<-read.table("Data_Figure_4_Zm_WL_known.txt",head=T)
At<-read.table("Data_Figure_4_At_WL_known.txt",head=T)
#Select relevant columns
Bd<-Bd[,c("Var1","Rank_Bd")]
Sb<-Sb[,c("Var1","Rank_Sb")]
Si<-Si[,c("Var1","Rank_Si")]
Zm<-Zm[,c("Var1","Rank_Zm")]
At<-At[,c("Var1","Rank_At")]
#Merge data for comparisons
BdSb<-merge(Bd,Sb)
SbSi<-merge(Sb,Si)
SbZm<-merge(Sb,Zm)
SbAt<-merge(Sb,At)
#plot data
plot(y=BdSb$Rank_Bd,x=BdSb$Rank_Sb)
plot(y=SbSi$Rank_Si,x=SbSi$Rank_Sb)
plot(y=SbZm$Rank_Zm,x=SbZm$Rank_Sb)
plot(y=SbAt$Rank_At,x=SbAt$Rank_Sb)
#test correlation between species
WL_SbSi<-cor.test(SbSi$Rank_Sb,SbSi$Rank_Si,method="kendall")
WL_SbZm<-cor.test(SbZm$Rank_Sb,SbZm$Rank_Zm,method="kendall")
WL_SbBd<-cor.test(BdSb$Rank_Sb,BdSb$Rank_Bd,method="kendall")
WL_SbAt<-cor.test(SbAt$Rank_Sb,SbAt$Rank_At,method="kendall")

#350x400

#Figure 4F
#Read in files
Sb_WL<-read.table("Data_Figure_4F_Sb_WL.txt",head=T)
Sb_BS<-read.table("Data_Figure_4F_Sb_BS.txt",head=T)

#Select relevant columns
Sb_WL<-Sb_WL[,c("Var1","Rank_Sb")]
Sb_BS<-Sb_BS[,c("Var1","Rank_Sb")]
colnames(Sb_WL)<-c("Var1","Rank_WL")
colnames(Sb_BS)<-c("Var1","Rank_BS")
Sb_BS_WL<-merge(Sb_WL,Sb_BS)

#Plot data
plot(x=Sb_BS_WL$Rank_WL,y=Sb_BS_WL$Rank_BS)

#test correlation
WLBS_Sb<-cor.test(Sb_BS_WL$Rank_WL,Sb_BS_WL$Rank_BS,method="kendall")