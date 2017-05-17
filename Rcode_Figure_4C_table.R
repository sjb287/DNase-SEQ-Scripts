#R code to make data for Figure 4C
#################################################
setwd("/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/")
BED_DGF_Si_WL<-read.table("DGF_Si_WL.bed")
BED_DGF_Si_BS<-read.table("DGF_Si_BS.bed")
BED_DGF_Sb_WL<-read.table("DGF_Sb_WL.bed")
BED_DGF_Sb_BS<-read.table("DGF_Sb_BS.bed")
BED_DGF_Zm_BS<-read.table("DGF_Zm_BS.bed")
BED_DGF_Zm_WL<-read.table("DGF_Zm_WL.bed")
BED_DGF_Bd_WL<-read.table("DGF_Bd_WL.bed")

setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP")
#Read in all TF binding sites known and de novo
IRL_DAP_Sb_WL<-read.table("Sb_WL_DGF_annotated_known_and_de_novo.bed")
IRL_DAP_Sb_BS<-read.table("Sb_BS_DGF_annotated_known_and_de_novo.bed")
IRL_DAP_Si_WL<-read.table("Si_WL_DGF_annotated_known_and_de_novo.bed")
IRL_DAP_Si_BS<-read.table("Si_BS_DGF_annotated_known_and_de_novo.bed")
IRL_DAP_Zm_WL<-read.table("Zm_WL_DGF_annotated_known_and_de_novo.bed")
IRL_DAP_Zm_BS<-read.table("Zm_BS_DGF_annotated_known_and_de_novo.bed")
IRL_DAP_Bd_WL<-read.table("Bd_WL_DGF_annotated_known_and_de_novo.bed")

#Read in all TF binding sites known
Si_WL_DAP<-read.table("Si_WL_DGF_annotated_DAP.bed")
Si_BS_DAP<-read.table("Si_BS_DGF_annotated_DAP.bed")
Sb_WL_DAP<-read.table("Sb_WL_DGF_annotated_DAP.bed")
Sb_BS_DAP<-read.table("Sb_BS_DGF_annotated_DAP.bed")
Zm_WL_DAP<-read.table("Zm_WL_DGF_annotated_DAP.bed")
Zm_BS_DAP<-read.table("Zm_BS_DGF_annotated_DAP.bed")
Bd_WL_DAP<-read.table("Bd_WL_DGF_annotated_DAP.bed")

#Read in all de novo binding sites
Si_WL_IRL<-read.table("Si_WL_DGF_annotated_IRL.bed")
Si_BS_IRL<-read.table("Si_BS_DGF_annotated_IRL.bed")
Sb_WL_IRL<-read.table("Sb_WL_DGF_annotated_IRL.bed")
Sb_BS_IRL<-read.table("Sb_BS_DGF_annotated_IRL.bed")
Zm_WL_IRL<-read.table("Zm_WL_DGF_annotated_IRL.bed")
Zm_BS_IRL<-read.table("Zm_BS_DGF_annotated_IRL.bed")
Bd_WL_IRL<-read.table("Bd_WL_DGF_annotated_IRL.bed")

#Make a new column with data on position summarized for merging
Zm_WL_DAP$site<-paste(Zm_WL_DAP[,1],Zm_WL_DAP[,2],Zm_WL_DAP[,3],sep=":")
Zm_BS_DAP$site<-paste(Zm_BS_DAP[,1],Zm_BS_DAP[,2],Zm_BS_DAP[,3],sep=":")
Si_WL_DAP$site<-paste(Si_WL_DAP[,1],Si_WL_DAP[,2],Si_WL_DAP[,3],sep=":")
Si_BS_DAP$site<-paste(Si_BS_DAP[,1],Si_BS_DAP[,2],Si_BS_DAP[,3],sep=":")
Sb_WL_DAP$site<-paste(Sb_WL_DAP[,1],Sb_WL_DAP[,2],Sb_WL_DAP[,3],sep=":")
Sb_BS_DAP$site<-paste(Sb_BS_DAP[,1],Sb_BS_DAP[,2],Sb_BS_DAP[,3],sep=":")
Bd_WL_DAP$site<-paste(Bd_WL_DAP[,1],Bd_WL_DAP[,2],Bd_WL_DAP[,3],sep=":")
Zm_WL_IRL$site<-paste(Zm_WL_IRL[,1],Zm_WL_IRL[,2],Zm_WL_IRL[,3],sep=":")
Zm_BS_IRL$site<-paste(Zm_BS_IRL[,1],Zm_BS_IRL[,2],Zm_BS_IRL[,3],sep=":")
Si_WL_IRL$site<-paste(Si_WL_IRL[,1],Si_WL_IRL[,2],Si_WL_IRL[,3],sep=":")
Si_BS_IRL$site<-paste(Si_BS_IRL[,1],Si_BS_IRL[,2],Si_BS_IRL[,3],sep=":")
Sb_WL_IRL$site<-paste(Sb_WL_IRL[,1],Sb_WL_IRL[,2],Sb_WL_IRL[,3],sep=":")
Sb_BS_IRL$site<-paste(Sb_BS_IRL[,1],Sb_BS_IRL[,2],Sb_BS_IRL[,3],sep=":")
Bd_WL_IRL$site<-paste(Bd_WL_IRL[,1],Bd_WL_IRL[,2],Bd_WL_IRL[,3],sep=":")
BED_DGF_Zm_WL$site<-paste(BED_DGF_Zm_WL[,1],BED_DGF_Zm_WL[,2],BED_DGF_Zm_WL[,3],sep=":")
BED_DGF_Zm_BS$site<-paste(BED_DGF_Zm_BS[,1],BED_DGF_Zm_BS[,2],BED_DGF_Zm_BS[,3],sep=":")
BED_DGF_Si_WL$site<-paste(BED_DGF_Si_WL[,1],BED_DGF_Si_WL[,2],BED_DGF_Si_WL[,3],sep=":")
BED_DGF_Si_BS$site<-paste(BED_DGF_Si_BS[,1],BED_DGF_Si_BS[,2],BED_DGF_Si_BS[,3],sep=":")
BED_DGF_Sb_WL$site<-paste(BED_DGF_Sb_WL[,1],BED_DGF_Sb_WL[,2],BED_DGF_Sb_WL[,3],sep=":")
BED_DGF_Sb_BS$site<-paste(BED_DGF_Sb_BS[,1],BED_DGF_Sb_BS[,2],BED_DGF_Sb_BS[,3],sep=":")
BED_DGF_Bd_WL$site<-paste(BED_DGF_Bd_WL[,1],BED_DGF_Bd_WL[,2],BED_DGF_Bd_WL[,3],sep=":")
#keep only the combined summary for merging
BED_DGF_Zm_WL<-unique(BED_DGF_Zm_WL[,"site"])
BED_DGF_Zm_BS<-unique(BED_DGF_Zm_BS[,"site"])
BED_DGF_Sb_WL<-unique(BED_DGF_Sb_WL[,"site"])
BED_DGF_Sb_BS<-unique(BED_DGF_Sb_BS[,"site"])
BED_DGF_Si_WL<-unique(BED_DGF_Si_WL[,"site"])
BED_DGF_Si_BS<-unique(BED_DGF_Si_BS[,"site"])
BED_DGF_Bd_WL<-unique(BED_DGF_Bd_WL[,"site"])
#Create a summary of known and de novo hits
Zm_WL_all_hits<-rbind(Zm_WL_DAP,Zm_WL_IRL)
Zm_BS_all_hits<-rbind(Zm_BS_DAP,Zm_BS_IRL)
Si_WL_all_hits<-rbind(Si_WL_DAP,Si_WL_IRL)
Si_BS_all_hits<-rbind(Si_BS_DAP,Si_BS_IRL)
Sb_WL_all_hits<-rbind(Sb_WL_DAP,Sb_WL_IRL)
Sb_BS_all_hits<-rbind(Sb_BS_DAP,Sb_BS_IRL)
Bd_WL_all_hits<-rbind(Bd_WL_DAP,Bd_WL_IRL)
#make sure there are no duplicates
Sb_BS_all_hits<-unique(Sb_BS_all_hits[,c("V4","site")])
Sb_WL_all_hits<-unique(Sb_WL_all_hits[,c("V4","site")])
Si_BS_all_hits<-unique(Si_BS_all_hits[,c("V4","site")])
Si_WL_all_hits<-unique(Si_WL_all_hits[,c("V4","site")])
Zm_BS_all_hits<-unique(Zm_BS_all_hits[,c("V4","site")])
Zm_WL_all_hits<-unique(Zm_WL_all_hits[,c("V4","site")])
Bd_WL_all_hits<-unique(Bd_WL_all_hits[,c("V4","site")])
Sb_BS_DAP<-Sb_BS_DAP[!duplicated(Sb_BS_DAP),]
Sb_WL_DAP<-Sb_WL_DAP[!duplicated(Sb_WL_DAP),]
Si_BS_DAP<-Si_BS_DAP[!duplicated(Si_BS_DAP),]
Si_WL_DAP<-Si_WL_DAP[!duplicated(Si_WL_DAP),]
Zm_BS_DAP<-Zm_BS_DAP[!duplicated(Zm_BS_DAP),]
Zm_WL_DAP<-Zm_WL_DAP[!duplicated(Zm_WL_DAP),]
Bd_WL_DAP<-Bd_WL_DAP[!duplicated(Bd_WL_DAP),]
Sb_BS_IRL<-Sb_BS_IRL[!duplicated(Sb_BS_IRL),]
Sb_WL_IRL<-Sb_WL_IRL[!duplicated(Sb_WL_IRL),]
Si_BS_IRL<-Si_BS_IRL[!duplicated(Si_BS_IRL),]
Si_WL_IRL<-Si_WL_IRL[!duplicated(Si_WL_IRL),]
Zm_BS_IRL<-Zm_BS_IRL[!duplicated(Zm_BS_IRL),]
Zm_WL_IRL<-Zm_WL_IRL[!duplicated(Zm_WL_IRL),]
Bd_WL_IRL<-Bd_WL_IRL[!duplicated(Bd_WL_IRL),]

#calculate how many DGFs are annotated
Sb_WL_covered<-length(unique(subset(Sb_WL_all_hits,Sb_WL_all_hits$site %in% BED_DGF_Sb_WL)["site"])[,1])
Sb_BS_covered<-length(unique(subset(Sb_BS_all_hits,Sb_BS_all_hits$site %in% BED_DGF_Sb_BS)["site"])[,1])
Si_WL_covered<-length(unique(subset(Si_WL_all_hits,Si_WL_all_hits$site %in% BED_DGF_Si_WL)["site"])[,1])
Si_BS_covered<-length(unique(subset(Si_BS_all_hits,Si_BS_all_hits$site %in% BED_DGF_Si_BS)["site"])[,1])
Zm_WL_covered<-length(unique(subset(Zm_WL_all_hits,Zm_WL_all_hits$site %in% BED_DGF_Zm_WL)["site"])[,1])
Zm_BS_covered<-length(unique(subset(Zm_BS_all_hits,Zm_BS_all_hits$site %in% BED_DGF_Zm_BS)["site"])[,1])
Bd_WL_covered<-length(unique(subset(Bd_WL_all_hits,Bd_WL_all_hits$site %in% BED_DGF_Bd_WL)["site"])[,1])
Sb_WL_covered_DAP<-length(unique(subset(Sb_WL_DAP,Sb_WL_DAP$site %in% BED_DGF_Sb_WL)["site"])[,1])
Sb_BS_covered_DAP<-length(unique(subset(Sb_BS_DAP,Sb_BS_DAP$site %in% BED_DGF_Sb_BS)["site"])[,1])
Si_WL_covered_DAP<-length(unique(subset(Si_WL_DAP,Si_WL_DAP$site %in% BED_DGF_Si_WL)["site"])[,1])
Si_BS_covered_DAP<-length(unique(subset(Si_BS_DAP,Si_BS_DAP$site %in% BED_DGF_Si_BS)["site"])[,1])
Zm_WL_covered_DAP<-length(unique(subset(Zm_WL_DAP,Zm_WL_DAP$site %in% BED_DGF_Zm_WL)["site"])[,1])
Zm_BS_covered_DAP<-length(unique(subset(Zm_BS_DAP,Zm_BS_DAP$site %in% BED_DGF_Zm_BS)["site"])[,1])
Bd_WL_covered_DAP<-length(unique(subset(Bd_WL_DAP,Bd_WL_DAP$site %in% BED_DGF_Bd_WL)["site"])[,1])
Sb_WL_covered_IRL<-length(unique(subset(Sb_WL_IRL,Sb_WL_IRL$site %in% BED_DGF_Sb_WL)["site"])[,1])
Sb_BS_covered_IRL<-length(unique(subset(Sb_BS_IRL,Sb_BS_IRL$site %in% BED_DGF_Sb_BS)["site"])[,1])
Si_WL_covered_IRL<-length(unique(subset(Si_WL_IRL,Si_WL_IRL$site %in% BED_DGF_Si_WL)["site"])[,1])
Si_BS_covered_IRL<-length(unique(subset(Si_BS_IRL,Si_BS_IRL$site %in% BED_DGF_Si_BS)["site"])[,1])
Zm_WL_covered_IRL<-length(unique(subset(Zm_WL_IRL,Zm_WL_IRL$site %in% BED_DGF_Zm_WL)["site"])[,1])
Zm_BS_covered_IRL<-length(unique(subset(Zm_BS_IRL,Zm_BS_IRL$site %in% BED_DGF_Zm_BS)["site"])[,1])
Bd_WL_covered_IRL<-length(unique(subset(Bd_WL_IRL,Bd_WL_IRL$site %in% BED_DGF_Bd_WL)["site"])[,1])
Sb_WL_covered_IRL_overlap<-length(unique(subset(Sb_WL_IRL,Sb_WL_IRL$site %in% Sb_WL_DAP$site)["site"])[,1])
Sb_BS_covered_IRL_overlap<-length(unique(subset(Sb_BS_IRL,Sb_BS_IRL$site %in% Sb_BS_DAP$site)["site"])[,1])
Si_WL_covered_IRL_overlap<-length(unique(subset(Si_WL_IRL,Si_WL_IRL$site %in% Si_WL_DAP$site)["site"])[,1])
Si_BS_covered_IRL_overlap<-length(unique(subset(Si_BS_IRL,Si_BS_IRL$site %in% Si_BS_DAP$site)["site"])[,1])
Zm_WL_covered_IRL_overlap<-length(unique(subset(Zm_WL_IRL,Zm_WL_IRL$site %in% Zm_WL_DAP$site)["site"])[,1])
Zm_BS_covered_IRL_overlap<-length(unique(subset(Zm_BS_IRL,Zm_BS_IRL$site %in% Zm_BS_DAP$site)["site"])[,1])
Bd_WL_covered_IRL_overlap<-length(unique(subset(Bd_WL_IRL,Bd_WL_IRL$site %in% Bd_WL_DAP$site)["site"])[,1])

#Calculate the total number of DGFs
num_Sb_WL_all<-length(unique(BED_DGF_Sb_WL))
num_Sb_BS_all<-length(unique(BED_DGF_Sb_BS))
num_Si_WL_all<-length(unique(BED_DGF_Si_WL))
num_Si_BS_all<-length(unique(BED_DGF_Si_BS))
num_Zm_WL_all<-length(unique(BED_DGF_Zm_WL))
num_Zm_BS_all<-length(unique(BED_DGF_Zm_BS))
num_Bd_WL_all<-length(unique(BED_DGF_Bd_WL))


#Create a summary table
names<-c("Sb_WL",
         "Sb_BS",
         "Si_WL",
         "Si_BS",
         "Zm_WL",
         "Zm_BS",
         "Bd_WL")
num_all<-c(num_Sb_WL_all,
           num_Sb_BS_all,
           num_Si_WL_all,
           num_Si_BS_all,
           num_Zm_WL_all,
           num_Zm_BS_all,
           num_Bd_WL_all)
num_covered<-c(Sb_WL_covered,
             Sb_BS_covered,
             Si_WL_covered,
             Si_BS_covered,
             Zm_WL_covered,
             Zm_BS_covered,
             Bd_WL_covered)
num_known<-c(Sb_WL_covered_DAP,
           Sb_BS_covered_DAP,
           Si_WL_covered_DAP,
           Si_BS_covered_DAP,
           Zm_WL_covered_DAP,
           Zm_BS_covered_DAP,
           Bd_WL_covered_DAP)
num_de_novo<-c(Sb_WL_covered_IRL,
             Sb_BS_covered_IRL,
             Si_WL_covered_IRL,
             Si_BS_covered_IRL,
             Zm_WL_covered_IRL,
             Zm_BS_covered_IRL,
             Bd_WL_covered_IRL)
num_known_and_de_novo<-c(Sb_WL_covered_IRL_overlap,
               Sb_BS_covered_IRL_overlap,
               Si_WL_covered_IRL_overlap,
               Si_BS_covered_IRL_overlap,
               Zm_WL_covered_IRL_overlap,
               Zm_BS_covered_IRL_overlap,
               Bd_WL_covered_IRL_overlap)

res<-data.frame("Sample"=names,
           "All_DGF"=num_all,
           "DGF_annotated"=num_covered,
           "DGF_hit_known"=num_known,
           "DGF_hit_de_novo"=num_de_novo,
           "DGF_hit_known_and_de_novo"=num_known_and_de_novo)

res$DGF_unannotated<-res$All_DGF-res$DGF_annotated
res$DGF_de_novo_only<-res$DGF_annotated-res$DGF_hit_known


#Calculate DGFs annotated by percent
p<-data.frame("Sample"=res$Sample)
p$known<-(res$DGF_hit_known/res$All_DGF)*100
p$de_novo_only<-(res$DGF_de_novo_only/res$All_DGF)*100
p$unannotated<-(res$DGF_unannotated/res$All_DGF)*100

#res_plot<-data.frame(t(p[,c(2:4)]))
#colnames(res_plot)<-colnames(res_plot)<-res$Sample

write.table(p,"Data_Figure_4C.txt",quote=F)
