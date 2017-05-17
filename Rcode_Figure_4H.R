#Rcode to make Figure 4H
#################################################
extractFamilyData<-function(ID){
  #DE
  Sb_DE_WL<-unique(subset(Sb_WL_DAP_DE_Families,Sb_WL_DAP_DE_Families[,5] %in% C4$Sb_ID & Sb_WL_DAP_DE_Families[,6] ==ID))
  Si_DE_WL<-unique(subset(Si_WL_DAP_DE_Families,Si_WL_DAP_DE_Families[,5] %in% C4$Si_ID & Si_WL_DAP_DE_Families[,6] ==ID))
  Zm_DE_WL<-unique(subset(Zm_WL_DAP_DE_Families,Zm_WL_DAP_DE_Families[,5] %in% C4$Zm_ID & Zm_WL_DAP_DE_Families[,6] ==ID))
  Sb_DE_BS<-unique(subset(Sb_BS_DAP_DE_Families,Sb_BS_DAP_DE_Families[,5] %in% C4$Sb_ID & Sb_BS_DAP_DE_Families[,6] ==ID))
  Si_DE_BS<-unique(subset(Si_BS_DAP_DE_Families,Si_BS_DAP_DE_Families[,5] %in% C4$Si_ID & Si_BS_DAP_DE_Families[,6] ==ID))
  Zm_DE_BS<-unique(subset(Zm_BS_DAP_DE_Families,Zm_BS_DAP_DE_Families[,5] %in% C4$Zm_ID & Zm_BS_DAP_DE_Families[,6] ==ID))
  #nonDE
  Sb_WL<-unique(subset(Sb_WL_DAP_Families,Sb_WL_DAP_Families[,5] %in% C4$Sb_ID & Sb_WL_DAP_Families[,6] ==ID))
  Si_WL<-unique(subset(Si_WL_DAP_Families,Si_WL_DAP_Families[,5] %in% C4$Si_ID & Si_WL_DAP_Families[,6] ==ID))
  Zm_WL<-unique(subset(Zm_WL_DAP_Families,Zm_WL_DAP_Families[,5] %in% C4$Zm_ID & Zm_WL_DAP_Families[,6] ==ID))
  Sb_BS<-unique(subset(Sb_BS_DAP_Families,Sb_BS_DAP_Families[,5] %in% C4$Sb_ID & Sb_BS_DAP_Families[,6] ==ID))
  Si_BS<-unique(subset(Si_BS_DAP_Families,Si_BS_DAP_Families[,5] %in% C4$Si_ID & Si_BS_DAP_Families[,6] ==ID))
  Zm_BS<-unique(subset(Zm_BS_DAP_Families,Zm_BS_DAP_Families[,5] %in% C4$Zm_ID & Zm_BS_DAP_Families[,6] ==ID))

  #annotate with gene ID
  C4_Sb<-C4[,c(2,4)]
  C4_Si<-C4[,c(2,3)]
  C4_Zm<-C4[,c(2,1)]
  colnames(C4_Sb)<-c("C4","Locus")
  colnames(C4_Si)<-c("C4","Locus")
  colnames(C4_Zm)<-c("C4","Locus")
  Sb_DE_WL<-merge(Sb_DE_WL,C4_Sb,by="Locus")
  Si_DE_WL<-merge(Si_DE_WL,C4_Si,by="Locus")
  Zm_DE_WL<-merge(Zm_DE_WL,C4_Zm,by="Locus")
  Sb_DE_BS<-merge(Sb_DE_BS,C4_Sb,by="Locus")
  Si_DE_BS<-merge(Si_DE_BS,C4_Si,by="Locus")
  Zm_DE_BS<-merge(Zm_DE_BS,C4_Zm,by="Locus")
  Sb_WL<-merge(Sb_WL,C4_Sb,by="Locus")
  Si_WL<-merge(Si_WL,C4_Si,by="Locus")
  Zm_WL<-merge(Zm_WL,C4_Zm,by="Locus")
  Sb_BS<-merge(Sb_BS,C4_Sb,by="Locus")
  Si_BS<-merge(Si_BS,C4_Si,by="Locus")
  Zm_BS<-merge(Zm_BS,C4_Zm,by="Locus")
  
  res<-list(Sb_WL,Si_WL,Zm_WL,Sb_BS,Si_BS,Zm_BS,Sb_DE_WL,Si_DE_WL,Zm_DE_WL,Sb_DE_BS,Si_DE_BS,Zm_DE_BS)
  return(res)
}

#################################################
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/")

#Read in files
#DE DGF only
Sb_BS_DAP_DE_Families<-read.table("Sb_BS_DGF_DE_annotated_DAP.bed")
Si_BS_DAP_DE_Families<-read.table("Si_BS_DGF_DE_annotated_DAP.bed")
Zm_BS_DAP_DE_Families<-read.table("Zm_BS_DGF_DE_annotated_DAP.bed")
Sb_WL_DAP_DE_Families<-read.table("Sb_WL_DGF_DE_annotated_DAP.bed")
Si_WL_DAP_DE_Families<-read.table("Si_WL_DGF_DE_annotated_DAP.bed")
Zm_WL_DAP_DE_Families<-read.table("Zm_WL_DGF_DE_annotated_DAP.bed")
#non-DE DGF
Sb_BS_DAP_Families<-read.table("Sb_BS_DGF_annotated_DAP.bed")
Si_BS_DAP_Families<-read.table("Si_BS_DGF_annotated_DAP.bed")
Zm_BS_DAP_Families<-read.table("Zm_BS_DGF_annotated_DAP.bed")
Sb_WL_DAP_Families<-read.table("Sb_WL_DGF_annotated_DAP.bed")
Si_WL_DAP_Families<-read.table("Si_WL_DGF_annotated_DAP.bed")
Zm_WL_DAP_Families<-read.table("Zm_WL_DGF_annotated_DAP.bed")
Bd_WL_DAP_Families<-read.table("Bd_WL_DGF_annotated_DAP.bed")
##################################################
#Look at enrichment of Motif Family Results
##################################################
#Calculate families
DAP_annotation<-read.table("/home/sjb287/documents/Motif_Searching/DAP_Motifs/DAP_motif_annotation.txt")
DAP_annotation<-DAP_annotation[,c(1,3)]
colnames(DAP_annotation)<-c("Motif","Family")

#Rename Columns for merging
colnames(Sb_BS_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_WL_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_BS_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_WL_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_BS_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_WL_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_BS_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_BS_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_BS_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Bd_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")

#Merge individual motifs with TF family annotations
Sb_BS_DAP_DE_Families<-unique(merge(Sb_BS_DAP_DE_Families,DAP_annotation,by="Motif"))
Sb_WL_DAP_DE_Families<-unique(merge(Sb_WL_DAP_DE_Families,DAP_annotation,by="Motif"))
Si_BS_DAP_DE_Families<-unique(merge(Si_BS_DAP_DE_Families,DAP_annotation,by="Motif"))
Si_WL_DAP_DE_Families<-unique(merge(Si_WL_DAP_DE_Families,DAP_annotation,by="Motif"))
Zm_BS_DAP_DE_Families<-unique(merge(Zm_BS_DAP_DE_Families,DAP_annotation,by="Motif"))
Zm_WL_DAP_DE_Families<-unique(merge(Zm_WL_DAP_DE_Families,DAP_annotation,by="Motif"))
Sb_WL_DAP_Families<-unique(merge(Sb_WL_DAP_Families,DAP_annotation,by="Motif"))
Si_WL_DAP_Families<-unique(merge(Si_WL_DAP_Families,DAP_annotation,by="Motif"))
Zm_WL_DAP_Families<-unique(merge(Zm_WL_DAP_Families,DAP_annotation,by="Motif"))
Sb_BS_DAP_Families<-unique(merge(Sb_BS_DAP_Families,DAP_annotation,by="Motif"))
Si_BS_DAP_Families<-unique(merge(Si_BS_DAP_Families,DAP_annotation,by="Motif"))
Zm_BS_DAP_Families<-unique(merge(Zm_BS_DAP_Families,DAP_annotation,by="Motif"))
Bd_WL_DAP_Families<-unique(merge(Bd_WL_DAP_Families,DAP_annotation,by="Motif"))

Sb_BS_DAP_DE_Families$Cell<-"BS_DE"
Sb_WL_DAP_DE_Families$Cell<-"WL_DE"
Si_BS_DAP_DE_Families$Cell<-"BS_DE"
Si_WL_DAP_DE_Families$Cell<-"WL_DE"
Zm_BS_DAP_DE_Families$Cell<-"BS_DE"
Zm_WL_DAP_DE_Families$Cell<-"WL_DE"
Sb_WL_DAP_Families$Cell<-"WL"
Si_WL_DAP_Families$Cell<-"WL"
Zm_WL_DAP_Families$Cell<-"WL"
Sb_BS_DAP_Families$Cell<-"BS"
Si_BS_DAP_Families$Cell<-"BS"
Zm_BS_DAP_Families$Cell<-"BS"
Bd_WL_DAP_Families$Cell<-"WL"

#################################################
#Extract C4 Genes annotated with C2C2-GATA TF Family
#################################################
#Read In annotation Files
#Note - instances where phytozome mapping did not provide an ortholog were corrected by manual addition of the nearest matching sequence
Si_C4_genes<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_Calvin_Zm_Si.csv",head=T,sep=",")
Sb_C4_genes<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_Calvin_Zm_Sb.csv",head=T,sep="\t")
Bd_C4_genes<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_Calvin_Zm_Bd.csv",head=T,sep="\t")
C4_IDs<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/Zm_C4_gene_IDs.csv",head=T,sep=",")
Si_C4_genes<-Si_C4_genes[,c(1,3)]
Sb_C4_genes<-Sb_C4_genes[,c(1,3)]
Bd_C4_genes<-Bd_C4_genes[,c(1,3)]
#Create Final annotation dataframe
colnames(Si_C4_genes)<-c("Zm_ID","Si_ID")
colnames(Sb_C4_genes)<-c("Zm_ID","Sb_ID")
colnames(Bd_C4_genes)<-c("Zm_ID","Bd_ID")
C4<-merge(C4_IDs,Si_C4_genes,by="Zm_ID")
C4<-merge(C4,Sb_C4_genes,by="Zm_ID")
C4<-merge(C4,Bd_C4_genes,by="Zm_ID")
#Extract from all footprints
C4_Sb<-unique(subset(Sb_WL_DAP_Families,Sb_WL_DAP_Families[,5] %in% C4$Sb_ID))
C4_Si<-unique(subset(Si_WL_DAP_Families,Si_WL_DAP_Families[,5] %in% C4$Si_ID))
C4_Zm<-unique(subset(Zm_WL_DAP_Families,Zm_WL_DAP_Families[,5] %in% C4$Zm_ID))

GATA<-extractFamilyData("C2C2-GATA")
BZR<-extractFamilyData("BZR")
TCP<-extractFamilyData("TCP")
bZIP<-extractFamilyData("bZIP")
bHLH<-extractFamilyData("bHLH")
HB<-extractFamilyData("HB")
LIM<-extractFamilyData("LIM")

#################################################
#TF Family
#################################################
#Read in annotations
Zm_GATA_TF<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/Zm_GATA_family.csv",head=T,sep="\t")
All_IDs<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/Cross_mapped_IDs_phytozome.csv",head=T,sep=" ")
GATA_TF<-merge(Zm_GATA_TF,All_IDs,by="Zm_ID")

######################
#Read in Sequence data
######################
Sb<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_DE/BaySEQ_Data/Sb_TPMs_DE.csv",sep=",",head=T)
Sb_Syn<-read.table("/data/genomes/sorghum/Sbicolor_255_v2.1.synonym_loci.txt")
colnames(Sb_Syn)<-c("ID","gene")
#re-annotate file with new IDs
Sb_TPMs<-merge(Sb_Syn,Sb,by="gene")
Sb_TPMs<-Sb_TPMs[,c(2:length(Sb_TPMs))]
colnames(Sb_TPMs)[1]<-"gene"
Si_TPMs<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_DE/BaySEQ_Data/Si_TPMs_DE_corrected.csv",sep=",",head=T)
Si_Syn<-read.table("/data/genomes/setaria/Sitalica/Sitalica_312_v2.2.synonym_loci.txt")
colnames(Si_Syn)<-c("ID","gene")
Si_Syn$gene<-paste(Si_Syn$gene,".g",sep="")
#re-annotate file with new IDs
Si_TPMs<-merge(Si_Syn,Si_TPMs,by="gene")
Si_TPMs<-Si_TPMs[,c(2:length(Si_TPMs))]
colnames(Si_TPMs)[1]<-"gene"
Zm_TPMs<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_DE/BaySEQ_Data/Zm_TPMs_DE.csv",sep=",",head=T)

#Subset genes for M and BS enrichment
Sb_BS_DE<-subset(Sb_TPMs,Sb_TPMs$DE=="bs>m" & Sb_TPMs$FDR.DE<0.01)$gene
Sb_M_DE<-subset(Sb_TPMs,Sb_TPMs$DE=="m>bs" & Sb_TPMs$FDR.DE<0.01)$gene
Si_BS_DE<-subset(Si_TPMs,Si_TPMs$DE=="bs>m" & Si_TPMs$FDR.DE<0.01)$gene
Si_M_DE<-subset(Si_TPMs,Si_TPMs$DE=="m>bs" & Si_TPMs$FDR.DE<0.01)$gene
Zm_BS_DE<-subset(Zm_TPMs,Zm_TPMs$DE=="bs>m" & Zm_TPMs$FDR.DE<0.01)$gene
Zm_M_DE<-subset(Zm_TPMs,Zm_TPMs$DE=="m>bs" & Zm_TPMs$FDR.DE<0.01)$gene
#Identify TF orthologs with consistent expression pattern
subset(GATA_TF,GATA_TF$Zm_ID %in% Zm_M_DE & GATA_TF$Sb_ID %in% Sb_M_DE & GATA_TF$Si_ID %in% Si_M_DE)
#Get expression data for all members of the family
GATA_Si<-subset(Si_TPMs,Si_TPMs$gene %in% GATA_TF$Si_ID)
GATA_Sb<-subset(Sb_TPMs,Sb_TPMs$gene %in% GATA_TF$Sb_ID)
GATA_Zm<-subset(Zm_TPMs,Zm_TPMs$gene %in% GATA_TF$Zm_ID)
#Calculate averages
GATA_Si$M_average<-apply(GATA_Si[,c(3,4,5)],1,mean)
GATA_Si$B_average<-apply(GATA_Si[,c(6,7,8)],1,mean)
GATA_Sb$M_average<-apply(GATA_Sb[,c(3,4,5)],1,mean)
GATA_Sb$B_average<-apply(GATA_Sb[,c(6,7,8)],1,mean)
GATA_Zm$M_average<-apply(GATA_Zm[,c(3,4)],1,mean)
GATA_Zm$B_average<-apply(GATA_Zm[,c(5,6)],1,mean)
#Plot Data
plot(x=GATA_Si$M_average,y=GATA_Si$B_average,xlim=c(0,60),ylim=c(0,60))
plot(x=GATA_Sb$M_average,y=GATA_Sb$B_average,xlim=c(0,60),ylim=c(0,60))
plot(x=GATA_Zm$M_average,y=GATA_Zm$B_average,xlim=c(0,60),ylim=c(0,60))

Z<-subset(GATA_Zm,GATA_Zm$gene=="GRMZM2G379005")
Sb<-subset(GATA_Sb,GATA_Sb$gene=="Sobic.004G337500")
Si<-subset(GATA_Si,GATA_Si$gene=="Seita.1G358400")
plot(x=Z$M_average,y=Z$B_average,xlim=c(0,60),ylim=c(0,60))
plot(x=Sb$M_average,y=Sb$B_average,xlim=c(0,60),ylim=c(0,60))
plot(x=Si$M_average,y=Si$B_average,xlim=c(0,60),ylim=c(0,60))

#350x400