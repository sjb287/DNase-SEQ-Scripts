#R code to make Table S10
#################################################
annotate<-function(inData,ID){
  temp<-merge(Orthologs,inData,by=ID)
  temp<-temp[!duplicated(temp),]
  temp<-temp[,c("Sb_ID","Si_ID","Zm_ID","Bd_ID","Motif")]
  return(temp)
}
#################################################
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/")
#Read in files
Sb_BS<-read.table("Sb_BS_DGF_annotated_known_and_de_novo.bed")
Si_BS<-read.table("Si_BS_DGF_annotated_known_and_de_novo.bed")
Zm_BS<-read.table("Zm_BS_DGF_annotated_known_and_de_novo.bed")
Sb_WL<-read.table("Sb_WL_DGF_annotated_known_and_de_novo.bed")
Si_WL<-read.table("Si_WL_DGF_annotated_known_and_de_novo.bed")
Zm_WL<-read.table("Zm_WL_DGF_annotated_known_and_de_novo.bed")
Bd_WL<-read.table("Bd_WL_DGF_annotated_known_and_de_novo.bed")
#Rename Columns for merging
colnames(Sb_BS)<-c("Chr","Start","Stop","Motif","Sb_ID")
colnames(Sb_WL)<-c("Chr","Start","Stop","Motif","Sb_ID")
colnames(Si_BS)<-c("Chr","Start","Stop","Motif","Si_ID")
colnames(Si_WL)<-c("Chr","Start","Stop","Motif","Si_ID")
colnames(Zm_BS)<-c("Chr","Start","Stop","Motif","Zm_ID")
colnames(Zm_WL)<-c("Chr","Start","Stop","Motif","Zm_ID")
colnames(Bd_WL)<-c("Chr","Start","Stop","Motif","Bd_ID")
Orthologs<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/Cross_mapped_IDs_phytozome.csv",head=T,sep=" ")
#annotate with orthologs
Sb_BS_orthologs<-annotate(Sb_BS,"Sb_ID")
Sb_WL_orthologs<-annotate(Sb_WL,"Sb_ID")
Si_BS_orthologs<-annotate(Si_BS,"Si_ID")
Si_WL_orthologs<-annotate(Si_WL,"Si_ID")
Zm_BS_orthologs<-annotate(Zm_BS,"Zm_ID")
Zm_WL_orthologs<-annotate(Zm_WL,"Zm_ID")
Bd_WL_orthologs<-annotate(Bd_WL,"Bd_ID")
#calculate lineage shared motifs
SbZm_BS<-merge(Sb_BS_orthologs,Zm_BS_orthologs)
SbZm_BS<-SbZm_BS[!duplicated(SbZm_BS),]
#calculate C4 shared motifs
SbZmSi_BS<-merge(SbZm_BS,Si_BS_orthologs)
SbZmSi_BS<-SbZmSi_BS[!duplicated(SbZmSi_BS),]
#calculate 4 way shared motifs
BdSbZmSi_BS<-merge(SbZmSi_BS,Bd_WL_orthologs)
BdSbZmSi_BS<-BdSbZmSi_BS[!duplicated(BdSbZmSi_BS),]

#calculate lineage shared motifs
SbZm_WL<-merge(Sb_WL_orthologs,Zm_WL_orthologs)
SbZm_WL<-SbZm_WL[!duplicated(SbZm_WL),]
#calculate C4 shared motifs
SbZmSi_WL<-merge(SbZm_WL,Si_WL_orthologs)
SbZmSi_WL<-SbZmSi_WL[!duplicated(SbZmSi_WL),]
#calculate 4 way shared motifs
BdSbZmSi_WL<-merge(SbZmSi_WL,Bd_WL_orthologs)
BdSbZmSi_WL<-BdSbZmSi_WL[!duplicated(BdSbZmSi_WL),]


#################################################
#Extract C4 Gene Data
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

WL_3_WAY<-subset(SbZmSi_WL,SbZmSi_WL$Zm_ID %in% C4$Zm_ID)
BS_3_WAY<-subset(SbZmSi_BS,SbZmSi_BS$Zm_ID %in% C4$Zm_ID)
BS_3_WAY<-merge(BS_3_WAY,Sb_BS)

WL_4_WAY<-subset(BdSbZmSi_WL,BdSbZmSi_WL$Zm_ID %in% C4$Zm_ID)
BS_4_WAY<-subset(BdSbZmSi_BS,BdSbZmSi_BS$Zm_ID %in% C4$Zm_ID)

#Get lists of all motifs in C4 genes
Sb_BS_C4_motifs<-subset(Sb_BS,Sb_BS$Sb_ID %in% C4$Sb_ID)
Si_BS_C4_motifs<-subset(Si_BS,Si_BS$Si_ID %in% C4$Si_ID)
Zm_BS_C4_motifs<-subset(Zm_BS,Zm_BS$Zm_ID %in% C4$Zm_ID)
colnames(Sb_BS_C4_motifs)[5]<-"ID"
colnames(Si_BS_C4_motifs)[5]<-"ID"
colnames(Zm_BS_C4_motifs)[5]<-"ID"
C4_BS_motifs<-rbind(Sb_BS_C4_motifs,Si_BS_C4_motifs,Zm_BS_C4_motifs)

Sb_WL_C4_motifs<-subset(Sb_WL,Sb_WL$Sb_ID %in% C4$Sb_ID)
Si_WL_C4_motifs<-subset(Si_WL,Si_WL$Si_ID %in% C4$Si_ID)
Zm_WL_C4_motifs<-subset(Zm_WL,Zm_WL$Zm_ID %in% C4$Zm_ID)
Bd_WL_C4_motifs<-subset(Bd_WL,Bd_WL$Bd_ID %in% C4$Bd_ID)


lineage_BS_motifs<-subset(Sb_BS_C4_motifs,Sb_BS_C4_motifs$Motif %in% Zm_BS_C4_motifs$Motif)
C4_BS_motifs<-subset(lineage_BS_motifs,lineage_BS_motifs$Motif %in% Si_BS_C4_motifs$Motif)
BS_motifs<-subset(C4_BS_motifs,C4_BS_motifs$Motif %in% Bd_WL_C4_motifs$Motif)