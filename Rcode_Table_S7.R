#R script to generate Table S7

#Read in Annotated DE DGFs
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/")
DGF_Sb_BS_DE<-read.table("DGF_DE_Sb_BS_annotated_locus.bed",head=T)
DGF_Si_BS_DE<-read.table("DGF_DE_Si_BS_annotated_locus.bed",head=T)
DGF_Zm_BS_DE<-read.table("DGF_DE_Zm_BS_annotated_locus.bed",head=T)
DGF_Sb_WL_DE<-read.table("DGF_DE_Sb_WL_annotated_locus.bed",head=T)
DGF_Si_WL_DE<-read.table("DGF_DE_Si_WL_annotated_locus.bed",head=T)
DGF_Zm_WL_DE<-read.table("DGF_DE_Zm_WL_annotated_locus.bed",head=T)

#Renames columns for downstream processing
colnames(DGF_Sb_BS_DE)<-c("Chr","Start","End","Gene")
colnames(DGF_Sb_WL_DE)<-c("Chr","Start","End","Gene")
colnames(DGF_Si_BS_DE)<-c("Chr","Start","End","Gene")
colnames(DGF_Si_WL_DE)<-c("Chr","Start","End","Gene")
colnames(DGF_Zm_BS_DE)<-c("Chr","Start","End","Gene")
colnames(DGF_Zm_WL_DE)<-c("Chr","Start","End","Gene")

#Read in DE Gene lists
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_DE")
Gene_Sb_BS_DE<-data.frame("Gene"=read.table("Sb_DE_FDR_01_BS.csv",head=T)$ID,"Cell"="BS")
Gene_Sb_M_DE<-data.frame("Gene"=read.table("Sb_DE_FDR_01_M.csv",head=F)[,1],"Cell"="M")
Gene_Si_BS_DE<-data.frame("Gene"=as.factor(paste(read.table("Si_DE_FDR_01_BS.csv",sep=",",head=T)$gene,".g",sep="")),"Cell"="BS")
Gene_Si_M_DE<-data.frame("Gene"=as.factor(paste(read.table("Si_DE_FDR_01_M.csv",sep=",",head=T)$gene,".g",sep="")),"Cell"="M")
Gene_Zm_BS_DE<-data.frame("Gene"=read.table("Zm_DE_FDR_01_BS.csv",sep=",",head=T)$gene,"Cell"="BS")
Gene_Zm_M_DE<-data.frame("Gene"=read.table("Zm_DE_FDR_01_M.csv",sep=",",head=T)$gene,"Cell"="M")

CalcStats<-function(Gene_BS,Gene_M,DGF_BS,DGF_WL,Species){
  
  #First create a dataframe of all DE genes
  Genes_all_DE<-rbind(Gene_BS,Gene_M)
  #Create a dataframe of all DE footprints
  DGF_all_DE<-unique(rbind(DGF_BS,DGF_WL))

  #Make a dataframe with all DE footprints in DE genes
  Genes_DE_with_DE_DGF<-unique(merge(Genes_all_DE,DGF_all_DE,by="Gene"))

  #Calculate the number of DE Genes
  num_Gene_DE<-length(unique(Genes_all_DE$Gene))
  #Calculate the number DE DGFs
  num_DGF_DE<-length(DGF_all_DE[,1])
  #Calculate the number DE DGFs associated with DE genes
  num_DGF_DE_in_DE_Genes<-length(Genes_DE_with_DE_DGF[,1])
  #Calculate the number of DE genes with at least one DE footprint
  num_Gene_DE_with_DE_DGF<-length(unique(Genes_DE_with_DE_DGF$Gene))
  #Caluclate Percentages
  percent_Gene_DE_with_DE_DGF<-round((num_Gene_DE_with_DE_DGF/num_Gene_DE)*100)
  percent_DGF_DE_in_DE_Genes<-round((num_DGF_DE_in_DE_Genes/num_DGF_DE)*100)
  
  #Create a Results Vector
  Results<-data.frame("Species"=Species,
           "Number.of.DE.DGFs"=num_DGF_DE,
           "Number.of.DE.Genes"=num_Gene_DE,
           "Number.of.DE.Genes.with.DE.DGF"=num_Gene_DE_with_DE_DGF,
           "Number.of.DE.DGF.Associated.With.DE.Genes"=num_DGF_DE_in_DE_Genes,
           "Percent.of.DE.Genes.with.DE.DGF"=percent_Gene_DE_with_DE_DGF,
           "Percent.of.DE.DGF.in.DE.Genes"=percent_DGF_DE_in_DE_Genes)

 return(Results)
}

Si_Results<-CalcStats(Gene_Si_BS_DE,Gene_Si_M_DE,DGF_Si_BS_DE,DGF_Si_WL_DE,"S. italica")
Sb_Results<-CalcStats(Gene_Sb_BS_DE,Gene_Sb_M_DE,DGF_Sb_BS_DE,DGF_Sb_WL_DE,"S. bicolor")
Zm_Results<-CalcStats(Gene_Zm_BS_DE,Gene_Zm_M_DE,DGF_Zm_BS_DE,DGF_Zm_WL_DE,"Z. mays")

#Write Results to file
Out_table<-rbind(Sb_Results,Zm_Results,Si_Results)
write.table(Out_table,"/home/sjb287/Table - S7.csv",quote=F,row.name=F,sep=",")