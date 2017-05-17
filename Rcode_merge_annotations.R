
annotate_hits<-function(motif_file,annotation,outfile,ID){

  #read in motif annotations   
  temp_1<-read.table(motif_file,head=T)
  
  #rename columns for merging
  colnames(temp_1)<-c("Chr","start","stop","AGI")
  colnames(annotation)<-c("Chr","start","stop",ID)
  
  #merge DGF:Motif annotations with DGF:gene annotations
  #creating an annotation of DGF:Gene:motif 
  results_1<-merge(temp_1,annotation)
  
  #write results to file
  write.table(results_1,outfile,quote=F,row.name=F,col.names=F,sep="\t")
  
}

#Set destination of files where DGFs are annotated with motifs
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/")

#Files DGFs annotated with motifs
DAP_Sb_WL<-"DAP_Sb_WL/Sb_WL_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Sb_BS<-"DAP_Sb_BS/Sb_BS_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Si_WL<-"DAP_Si_WL/Si_WL_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Si_BS<-"DAP_Si_BS/Si_BS_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Zm_WL<-"DAP_Zm_WL/Zm_WL_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Zm_BS<-"DAP_Zm_BS/Zm_BS_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Bd_WL<-"DAP_Bd_WL/Bd_WL_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Sb_WL<-"IRL_Sb_WL/Sb_WL_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Sb_BS<-"IRL_Sb_BS/Sb_BS_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Si_WL<-"IRL_Si_WL/Si_WL_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Si_BS<-"IRL_Si_BS/Si_BS_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Zm_WL<-"IRL_Zm_WL/Zm_WL_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Zm_BS<-"IRL_Zm_BS/Zm_BS_DGF_fdr_0.01_fimo_hits_all.bed"
IRL_Bd_WL<-"IRL_Bd_WL/Bd_WL_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Sb_WL_DE<-"DGF_DE_annotated_Sb_WL/Sb_WL_DE_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Sb_BS_DE<-"DGF_DE_annotated_Sb_BS/Sb_BS_DE_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Si_WL_DE<-"DGF_DE_annotated_Si_WL/Si_WL_DE_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Si_BS_DE<-"DGF_DE_annotated_Si_BS/Si_BS_DE_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Zm_WL_DE<-"DGF_DE_annotated_Zm_WL/Zm_WL_DE_DGF_fdr_0.01_fimo_hits_all.bed"
DAP_Zm_BS_DE<-"DGF_DE_annotated_Zm_BS/Zm_BS_DE_DGF_fdr_0.01_fimo_hits_all.bed"



#Files with DGF annotated with genomic loci
Bd_WL_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Bd_WL_annotated_locus.bed")
Sb_WL_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Sb_WL_annotated_locus.bed")
Si_WL_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Si_WL_annotated_locus_new_v.bed")
Zm_WL_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Zm_WL_annotated_locus.bed")
Sb_BS_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Sb_BS_annotated_locus.bed")
Si_BS_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Si_BS_annotated_locus_new_v.bed")
Zm_BS_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_Zm_BS_annotated_locus.bed")
Sb_WL_DE_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_DE_Sb_WL_annotated_locus.bed")
Si_WL_DE_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_DE_Si_WL_annotated_locus_new_v.bed")
Zm_WL_DE_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_DE_Zm_WL_annotated_locus.bed")
Sb_BS_DE_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_DE_Sb_BS_annotated_locus.bed")
Si_BS_DE_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_DE_Si_BS_annotated_locus_new_v.bed")
Zm_BS_DE_A<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_DE_Zm_BS_annotated_locus.bed")

#Run Script
annotate_hits(DAP_Sb_WL,Sb_WL_A,"Sb_WL_DGF_annotated_DAP.bed","Sb_ID")
annotate_hits(DAP_Si_WL,Si_WL_A,"Si_WL_DGF_annotated_DAP.bed","Si_ID")
annotate_hits(DAP_Zm_WL,Zm_WL_A,"Zm_WL_DGF_annotated_DAP.bed","Zm_ID")
annotate_hits(DAP_Bd_WL,Bd_WL_A,"Bd_WL_DGF_annotated_DAP.bed","Bd_ID")
annotate_hits(DAP_Zm_BS,Zm_BS_A,"Zm_BS_DGF_annotated_DAP.bed","Zm_ID")
annotate_hits(DAP_Sb_BS,Sb_BS_A,"Sb_BS_DGF_annotated_DAP.bed","Sb_ID")
annotate_hits(DAP_Si_BS,Si_BS_A,"Si_BS_DGF_annotated_DAP.bed","Si_ID")
annotate_hits(IRL_Sb_WL,Sb_WL_A,"Sb_WL_DGF_annotated_IRL.bed","Sb_ID")
annotate_hits(IRL_Si_WL,Si_WL_A,"Si_WL_DGF_annotated_IRL.bed","Si_ID")
annotate_hits(IRL_Zm_WL,Zm_WL_A,"Zm_WL_DGF_annotated_IRL.bed","Zm_ID")
annotate_hits(IRL_Bd_WL,Bd_WL_A,"Bd_WL_DGF_annotated_IRL.bed","Bd_ID")
annotate_hits(IRL_Sb_BS,Sb_BS_A,"Sb_BS_DGF_annotated_IRL.bed","Sb_ID")
annotate_hits(IRL_Si_BS,Si_BS_A,"Si_BS_DGF_annotated_IRL.bed","Si_ID")
annotate_hits(IRL_Zm_BS,Zm_BS_A,"Zm_BS_DGF_annotated_IRL.bed","Zm_ID")
annotate_hits(DAP_Sb_WL_DE,Sb_WL_DE_A,"Sb_WL_DGF_DE_annotated_DAP.bed","Sb_ID")
annotate_hits(DAP_Si_WL_DE,Si_WL_DE_A,"Si_WL_DGF_DE_annotated_DAP.bed","Si_ID")
annotate_hits(DAP_Zm_WL_DE,Zm_WL_DE_A,"Zm_WL_DGF_DE_annotated_DAP.bed","Zm_ID")
annotate_hits(DAP_Zm_BS_DE,Zm_BS_DE_A,"Zm_BS_DGF_DE_annotated_DAP.bed","Zm_ID")
annotate_hits(DAP_Sb_BS_DE,Sb_BS_DE_A,"Sb_BS_DGF_DE_annotated_DAP.bed","Sb_ID")
annotate_hits(DAP_Si_BS_DE,Si_BS_DE_A,"Si_BS_DGF_DE_annotated_DAP.bed","Si_ID")
