#Rcode to create Figure 4G
#################################################

#Looking to see if motif instances are enriched in one cell-type compared to another
# all genes
runTest<-function(hitInSample, hitInPop, failInPop, sampleSize){
  
  #Test for over representation
  over<-phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  
  return(over)
  
}

processData<-function(x){
  
  n_hitInSample<-as.numeric(x[2])
  n_hitInPop<-as.numeric(x[4])
  n_failInPop<-as.numeric(x[5])-n_hitInPop

  #Run Test on BS enriched
  results<-runTest(n_hitInSample, n_hitInPop, n_failInPop, as.numeric(x[6]))
  return(results)
}

calculateFamilies<-function(in_list){
  #extract the relevant data files for footprints by family
  Sb_BS_DAP_DE_F<-in_list[[1]][,c("Chr","Start","Stop","Family","Locus")]
  Sb_DAP_minus_BS_DE_F<-in_list[[2]][,c("Chr","Start","Stop","Family","Locus")]
  Sb_WL_DAP_DE_F<-in_list[[3]][,c("Chr","Start","Stop","Family","Locus")]
  Sb_DAP_minus_WL_DE_F<-in_list[[4]][,c("Chr","Start","Stop","Family","Locus")]
  Si_BS_DAP_DE_F<-in_list[[5]][,c("Chr","Start","Stop","Family","Locus")]
  Si_DAP_minus_BS_DE_F<-in_list[[6]][,c("Chr","Start","Stop","Family","Locus")]
  Si_WL_DAP_DE_F<-in_list[[7]][,c("Chr","Start","Stop","Family","Locus")]
  Si_DAP_minus_WL_DE_F<-in_list[[8]][,c("Chr","Start","Stop","Family","Locus")]
  Zm_BS_DAP_DE_F<-in_list[[9]][,c("Chr","Start","Stop","Family","Locus")]
  Zm_DAP_minus_BS_DE_F<-in_list[[10]][,c("Chr","Start","Stop","Family","Locus")]
  Zm_WL_DAP_DE_F<-in_list[[11]][,c("Chr","Start","Stop","Family","Locus")]
  Zm_DAP_minus_WL_DE_F<-in_list[[12]][,c("Chr","Start","Stop","Family","Locus")]
  
  #List of differentially expressed genes for enrichment test 
  Sb_BS_DE<-in_list[[13]]
  Sb_M_DE<-in_list[[14]]
  Si_BS_DE<-in_list[[15]]
  Si_M_DE<-in_list[[16]]
  Zm_BS_DE<-in_list[[17]]
  Zm_M_DE<-in_list[[18]]
  
  #Calculate whether motifs are enriched in DE footprints compared to all remaining footprints
  Sb_BS_vs_all_F<-calculateEnrichment(Sb_BS_DAP_DE_F,Sb_DAP_minus_BS_DE_F)
  Sb_WL_vs_all_F<-calculateEnrichment(Sb_WL_DAP_DE_F,Sb_DAP_minus_WL_DE_F)
  Si_BS_vs_all_F<-calculateEnrichment(Si_BS_DAP_DE_F,Si_DAP_minus_BS_DE_F)
  Si_WL_vs_all_F<-calculateEnrichment(Si_WL_DAP_DE_F,Si_DAP_minus_WL_DE_F)
  Zm_BS_vs_all_F<-calculateEnrichment(Zm_BS_DAP_DE_F,Zm_DAP_minus_BS_DE_F)
  Zm_WL_vs_all_F<-calculateEnrichment(Zm_WL_DAP_DE_F,Zm_DAP_minus_WL_DE_F)
  
  #Calculate whether motifs are enriched in DE genes in DE foorprints compared to all
  Sb_DE_F<-Calc_DE_DGF_enrichment(Sb_BS_DE,Sb_M_DE,Sb_BS_DAP_DE_F,Sb_WL_DAP_DE_F,Sb_DAP_minus_BS_DE_F,Sb_DAP_minus_WL_DE_F)
  Si_DE_F<-Calc_DE_DGF_enrichment(Si_BS_DE,Si_M_DE,Si_BS_DAP_DE_F,Si_WL_DAP_DE_F,Si_DAP_minus_BS_DE_F,Si_DAP_minus_WL_DE_F)
  Zm_DE_F<-Calc_DE_DGF_enrichment(Zm_BS_DE,Zm_M_DE,Zm_BS_DAP_DE_F,Zm_WL_DAP_DE_F,Zm_DAP_minus_BS_DE_F,Zm_DAP_minus_WL_DE_F)
  
  return(list(Sb_BS_vs_all_F,
              Sb_WL_vs_all_F,
              Si_BS_vs_all_F,
              Si_WL_vs_all_F,
              Zm_BS_vs_all_F,
              Zm_WL_vs_all_F,
              Sb_DE_F,
              Si_DE_F,
              Zm_DE_F))
}

calculateEnrichment<-function(file_1,file_2){
  #only keep info on DGF TF interactions
  file_1<-file_1[,c(1,2,3,4)]
  file_2<-file_2[,c(1,2,3,4)]
  
  #remove duplicates
  file_1<-file_1[!duplicated(file_1),]
  file_2<-file_2[!duplicated(file_2),]
  file_1_DGFs<-unique(file_1[,c(1,2,3)])
  file_2_DGFs<-unique(file_2[,c(1,2,3)])
  
  #rename columns for merging
  colnames(file_1)<-c("Chr","Start","Stop","TF")
  colnames(file_2)<-c("Chr","Start","Stop","TF")
  
  #create combined set
  all<-rbind(file_1,file_2)
  all_DGFs<-unique(all[,c(1,2,3)])
  
  #Calculate the frequency of TF interactions in both datasets
  Freq_file_1_TF<-data.frame(table(file_1$TF))
  Freq_file_2_TF<-data.frame(table(file_2$TF))
  
  #rename columns for merging
  colnames(Freq_file_1_TF)<-c("Motif","file_1_Freq")
  colnames(Freq_file_2_TF)<-c("Motif","file_2_Freq")
  Total<-merge(Freq_file_1_TF,Freq_file_2_TF)
  
  #Calculate total number of hits in the population
  Total$Sum<-Total$file_1_Freq+Total$file_2_Freq
  
  #Caluclate population Size
  Total$PopSize<-length(all_DGFs[,1])
  
  #Caluclate population Size
  Total$sampleSize<-length(file_1_DGFs[,1])

  Total$pvalue<-apply(X=Total,MARGIN=1,FUN=processData)
  Total$adj.pvalue<-p.adjust(Total$pvalue,method="BH")
  Total<-Total[order(Total$adj.pvalue),]
  return(Total)
}


Calc_DE_DGF_enrichment<-function(BS_DE,M_DE,BS_DAP_DE,WL_DAP_DE,DAP_minus_BS_DE,DAP_minus_WL_DE){
    
  #Extract DE footprints in DE genes for comparison
  BS_DE_DGF_DE_genes_Positive<-subset(BS_DAP_DE,BS_DAP_DE[,5] %in% BS_DE)
  BS_DE_DGF_DE_genes_Negative<-subset(BS_DAP_DE,BS_DAP_DE[,5] %in% M_DE)
  WL_DE_DGF_DE_genes_Positive<-subset(WL_DAP_DE,WL_DAP_DE[,5] %in% M_DE)
  WL_DE_DGF_DE_genes_Negative<-subset(WL_DAP_DE,WL_DAP_DE[,5] %in% BS_DE)

  #Create sets for comparison
  #BS 
  #subset all DE DGFs not in DE genes
  temp1<-BS_DAP_DE[!(BS_DAP_DE[,5] %in% BS_DE_DGF_DE_genes_Positive[,5]),]
  #join to the remainder of DGFs
  BS_positive<-rbind(temp1,DAP_minus_BS_DE)
  #subset all DE DGFs not in DE genes
  temp2<-BS_DAP_DE[!(BS_DAP_DE[,5] %in% BS_DE_DGF_DE_genes_Negative[,5]),]
  #join to the remainder of DGFs
  BS_negative<-rbind(temp2,DAP_minus_BS_DE)

  #WL
  #subset all DE DGFs not in DE genes
  temp3<-WL_DAP_DE[!(WL_DAP_DE[,5] %in% WL_DE_DGF_DE_genes_Positive[,5]),]
  WL_positive<-rbind(temp3,DAP_minus_WL_DE)
  temp4<-WL_DAP_DE[!(WL_DAP_DE[,5] %in% WL_DE_DGF_DE_genes_Negative[,5]),]
  WL_negative<-rbind(temp4,DAP_minus_WL_DE)
  
  #Run Test 
  WL_DE_DGF_positive_vs_all<-calculateEnrichment(WL_DE_DGF_DE_genes_Positive,WL_positive)
  WL_DE_DGF_negative_vs_all<-calculateEnrichment(WL_DE_DGF_DE_genes_Negative,WL_negative)
  BS_DE_DGF_positive_vs_all<-calculateEnrichment(BS_DE_DGF_DE_genes_Positive,BS_positive)
  BS_DE_DGF_negative_vs_all<-calculateEnrichment(BS_DE_DGF_DE_genes_Negative,BS_negative)
  

  #Return Data
  res<-list(WL_DE_DGF_positive_vs_all,
            WL_DE_DGF_negative_vs_all,
            BS_DE_DGF_positive_vs_all,
            BS_DE_DGF_negative_vs_all)
  
  return(res)
}

#################################################
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/")

#Read in files
#DE DGF only
Sb_BS_DAP_DE<-read.table("Sb_BS_DGF_DE_annotated_DAP.bed")
Si_BS_DAP_DE<-read.table("Si_BS_DGF_DE_annotated_DAP.bed")
Zm_BS_DAP_DE<-read.table("Zm_BS_DGF_DE_annotated_DAP.bed")
Sb_WL_DAP_DE<-read.table("Sb_WL_DGF_DE_annotated_DAP.bed")
Si_WL_DAP_DE<-read.table("Si_WL_DGF_DE_annotated_DAP.bed")
Zm_WL_DAP_DE<-read.table("Zm_WL_DGF_DE_annotated_DAP.bed")
#non-DE DGF
Sb_BS_DAP<-read.table("Sb_BS_DGF_annotated_DAP.bed")
Si_BS_DAP<-read.table("Si_BS_DGF_annotated_DAP.bed")
Zm_BS_DAP<-read.table("Zm_BS_DGF_annotated_DAP.bed")
Sb_WL_DAP<-read.table("Sb_WL_DGF_annotated_DAP.bed")
Si_WL_DAP<-read.table("Si_WL_DGF_annotated_DAP.bed")
Zm_WL_DAP<-read.table("Zm_WL_DGF_annotated_DAP.bed")
Bd_WL_DAP<-read.table("Bd_WL_DGF_annotated_DAP.bed")
#All DGF minus DE
Sb_DAP_minus_BS_DE<-read.table("Sb_DGF_all_minus_BS_DE_annotated.bed")
Sb_DAP_minus_WL_DE<-read.table("Sb_DGF_all_minus_WL_DE_annotated.bed")
Si_DAP_minus_BS_DE<-read.table("Si_DGF_all_minus_BS_DE_annotated.bed")
Si_DAP_minus_WL_DE<-read.table("Si_DGF_all_minus_WL_DE_annotated.bed")
Zm_DAP_minus_BS_DE<-read.table("Zm_DGF_all_minus_BS_DE_annotated.bed")
Zm_DAP_minus_WL_DE<-read.table("Zm_DGF_all_minus_WL_DE_annotated.bed")
#All DGF
Sb_DAP<-read.table("Sb_DGF_all_annotated.bed")
Si_DAP<-read.table("Si_DGF_all_annotated.bed")
Zm_DAP<-read.table("Zm_DGF_all_annotated.bed")
#All DGF in C4 species
SbSiZm_DAP<-read.table("SbSiZm_DGF_all_annotated.bed")
ZmSb_DAP<-rbind(Sb_DAP,Zm_DAP)

#####################
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

#################################################
#Read In annotation Files
#Note - instances where phytozome mapping did not provide an ortholog were corrected by manual addition of the nearest matching sequence
Si_C4_genes<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_and_calvin.csv",head=F,sep=",")
Sb_C4_genes<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_Calvin_Zm_Sb.csv",head=T,sep="\t")
Bd_C4_genes<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_Calvin_Zm_Bd.csv",head=T,sep="\t")
Sb_C4_genes<-Sb_C4_genes[,c(1,3)]
Bd_C4_genes<-Bd_C4_genes[,c(1,3)]
#Create Final annotation dataframe
colnames(Si_C4_genes)<-c("Gene","Zm_ID","Si_ID")
colnames(Sb_C4_genes)<-c("Zm_ID","Sb_ID")
colnames(Bd_C4_genes)<-c("Zm_ID","Bd_ID")
C4<-merge(Si_C4_genes,Sb_C4_genes,by="Zm_ID")
C4<-merge(C4,Bd_C4_genes,by="Zm_ID")

##################################################
#Look at enrichment of Motif Family Results
##################################################

#Calculate families
DAP_annotation<-read.table("/home/sjb287/documents/Motif_Searching/DAP_Motifs/DAP_motif_annotation.txt")
DAP_annotation<-DAP_annotation[,c(1,3)]
colnames(DAP_annotation)<-c("Motif","Family")

#Use the same input files as for individual motifs
Sb_BS_DAP_DE_Families<-Sb_BS_DAP_DE
Sb_WL_DAP_DE_Families<-Sb_WL_DAP_DE
Si_BS_DAP_DE_Families<-Si_BS_DAP_DE
Si_WL_DAP_DE_Families<-Si_WL_DAP_DE
Zm_BS_DAP_DE_Families<-Zm_BS_DAP_DE
Zm_WL_DAP_DE_Families<-Zm_WL_DAP_DE
Bd_WL_DAP_Families<-Bd_WL_DAP
Si_WL_DAP_Families<-Si_WL_DAP
Sb_WL_DAP_Families<-Sb_WL_DAP
Zm_WL_DAP_Families<-Zm_WL_DAP
Sb_DAP_minus_BS_DE_Families<-Sb_DAP_minus_BS_DE
Sb_DAP_minus_WL_DE_Families<-Sb_DAP_minus_WL_DE
Si_DAP_minus_BS_DE_Families<-Si_DAP_minus_BS_DE
Si_DAP_minus_WL_DE_Families<-Si_DAP_minus_WL_DE
Zm_DAP_minus_BS_DE_Families<-Zm_DAP_minus_BS_DE
Zm_DAP_minus_WL_DE_Families<-Zm_DAP_minus_WL_DE

#Rename Columns for merging
colnames(Sb_BS_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_WL_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_BS_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_WL_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_BS_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_WL_DAP_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Bd_WL_DAP_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_DAP_minus_BS_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Sb_DAP_minus_WL_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_DAP_minus_BS_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Si_DAP_minus_WL_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_DAP_minus_BS_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")
colnames(Zm_DAP_minus_WL_DE_Families)<-c("Chr","Start","Stop","Motif","Locus")

#Merge individual motifs with TF family annotations
Sb_BS_DAP_DE_Families<-unique(merge(Sb_BS_DAP_DE_Families,DAP_annotation,by="Motif"))
Sb_WL_DAP_DE_Families<-unique(merge(Sb_WL_DAP_DE_Families,DAP_annotation,by="Motif"))
Si_BS_DAP_DE_Families<-unique(merge(Si_BS_DAP_DE_Families,DAP_annotation,by="Motif"))
Si_WL_DAP_DE_Families<-unique(merge(Si_WL_DAP_DE_Families,DAP_annotation,by="Motif"))
Zm_BS_DAP_DE_Families<-unique(merge(Zm_BS_DAP_DE_Families,DAP_annotation,by="Motif"))
Zm_WL_DAP_DE_Families<-unique(merge(Zm_WL_DAP_DE_Families,DAP_annotation,by="Motif"))
Sb_DAP_minus_BS_DE_Families<-unique(merge(Sb_DAP_minus_BS_DE_Families,DAP_annotation,by="Motif"))
Sb_DAP_minus_WL_DE_Families<-unique(merge(Sb_DAP_minus_WL_DE_Families,DAP_annotation,by="Motif"))
Si_DAP_minus_BS_DE_Families<-unique(merge(Si_DAP_minus_BS_DE_Families,DAP_annotation,by="Motif"))
Si_DAP_minus_WL_DE_Families<-unique(merge(Si_DAP_minus_WL_DE_Families,DAP_annotation,by="Motif"))
Zm_DAP_minus_BS_DE_Families<-unique(merge(Zm_DAP_minus_BS_DE_Families,DAP_annotation,by="Motif"))
Zm_DAP_minus_WL_DE_Families<-unique(merge(Zm_DAP_minus_WL_DE_Families,DAP_annotation,by="Motif"))
Sb_WL_DAP_Families<-unique(merge(Sb_WL_DAP_Families,DAP_annotation,by="Motif"))
Si_WL_DAP_Families<-unique(merge(Si_WL_DAP_Families,DAP_annotation,by="Motif"))
Zm_WL_DAP_Families<-unique(merge(Zm_WL_DAP_Families,DAP_annotation,by="Motif"))
Bd_WL_DAP_Families<-unique(merge(Bd_WL_DAP_Families,DAP_annotation,by="Motif"))

#Perform enrichment tests 
families_res<-calculateFamilies(list(Sb_BS_DAP_DE_Families,
                                Sb_DAP_minus_BS_DE_Families,
                                Sb_WL_DAP_DE_Families,
                                Sb_DAP_minus_WL_DE_Families,
                                Si_BS_DAP_DE_Families,
                                Si_DAP_minus_BS_DE_Families,
                                Si_WL_DAP_DE_Families,
                                Si_DAP_minus_WL_DE_Families,
                                Zm_BS_DAP_DE_Families,
                                Zm_DAP_minus_BS_DE_Families,
                                Zm_WL_DAP_DE_Families,
                                Zm_DAP_minus_WL_DE_Families,
                                Sb_BS_DE,
                                Sb_M_DE,
                                Si_BS_DE,
                                Si_M_DE,
                                Zm_BS_DE,
                                Zm_M_DE))
#Results generated:
#In DE footprints compared to all footprints
Sb_BS_vs_all_F<-families_res[[1]]
Sb_WL_vs_all_F<-families_res[[2]]
Si_BS_vs_all_F<-families_res[[3]]
Si_WL_vs_all_F<-families_res[[4]]
Zm_BS_vs_all_F<-families_res[[5]]
Zm_WL_vs_all_F<-families_res[[6]]
#In DE genes in DE footprints compared to all
Sb_BS_F<-families_res[[7]]
Si_BS_F<-families_res[[8]]
Zm_BS_F<-families_res[[9]]
#Subset to identify motif families over represented in DE footprints
#this data is what goes into making figure 4G
Sb_BS_vs_all_F<-subset(Sb_BS_vs_all_F,Sb_BS_vs_all_F$adj.pvalue<0.05)
Sb_WL_vs_all_F<-subset(Sb_WL_vs_all_F,Sb_WL_vs_all_F$adj.pvalue<0.05)
Si_BS_vs_all_F<-subset(Si_BS_vs_all_F,Si_BS_vs_all_F$adj.pvalue<0.05)
Si_WL_vs_all_F<-subset(Si_WL_vs_all_F,Si_WL_vs_all_F$adj.pvalue<0.05)
Zm_BS_vs_all_F<-subset(Zm_BS_vs_all_F,Zm_BS_vs_all_F$adj.pvalue<0.05)
Zm_WL_vs_all_F<-subset(Zm_WL_vs_all_F,Zm_WL_vs_all_F$adj.pvalue<0.05)

Bd_WL_DAP_F<-Bd_WL_DAP_Families[,c("Chr","Start","Stop","Family","Locus")]
Sb_WL_DAP_F<-Sb_WL_DAP_Families[,c("Chr","Start","Stop","Family","Locus")]
Si_WL_DAP_F<-Si_WL_DAP_Families[,c("Chr","Start","Stop","Family","Locus")]
Zm_WL_DAP_F<-Zm_WL_DAP_Families[,c("Chr","Start","Stop","Family","Locus")]
C4_WL_F<-rbind(Sb_WL_DAP_F,Si_WL_DAP_F,Zm_WL_DAP_F)

Bd_WL_vs_C4_F<-calculateEnrichment(Bd_WL_DAP_F,C4_WL_F)
C4_vs_Bd_WL_F<-calculateEnrichment(C4_WL_F,Bd_WL_DAP_F)
C4_vs_Bd_WL<-calculateEnrichment(SbSiZm_DAP,Bd_WL_DAP)
C4_vs_Bd_WL<-subset(C4_vs_Bd_WL,C4_vs_Bd_WL$adj.pvalue<0.05)

C4_vs_Bd_WL_F<-subset(C4_vs_Bd_WL_F,C4_vs_Bd_WL_F$adj.pvalue<0.05)
Bd_WL_vs_C4_F<-subset(Bd_WL_vs_C4_F,Bd_WL_vs_C4_F$adj.pvalue<0.05)
