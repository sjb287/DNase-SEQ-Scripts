#Code required to make Figure 4D 
#Figure 4D
#################################################
#Calculate the number of hits per family
FamilyHits<-function(infile,ID){
  
  #Read in hits file and rename columns for merging
  temp1<-infile
  temp1<-temp1[!duplicated(temp1),]
  colnames(temp1)<-c("Chr","Start","Stop","TF","Locus")
  
  #annotate the hits with TF family
  temp2<-merge(temp1,DAP,by="TF")
  
  #select only columns of interest
  temp2<-temp2[,c("Chr","Start","Stop","Family","Locus")]
  
  #the columns are than renamed again, this is done so the return vector can be used
  #for calculations using functions written for TFs
  colnames(temp2)<-c("Chr","start","stop","TF",ID)
  
  #return results
  temp2<-temp2[!duplicated(temp2),]
  return(temp2)
}

#################################################
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP")
#Define files for TF family annotations
DAP_annotation="Data_Figure_4_motif_annotation.txt"
DAP<-read.table(DAP_annotation,head=F)
colnames(DAP)<-c("TF","AGI","Family")

#read in hits
At_WL_hits<-read.table("/data/reads/DNase_data/2012_DNase_Zhang/wdleaf_FDR/At_WL_DGF_fdr_0.01_fimo_hits_all.bed", head=F)
colnames(At_WL_hits)<-c("Chr","Start","Stop","TF")
At_WL_hits<-At_WL_hits[!duplicated(At_WL_hits),]

#Annotate hits
At_WL_Locus<-read.table("/data/reads/DNase_data/2012_DNase_Zhang/wdleaf_FDR/temp_2.bed", head=F)
colnames(At_WL_Locus)<-c("Chr","Start","Stop","Locus")
At_WL_Locus<-At_WL_Locus[!duplicated(At_WL_Locus),]

#Remove duplicate entries
At_WL<-merge(At_WL_hits,At_WL_Locus)
At_WL<-At_WL[!duplicated(At_WL),]

#Get Hits by Family
At_Family<-FamilyHits(At_WL,"At_ID")

#Calculate Frequencies
Freq_At_WL_TF<-data.frame(table(At_WL$TF))
#Calculate Frequencies Family
Freq_At_WL_Families<-data.frame(table(At_Family$TF))

#Order
Freq_At_WL_TF<-Freq_At_WL_TF[order(Freq_At_WL_TF$Freq),]
Freq_At_WL_Families<-Freq_At_WL_Families[order(Freq_At_WL_Families$Freq),]

#Calculate Ranks
Freq_At_WL_TF$Rank_At<-rank(-Freq_At_WL_TF$Freq)
Freq_At_WL_Families$Rank_At<-rank(-Freq_At_WL_Families$Freq)

#Specify outfile name
outfile<-"Data_Figure_4_At_WL_known.txt"
outfile2<-"Data_Figure_4_At_WL_Families.txt"

#Write resuts to file
write.table(Freq_At_WL_TF,outfile,quote=F,row.names=F)  
write.table(Freq_At_WL_Families,outfile2,quote=F,row.names=F)  