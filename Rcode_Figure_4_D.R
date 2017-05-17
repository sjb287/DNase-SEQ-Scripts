#code to create data for rank plots
#################################################
#Functions
#################################################
CreateTable<-function(Bd_WL,Sb_BS,Sb_WL,Si_BS,Si_WL,Zm_BS,Zm_WL,suffix){

  #Remove duplicate entries
  Bd_WL<-Bd_WL[!duplicated(Bd_WL),]
  Sb_WL<-Sb_WL[!duplicated(Sb_WL),]
  Sb_BS<-Sb_BS[!duplicated(Sb_BS),]
  Si_WL<-Si_WL[!duplicated(Si_WL),]
  Si_BS<-Si_BS[!duplicated(Si_BS),]
  Zm_WL<-Zm_WL[!duplicated(Zm_WL),]
  Zm_BS<-Zm_BS[!duplicated(Zm_BS),]

  #Calculate Frequencies
  Freq_Bd_WL_TF<-data.frame(table(Bd_WL$TF))
  Freq_Si_BS_TF<-data.frame(table(Si_BS$TF))
  Freq_Si_WL_TF<-data.frame(table(Si_WL$TF))
  Freq_Sb_BS_TF<-data.frame(table(Sb_BS$TF))
  Freq_Sb_WL_TF<-data.frame(table(Sb_WL$TF))
  Freq_Zm_BS_TF<-data.frame(table(Zm_BS$TF))
  Freq_Zm_WL_TF<-data.frame(table(Zm_WL$TF))

  #Order
  Freq_Bd_WL_TF<-Freq_Bd_WL_TF[order(Freq_Bd_WL_TF$Freq),]
  Freq_Si_WL_TF<-Freq_Si_WL_TF[order(Freq_Si_WL_TF$Freq),]
  Freq_Si_BS_TF<-Freq_Si_BS_TF[order(Freq_Si_BS_TF$Freq),]
  Freq_Sb_WL_TF<-Freq_Sb_WL_TF[order(Freq_Sb_WL_TF$Freq),]
  Freq_Sb_BS_TF<-Freq_Sb_BS_TF[order(Freq_Sb_BS_TF$Freq),]
  Freq_Zm_WL_TF<-Freq_Zm_WL_TF[order(Freq_Zm_WL_TF$Freq),]
  Freq_Zm_BS_TF<-Freq_Zm_BS_TF[order(Freq_Zm_BS_TF$Freq),]

  #Calculate Ranks
  Freq_Bd_WL_TF$Rank_Bd<-rank(-Freq_Bd_WL_TF$Freq)
  Freq_Si_WL_TF$Rank_Si<-rank(-Freq_Si_WL_TF$Freq)
  Freq_Si_BS_TF$Rank_Si<-rank(-Freq_Si_BS_TF$Freq)
  Freq_Sb_WL_TF$Rank_Sb<-rank(-Freq_Sb_WL_TF$Freq)
  Freq_Sb_BS_TF$Rank_Sb<-rank(-Freq_Sb_BS_TF$Freq)
  Freq_Zm_WL_TF$Rank_Zm<-rank(-Freq_Zm_WL_TF$Freq)
  Freq_Zm_BS_TF$Rank_Zm<-rank(-Freq_Zm_BS_TF$Freq)

  #Combine Data
  FreqSbZm_WL<-merge(Freq_Sb_WL_TF,Freq_Zm_WL_TF,by="Var1")
  colnames(FreqSbZm_WL)<-c("Var1","Freq_Sb","Rank_Sb","Freq_Zm","Rank_Zm")
  FreqSbZmSi_WL<-merge(FreqSbZm_WL,Freq_Si_WL_TF,by="Var1")
  colnames(FreqSbZmSi_WL)<-c("Var1","Freq_Sb","Rank_Sb","Freq_Zm","Rank_Zm","Freq_Si","Rank_Si")
  FreqSbZmSiBd_WL<-merge(FreqSbZmSi_WL,Freq_Bd_WL_TF,by="Var1")
  colnames(FreqSbZmSiBd_WL)<-c("TF",
                               "Freq_Sb",
                               "Rank_Sb",
                               "Freq_Zm",
                               "Rank_Zm",
                               "Freq_Si",
                               "Rank_Si",
                               "Freq_Bd",
                               "Rank_Bd"
                               )
  FreqSbZm_BS<-merge(Freq_Sb_BS_TF,Freq_Zm_BS_TF,by="Var1")
  colnames(FreqSbZm_BS)<-c("Var1","Freq_Sb","Rank_Sb","Freq_Zm","Rank_Zm")
  FreqSbZmSi_BS<-merge(FreqSbZm_BS,Freq_Si_BS_TF,by="Var1")
  colnames(FreqSbZmSi_BS)<-c("Var1","Freq_Sb","Rank_Sb","Freq_Zm","Rank_Zm","Freq_Si","Rank_Si")
  FreqSbZmSiBd_BS<-merge(FreqSbZmSi_BS,Freq_Bd_WL_TF,by="Var1")
  colnames(FreqSbZmSiBd_BS)<-c("TF",
                             "Freq_Sb",
                             "Rank_Sb",
                             "Freq_Zm",
                             "Rank_Zm",
                             "Freq_Si",
                             "Rank_Si",
                             "Freq_Bd",
                             "Rank_Bd"
                            )

  outfile_1<-paste("Summary_Frequency_WL_Motif",suffix,sep="_")
  outfile_2<-paste("Summary_Frequency_BS_Motif",suffix,sep="_")
  outfile_3<-paste("Summary_Frequency_3_way_WL_Motif",suffix,sep="_")
  outfile_4<-paste("Summary_Frequency_3_way_BS_Motif",suffix,sep="_")
  outfile_5<-paste("Data_Figure_4_Bd_WL",suffix,sep="_")
  outfile_6<-paste("Data_Figure_4_Sb_WL",suffix,sep="_")
  outfile_7<-paste("Data_Figure_4_Si_WL",suffix,sep="_")
  outfile_8<-paste("Data_Figure_4_Zm_WL",suffix,sep="_")
  outfile_9<-paste("Data_Figure_4_Sb_BS",suffix,sep="_")
  outfile_10<-paste("Data_Figure_4_Si_BS",suffix,sep="_")
  outfile_11<-paste("Data_Figure_4_Zm_BS",suffix,sep="_")
  
  write.table(FreqSbZmSiBd_WL,outfile_1,quote=F,row.names=F)
  write.table(FreqSbZmSiBd_BS,outfile_2,quote=F,row.names=F)
  write.table(FreqSbZmSi_WL,outfile_3,quote=F,row.names=F)
  write.table(FreqSbZmSi_BS,outfile_4,quote=F,row.names=F)
  write.table(Freq_Bd_WL_TF,outfile_5,quote=F,row.names=F)  
  write.table(Freq_Sb_WL_TF,outfile_6,quote=F,row.names=F)  
  write.table(Freq_Si_WL_TF,outfile_7,quote=F,row.names=F)  
  write.table(Freq_Zm_WL_TF,outfile_8,quote=F,row.names=F)  
  write.table(Freq_Sb_BS_TF,outfile_9,quote=F,row.names=F)  
  write.table(Freq_Si_BS_TF,outfile_10,quote=F,row.names=F)  
  write.table(Freq_Zm_BS_TF,outfile_11,quote=F,row.names=F)  
  
}

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

runProgram<-function(Bd_WL,
                     Sb_BS,
                     Sb_WL,
                     Si_BS,
                     Si_WL,
                     Zm_BS,
                     Zm_WL,suffix){
  
  #Read in motif files
  Si_WL<-read.table(Si_WL,head=F)
  Si_BS<-read.table(Si_BS,head=F)
  Sb_WL<-read.table(Sb_WL,head=F)
  Sb_BS<-read.table(Sb_BS,head=F)
  Zm_WL<-read.table(Zm_WL,head=F)
  Zm_BS<-read.table(Zm_BS,head=F)
  Bd_WL<-read.table(Bd_WL,head=F)
  colnames(Si_WL)<-c("Chr","start","stop","TF","Si_ID")
  colnames(Si_BS)<-c("Chr","start","stop","TF","Si_ID")
  colnames(Sb_WL)<-c("Chr","start","stop","TF","Sb_ID")
  colnames(Sb_BS)<-c("Chr","start","stop","TF","Sb_ID")
  colnames(Zm_WL)<-c("Chr","start","stop","TF","Zm_ID")
  colnames(Zm_BS)<-c("Chr","start","stop","TF","Zm_ID")
  colnames(Bd_WL)<-c("Chr","start","stop","TF","Bd_ID")
  
  CreateTable(Bd_WL,
              Sb_BS,
              Sb_WL,
              Si_BS,
              Si_WL,
              Zm_BS,
              Zm_WL,suffix
  )
}

runProgramTFfamilies<-function(Si_WL,Si_BS,Sb_WL,Sb_BS,Zm_WL,Zm_BS,Bd_WL,suffix){
  
  #read in files
  Si_WL<-read.table(Si_WL,head=F)
  Si_BS<-read.table(Si_BS,head=F)
  Sb_WL<-read.table(Sb_WL,head=F)
  Sb_BS<-read.table(Sb_BS,head=F)
  Zm_WL<-read.table(Zm_WL,head=F)
  Zm_BS<-read.table(Zm_BS,head=F)
  Bd_WL<-read.table(Bd_WL,head=F)
  
  #Count the number of hits by TF family rather than individual IDs
  Si_WL_families<-FamilyHits(Si_WL,"Si_ID")
  Si_BS_families<-FamilyHits(Si_BS,"Si_ID")
  Sb_WL_families<-FamilyHits(Sb_WL,"Sb_ID")
  Sb_BS_families<-FamilyHits(Sb_BS,"Sb_ID")
  Zm_WL_families<-FamilyHits(Zm_WL,"Zm_ID")
  Zm_BS_families<-FamilyHits(Zm_BS,"Zm_ID")
  Bd_WL_families<-FamilyHits(Bd_WL,"Bd_ID")
  
  #Run the program with Family hits
  CreateTable(Sb_WL_families,
             Sb_BS_families,
             Si_WL_families,
             Si_BS_families,
             Zm_WL_families,
             Zm_BS_families,
             Bd_WL_families,
             suffix
  )
}
#################################################
#Define files for TF family annotations
DAP_annotation="Data_Figure_4_motif_annotation.txt"
DAP<-read.table(DAP_annotation,head=F)
colnames(DAP)<-c("TF","AGI","Family")

setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP")
#Files containing DGFs annotated with known motifs
Bd_WL_DAP<-"Bd_WL_DGF_annotated_DAP.bed"
Sb_BS_DAP<-"Sb_BS_DGF_annotated_DAP.bed"
Sb_WL_DAP<-"Sb_WL_DGF_annotated_DAP.bed"
Si_BS_DAP<-"Si_BS_DGF_annotated_DAP.bed"
Si_WL_DAP<-"Si_WL_DGF_annotated_DAP.bed"
Zm_BS_DAP<-"Zm_BS_DGF_annotated_DAP.bed"
Zm_WL_DAP<-"Zm_WL_DGF_annotated_DAP.bed"

#DE
Si_BS_DE_DAP="Si_BS_DGF_DE_annotated_DAP.bed"
Sb_BS_DE_DAP="Sb_BS_DGF_DE_annotated_DAP.bed"
Zm_BS_DE_DAP="Zm_BS_DGF_DE_annotated_DAP.bed"
Si_WL_DE_DAP="Si_WL_DGF_DE_annotated_DAP.bed"
Sb_WL_DE_DAP="Sb_WL_DGF_DE_annotated_DAP.bed"
Zm_WL_DE_DAP="Zm_WL_DGF_DE_annotated_DAP.bed"

runProgramTFfamilies(Si_WL_DAP,
                     Si_BS_DAP,
                     Sb_WL_DAP,
                     Sb_BS_DAP,
                     Zm_WL_DAP,
                     Zm_BS_DAP,
                     Bd_WL_DAP,
                     "known_Families.txt"
)

runProgramTFfamilies(Si_WL_DE_DAP,
                     Si_BS_DE_DAP,
                     Sb_WL_DE_DAP,
                     Sb_BS_DE_DAP,
                     Zm_WL_DE_DAP,
                     Zm_BS_DE_DAP,
                     Bd_WL_DAP,
                     "DAP_DE_families.txt"
)

runProgram(Bd_WL_DAP,
            Sb_BS_DAP,
            Sb_WL_DAP,
            Si_BS_DAP,
            Si_WL_DAP,
            Zm_BS_DAP,
            Zm_WL_DAP,"known.txt"
)

runProgram(Bd_WL_DAP,
           Sb_BS_DE_DAP,
           Sb_WL_DE_DAP,
           Si_BS_DE_DAP,
           Si_WL_DE_DAP,
           Zm_BS_DE_DAP,
           Zm_WL_DE_DAP,"DAP_DE.txt"
)