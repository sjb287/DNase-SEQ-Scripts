#Code to rank motif instances for Figure 4

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

  #define outfile names
  outfile_1<-paste("Summary_Frequency_WL_Motif",suffix,sep="_")
  outfile_2<-paste("Summary_Frequency_BS_Motif",suffix,sep="_")
  outfile_3<-paste("Summary_Frequency_3_way_WL_Motif",suffix,sep="_")
  outfile_4<-paste("Summary_Frequency_3_way_BS_Motif",suffix,sep="_")
  outfile_5<-"Data_Figure_4_Bd_WL.txt"
  outfile_6<-"Data_Figure_4_Sb_WL.txt"
  outfile_7<-"Data_Figure_4_Si_WL.txt"
  outfile_8<-"Data_Figure_4_Zm_WL.txt"
  outfile_9<-"Data_Figure_4_Sb_BS.txt"
  outfile_10<-"Data_Figure_4_Si_BS.txt"
  outfile_11<-"Data_Figure_4_Zm_BS.txt"
  
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

runProgram<-function(Bd_WL,
                     Sb_BS,
                     Sb_WL,
                     Si_BS,
                     Si_WL,
                     Zm_BS,
                     Zm_WL,suffix){
  
  #Read in motif files
  Bd_WL<-read.table(Bd_WL,head=T)
  Sb_BS<-read.table(Sb_BS,head=T)
  Sb_WL<-read.table(Sb_WL,head=T)
  Si_BS<-read.table(Si_BS,head=T)
  Si_WL<-read.table(Si_WL,head=T)
  Zm_BS<-read.table(Zm_BS,head=T)
  Zm_WL<-read.table(Zm_WL,head=T)

  CreateTable(Bd_WL,
              Sb_BS,
              Sb_WL,
              Si_BS,
              Si_WL,
              Zm_BS,
              Zm_WL,suffix
  )
}


#Define files for TF family annotations
Motif_annotation="Data_Figure_4_motif_annotation.txt"
Motifs<-read.table(Motif_annotation,head=F)
colnames(Motifs)<-c("TF","AGI","Family")

#known motifs
Bd_WL_motifs<-"DGF_known_annotated_Bd_WL.bed"
Sb_BS_motifs<-"DGF_known_annotated_Sb_BS.bed"
Sb_WL_motifs<-"DGF_known_annotated_Sb_WL.bed"
Si_BS_motifs<-"DGF_known_annotated_Si_BS.bed"
Si_WL_motifs<-"DGF_known_annotated_Si_WL.bed"
Zm_BS_motifs<-"DGF_known_annotated_Zm_BS.bed"
Zm_WL_motifs<-"DGF_known_annotated_Zm_WL.bed"

runProgram(Bd_WL_motifs,
            Sb_BS_motifs,
            Sb_WL_motifs,
            Si_BS_motifs,
            Si_WL_motifs,
            Zm_BS_motifs,
            Zm_WL_motifs,"known.txt"
)