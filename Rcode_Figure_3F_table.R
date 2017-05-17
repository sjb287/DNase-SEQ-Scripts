extractSubset<-function(input,ID,outfile){
  
  all_data<-read.table(input,head=F)
  subset_data<-subset(all_data,all_data[,4] %in% annotations[,c(ID)])
  subset_data<-subset_data[!duplicated(subset_data),]
  return(subset_data)  
}



setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/")
#Files with all DGF annotated
DGF_Sb_WL_A="DGF_Sb_WL_annotated_locus.bed"
DGF_Sb_BS_A="DGF_Sb_BS_annotated_locus.bed"
DGF_Zm_WL_A="DGF_Zm_WL_annotated_locus.bed"
DGF_Zm_BS_A="DGF_Zm_BS_annotated_locus.bed"
#Files with crossmapped DGFs annotated
DGF_BS_Sb_Zm_occupied="DGF_BS_S.bicolor_to_Z.mays_conserved_and_occupied_annotated.bed"
DGF_WL_Sb_Zm_occupied="DGF_WL_S.bicolor_to_Z.mays_conserved_and_occupied_annotated.bed"
DGF_BS_Zm_Sb_occupied="DGF_BS_Z.mays_to_S.bicolor_conserved_and_occupied_annotated.bed"
DGF_WL_Zm_Sb_occupied="DGF_WL_Z.mays_to_S.bicolor_conserved_and_occupied_annotated.bed"
DGF_BS_Sb_Zm_conserved="DGF_BS_S.bicolor_to_Z.mays_conserved_annotated.bed"
DGF_WL_Sb_Zm_conserved="DGF_WL_S.bicolor_to_Z.mays_conserved_annotated.bed"
DGF_BS_Zm_Sb_conserved="DGF_BS_Z.mays_to_S.bicolor_conserved_annotated.bed"
DGF_WL_Zm_Sb_conserved="DGF_WL_Z.mays_to_S.bicolor_conserved_annotated.bed"

#read in annotations
annotations<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_annotations/C4_Calvin_Zm_Sb.csv",head=T)
annotations<-annotations[,c(1,3)]

#rename columns for processing
colnames(annotations)<-c("Zm_ID","Sb_ID")

#Get subset of all C4 and Calvin Cycle DGFs
C4_Sb_WL<-extractSubset(DGF_Sb_WL_A,"Sb_ID")
C4_Sb_BS<-extractSubset(DGF_Sb_BS_A,"Sb_ID")
C4_Zm_WL<-extractSubset(DGF_Zm_WL_A,"Zm_ID")
C4_Zm_BS<-extractSubset(DGF_Zm_BS_A,"Zm_ID")

#Get subset of all C4 and Calvin Cycle DGFs that are conserved
C4_Sb_WL_conserved<-extractSubset(DGF_WL_Sb_Zm_conserved,"Sb_ID")
C4_Sb_BS_conserved<-extractSubset(DGF_BS_Sb_Zm_conserved,"Sb_ID")
C4_Zm_WL_conserved<-extractSubset(DGF_WL_Zm_Sb_conserved,"Zm_ID")
C4_Zm_BS_conserved<-extractSubset(DGF_WL_Zm_Sb_conserved,"Zm_ID")

#Get subset of all C4 and Calvin Cycle DGFs that are conserved and occupied
C4_Sb_WL_occupied<-extractSubset(DGF_WL_Sb_Zm_occupied,"Sb_ID")
C4_Sb_BS_occupied<-extractSubset(DGF_BS_Sb_Zm_occupied,"Sb_ID")
C4_Zm_WL_occupied<-extractSubset(DGF_WL_Zm_Sb_occupied,"Zm_ID")
C4_Zm_BS_occupied<-extractSubset(DGF_BS_Zm_Sb_occupied,"Zm_ID")

#create a results table, using the length function to count the number of DGF
results<-data.frame("Species"=c("S. bicolor",
                                "S. bicolor",
                                "Z. mays",
                                "Z. mays",
                                "S. bicolor",
                                "Z. mays"),
                    "Tissue"=c("WL",
                                "BS",
                                "WL",
                                "BS",
                               "Total",
                               "Total"),
                    "Num_DGF"=c(length(C4_Sb_WL[,1]),
                                length(C4_Sb_BS[,1]),
                                length(C4_Zm_WL[,1]),
                                length(C4_Zm_BS[,1]),
                                length(C4_Sb_WL[,1])+length(C4_Sb_BS[,1]),
                                length(C4_Zm_WL[,1])+length(C4_Zm_BS[,1])),
                    "Num_DGF_conserved"=c(length(C4_Sb_WL_conserved[,1]),
                                          length(C4_Sb_BS_conserved[,1]),
                                          length(C4_Zm_WL_conserved[,1]),
                                          length(C4_Zm_BS_conserved[,1]),
                                          length(C4_Sb_WL_conserved[,1])+length(C4_Sb_BS_conserved[,1]),
                                          length(C4_Zm_WL_conserved[,1])+length(C4_Zm_BS_conserved[,1])),
                    "Num_DGF_diverged"=c(length(C4_Sb_WL[,1])-length(C4_Sb_WL_conserved[,1]),
                                         length(C4_Sb_BS[,1])-length(C4_Sb_BS_conserved[,1]),
                                         length(C4_Zm_WL[,1])-length(C4_Zm_WL_conserved[,1]),
                                         length(C4_Zm_BS[,1])-length(C4_Zm_BS_conserved[,1]),
                                         length(C4_Sb_WL[,1])-length(C4_Sb_WL_conserved[,1])+length(C4_Sb_BS[,1])-length(C4_Sb_BS_conserved[,1]),
                                         length(C4_Zm_WL[,1])-length(C4_Zm_WL_conserved[,1])+length(C4_Zm_BS[,1])-length(C4_Zm_BS_conserved[,1])),
                    "Num_DGF_conserved_and_occupued"=c(length(C4_Sb_WL_occupied[,1]),
                                                       length(C4_Sb_BS_occupied[,1]),
                                                       length(C4_Zm_WL_occupied[,1]),
                                                       length(C4_Zm_BS_occupied[,1]),
                                                       length(C4_Sb_WL_occupied[,1])+length(C4_Sb_BS_occupied[,1]),
                                                       length(C4_Zm_WL_occupied[,1])+length(C4_Zm_BS_occupied[,1])))
