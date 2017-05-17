#R code to make table for use in figure 3C

#################################################
#Function to Calculate the mean number of mapped and mapped and conserved DGFs for each randomization trial
calcMeans<-function(df,Sp,Tissue){
  
  #calculate the mean values from 100 random trials  
  df$mean<-apply(df[,c(3:102)],1,function(x) mean(x))
  #Add a column specifying the species DGF data that was randomized
  df$Species<-Sp
  #Add a colum specifying the tissue 
  df$Tissue<-Tissue
  #Split the summary table into two df, one for conserved and one for conserved and occupied
  temp_1<-split(df,df$V2)
  #create a results dataframe
  temp_2<-data.frame("Species_1"=temp_1[[1]]$Species,
                     "Species_2"=temp_1[[1]][,1],
                     "Tissue"=temp_1[[1]]$Tissue,
                     "Comparison"=c(rep("Positional Conservation",length(temp_1[[2]][,1]))),
                     "Random_average_conserved"=temp_1[[2]]$mean,
                     "Random_average_conserved_and_occupied"=temp_1[[1]]$mean)
  temp_2$Samples<-paste(temp_2$Species_1,temp_2$Species_2,temp_2$Tissue,sep="_")
  return(temp_2)
  
}

readData<-function(x){
  #open file
  temp<-read.table(x,head=F)
  #calculate the number of footprints
  result<-length(unique(temp[,4]))
  #return the result
  return(result)
  
}


#################################################
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/Random/")
#File locations for the random cross map comparisons
R_Sb_WL=read.table("Random_Summary_Sb_WL.txt",head=F)
R_Sb_BS=read.table("Random_Summary_Sb_BS.txt",head=F)

setwd("/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/")
#Files with all DGFs
DGF_Sb_WL=read.table("DGF_Sb_WL.bed")
DGF_Sb_BS=read.table("DGF_Sb_BS.bed")

#Calculate the number of DGFs
Num_Sb_WL<-length(DGF_Sb_WL[,1])
Num_Sb_BS<-length(DGF_Sb_BS[,1])
Num_DGF<-data.frame("Tissue"=c("WL","BS"),"Num_DGF"=c(Num_Sb_WL,Num_Sb_BS))

setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/")

#Files with crossmapped DGFs (conserved)
conserved<-data.frame("Conserved"=
            c("DGF_BS_S.bicolor_to_Z.mays_conserved.bed",
             "DGF_BS_S.bicolor_to_S.italica_conserved.bed",
             "DGF_BS_S.bicolor_to_B.distachyon_conserved.bed",
             "DGF_S.bicolor_to_Z.mays_conserved.bed",
             "DGF_S.bicolor_to_S.italica_conserved.bed",
             "DGF_S.bicolor_to_B.distachyon_conserved.bed"))
             
#Files with crossmapped DGFs (occupied)
occupied<-data.frame("Occupied"=
          c("DGF_BS_S.bicolor_to_Z.mays_conserved_and_occupied.bed",
            "DGF_BS_S.bicolor_to_S.italica_conserved_and_occupied.bed",
            "DGF_BS_S.bicolor_to_B.distachyon_conserved_and_occupied.bed",
            "DGF_S.bicolor_to_Z.mays_conserved_and_occupied.bed",
            "DGF_S.bicolor_to_S.italica_conserved_and_occupied.bed",
            "DGF_S.bicolor_to_B.distachyon_conserved_and_occupied.bed"))

samples<-data.frame("Samples"=c("S. bicolor_Zm_BS",
                                "S. bicolor_Si_BS",
                                "S. bicolor_Bd_BS",
                                "S. bicolor_Zm_WL",
                                "S. bicolor_Si_WL",
                                "S. bicolor_Bd_WL"))

temp_df<-data.frame("Samples"=samples,
           "Observed_conserved"=apply(conserved,1,function(x) readData(x)),
           "Observed_conserved_and_occupied"=apply(occupied,1,function(x) readData(x)))

#Process the random data
tab_Sb_WL<-calcMeans(R_Sb_WL,"S. bicolor","WL")
tab_Sb_BS<-calcMeans(R_Sb_BS,"S. bicolor","BS")

#Combine data from the two cell types into a single df
tab_Sb<-rbind(tab_Sb_WL,tab_Sb_BS)

#Combined Observed and random data
tab_all<-merge(tab_Sb,temp_df,by="Samples")
tab_all<-merge(tab_all,Num_DGF,by="Tissue")
tab_all<-tab_all[,c("Species_1",
                    "Species_2",
                    "Tissue",
                    "Num_DGF",
                    "Comparison",
                    "Observed_conserved",
                    "Random_average_conserved",
                    "Observed_conserved_and_occupied",
                    "Random_average_conserved_and_occupied")]

#Calculate percentage conservation
tab_all$Observed_conserved_percent<-(tab_all$Observed_conserved/tab_all$Num_DGF)*100
tab_all$Random_conserved_percent<-(tab_all$Random_average_conserved/tab_all$Num_DGF)*100
tab_all$Observed_conserved_and_occupied_percent<-(tab_all$Observed_conserved_and_occupied/tab_all$Num_DGF)*100
tab_all$Random_conserved_and_occupied_percent<-(tab_all$Random_average_conserved_and_occupied/tab_all$Num_DGF)*100
#replace abbreviations for final table
tab_all<-data.frame(lapply(tab_all,function(x) {gsub("Bd", "B. distachyon", x)}))
tab_all<-data.frame(lapply(tab_all,function(x) {gsub("Si", "S. italica", x)}))
tab_all<-data.frame(lapply(tab_all,function(x) {gsub("Zm", "Z. mays", x)}))
#write results to file
write.table(tab_all,"Data_Figure_3C_table.csv",sep=",",row.names=F,quote=F)