#R code to make Figure S5C
library(gplots)

#Function to calculate enrichment
CalcEnrich<-function(DAP,ID){
  Comparison=DAP$Comparison
  five_prime_UTR=c(log(DAP$P.5.UTR/Distribution[c(ID),c("P.5UTR")],2))
  three_prime_UTR=c(log(DAP$P.3.UTR/Distribution[c(ID),c("P.3UTR")],2))
  CDS=c(log(DAP$P.CDS/Distribution[c(ID),c("P.CDS")],2))
  Promoter=c(log(DAP$P.Promoter/Distribution[c(ID),c("P.Promoter")],2))
  Intron=c(log(DAP$P.Intron/Distribution[c(ID),c("P.Intron")],2))
  Downstream=c(log(DAP$P.Downstream/Distribution[c(ID),c("P.Downstream")],2))
  Intergenic=c(log(DAP$P.Intergenic/Distribution[c(ID),c("P.Intergenic")],2))

  res<-data.frame("Comparison"=Comparison,
             "five_prime_UTR"=five_prime_UTR,
             "three_prime_UTR"=three_prime_UTR,
             "CDS"=CDS,
             "Promoter"=Promoter,
             "Intron"=Intron,
             "Downstream"=Downstream,
             "Intergenic"=Intergenic)
  return(res)
            
}

#Function to create heatmaps
createPlot<-function(x){
  #remove inifite values
  x<-x[!(x$three_prime_UTR=="-Inf"|x$five_prime_UTR=="-Inf" | x$CDS=="-Inf"| x$Intron=="-Inf"),]
  #create matrix
  x<-as.matrix(x[,c(2:8)])
  #create heatmap
  heatmap.2(x,dendrogram='none', Colv=F,trace='none')
}


#Read in data about distribution of DGF between genomic features
Species_dist<-read.table("Data_Figure_S5_genomic_feature_distribution.csv",head=T)
names<-Species_dist$Species
Distribution<-data.frame(t(Species_dist[,-c(1)]))
colnames(Distribution)<-names

#Read in data about distribution of motifs between genomic features
Sb_WL_DAP<-read.table("Data_Figure_S5_genomic_distribution_Sb_WL.csv",head=T)

#Calculate enrichment of motifs within features
E_Sb_WL<-CalcEnrich(Sb_WL_DAP,"Sb")

#Create heatmaps of enrichment
createPlot(E_Sb_WL)
