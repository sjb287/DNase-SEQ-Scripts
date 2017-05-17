#Code to generate data for Figure 5
#################################################
#Define functions
#################################################
subSetGenes<-function(gene,x,n){
  
  #subset all interactions corresponding to gene
  interactions<-subset(x,x[,4]==gene)  
  #create outfile name
  ofile<-paste(n,gene,"bed",sep=".")
  #write list of interactions file
  write.table(interactions,ofile,row.names=F,quote=F,sep="\t")
  
}

extractGeneList<-function(GOI,prefix){
  #create directory to store results
  dir.create(prefix)
  #move into directory
  setwd(prefix)
  #extract list of unique TFs
  TFs<-data.frame(unique(GOI[,4]))
  #write outfile to each
  apply(TFs,1,function(x) subSetGenes(x,GOI,prefix))
  
}
#################################################
#read in files
Zm_WL<-read.table("Data_Figure_5_motifs_Zm_WL.bed")
Zm_BS<-read.table("Data_Figure_5_motifs_Zm_BS.bed")
Sb_WL<-read.table("Data_Figure_5_motifs_Sb_WL.bed")
Sb_BS<-read.table("Data_Figure_5_motifs_Sb_BS.bed")
Si_WL<-read.table("Data_Figure_5_motifs_Si_WL.bed")
Si_BS<-read.table("Data_Figure_5_motifs_Si_BS.bed")
Bd_WL<-read.table("Data_Figure_5_motifs_Bd_WL.bed")


#Run Script
extractGeneList(Zm_WL,"Zm_WL")
extractGeneList(Zm_BS,"Zm_BS")
extractGeneList(Sb_WL,"Sb_WL")
extractGeneList(Sb_BS,"Sb_BS")
extractGeneList(Bd_WL,"Bd_WL")
extractGeneList(Si_BS,"Si_BS")
extractGeneList(Si_WL,"Si_WL")