library(RColorBrewer)
require(plyr)
#################################################
createPieChart<-function(df,title){
  df$Feature<-factor(df$Feature, levels=c("Intergenic.Region","Promoters","fiveUTRs","Exons","Introns","threeUTRs","immediateDownstream"))
  df<-arrange(df, Feature)
  x<-round(df$Percentage,1)
  
  # Give the chart file a name.
  svg(file = title,width=4,height=4)
  # Plot the chart.
  pie(x, x, col = colPalette, clockwise=TRUE, init.angle=90)
  
  # Save the file.
  dev.off()
}
#################################################
#Set colour palette
colPalette <- brewer.pal(7, "Set1")

#Set working directory
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_distribution/Files_distribution_DHS/")

#Define Input Files
Sb_WL<-read.table("Sb_WL_DHS_localisation.sum",head=T)
Sb_BS<-read.table("Sb_BS_DHS_localisation.sum",head=T)
Si_WL<-read.table("Si_WL_DHS_localisation.sum",head=T)
Si_BS<-read.table("Si_BS_DHS_localisation.sum",head=T)
Bd_WL<-read.table("Bd_WL_DHS_localisation.sum",head=T)
Zm_BS<-read.table("Zm_BS_DHS_localisation.sum",head=T)
Zm_WL<-read.table("Zm_WL_DHS_localisation.sum",head=T)
BS_4_way<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/BS_4_WAY_DGF_localisation.sum", head=T)
WL_4_way<-read.table("/home/sjb287/documents/DNase/Files_Publication/Files_submission/WL_4_WAY_DGF_localisation.sum", head=T)

createPieChart(Sb_WL,"Sb_WL_DHS_localisation.svg")
createPieChart(Sb_BS,"Sb_BS_DHS_localisation.svg")
createPieChart(Bd_WL,"Bd_WL_DHS_localisation.svg")
createPieChart(Zm_WL,"Zm_WL_DHS_localisation.svg")
createPieChart(Zm_BS,"Zm_BS_DHS_localisation.svg")
createPieChart(Si_WL,"Si_WL_DHS_localisation.svg")
createPieChart(Si_BS,"Si_BS_DHS_localisation.svg")

createPieChart(BS_4_way,"BS_4_way_DGF_localisation.svg")
createPieChart(WL_4_way,"WL_4_way_DGF_localisation.svg")

setwd("/home/sjb287/documents/DNase/Files_Publication/Files_distribution/Files_distribution_DGF/")
Sb_WL<-read.table("Sb_WL_DGF_localisation.sum",head=T)
Sb_BS<-read.table("Sb_BS_DGF_localisation.sum",head=T)
Si_WL<-read.table("Si_WL_DGF_localisation.sum",head=T)
Si_BS<-read.table("Si_BS_DGF_localisation.sum",head=T)
Bd_WL<-read.table("Bd_WL_DGF_localisation.sum",head=T)
Zm_WL<-read.table("Zm_WL_DGF_localisation.sum",head=T)
Zm_BS<-read.table("Zm_BS_DGF_localisation.sum",head=T)
createPieChart(Sb_WL,"Sb_WL_DGF_localisation.svg")
createPieChart(Sb_BS,"Sb_BS_DGF_localisation.svg")
createPieChart(Bd_WL,"Bd_WL_DGF_localisation.svg")
createPieChart(Zm_WL,"Zm_WL_DGF_localisation.svg")
createPieChart(Zm_BS,"Zm_BS_DGF_localisation.svg")
createPieChart(Si_WL,"Si_WL_DGF_localisation.svg")
createPieChart(Si_BS,"Si_BS_DGF_localisation.svg")
