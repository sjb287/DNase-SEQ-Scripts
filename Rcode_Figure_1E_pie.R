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

#Set file destination
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/")

#read in files
Sb_WL<-read.table("Figure_1E_Sb_WL_DGF_localisation.sum",head=T)
Si_WL<-read.table("Figure_1E_Si_WL_DGF_localisation.sum",head=T)
Bd_WL<-read.table("Figure_1E_Bd_WL_DGF_localisation.sum",head=T)
Zm_WL<-read.table("Figure_1E_Zm_WL_DGF_localisation.sum",head=T)

#create pie charts
createPieChart(Sb_WL,"Sb_WL_DGF_localisation.svg")
createPieChart(Bd_WL,"Bd_WL_DGF_localisation.svg")
createPieChart(Zm_WL,"Zm_WL_DGF_localisation.svg")
createPieChart(Si_WL,"Si_WL_DGF_localisation.svg")
