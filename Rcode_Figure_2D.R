#Code Required to make Figure 2C
require(ggplot2)

#Read in data
data<-read.table("~/R/DNAse_publication/Data_Figure_2C.csv",head=T,sep=",")
data$Type<-factor(data$Type, levels=c("Unique"))
data$Species<-factor(data$Species,levels=c("S. bicolor","Z. mays","S. italica"))

#Create Plot
pp<-ggplot(data, aes(Species, DGFs)) +   
  geom_bar(aes(fill = Type), position = "dodge", stat="identity") + 
  theme_minimal()

svg("Figure 2C.svg")
print(pp)
dev.off()