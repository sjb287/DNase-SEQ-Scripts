#Rcode to make Figure 4C plot
require(ggplot2)
require(reshape2)
#load in data
sample_data<-read.table("~/R/DNAse_publication/Data_Figure_4C.txt",head=T,sep=" ")
sample_data<-melt(sample_data)
colnames(sample_data)<-c("Sample","Type","Percent")

sample_data$Type<-factor(sample_data$Type,levels=c("unannotated","de_novo_only","known"))

#plot data
ggplot(sample_data,aes(x = Sample,y=Percent, fill=Type)) + geom_bar(stat="identity",colour="black") + scale_fill_manual(values=c("bisque2","firebrick2","steelblue4")) +theme_minimal()

#500x350