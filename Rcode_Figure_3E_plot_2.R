#Code required to make images used in Figure 5
library(genoPlotR)
require(reshape2)
create_dna_segs_gene<-function(input){
  
  dna_segs_gene<-data.frame("name"=input[,1],
                            "start"=as.numeric(input[,2]),
                            "end"=as.numeric(input[,3]),
                            "strand"=input[,4],
                            "col"=rep("grey",length(input[,1])),
                            "lty"=rep(1,length(input[,1])),
                            "pch"=rep(1,length(input[,1])),                      
                            "cex"=rep(8,length(input[,1])),
                            "gene_type"=rep("arrows",length(input[,1])))
  return(dna_segs_gene)
}


create_dna_segs_DGFs<-function(input){
  
dna_segs_DGFs<-data.frame("name"=input$V4,
                          "start"=as.numeric(input$V2),
                          "end"=as.numeric(input$V3),
                          "strand"="+",
                          "col"=rep("grey",length(input[,1])),
                          "lty"=rep(1,length(input[,1])),
                          "pch"=rep(1,length(input[,1])),                      
                          "cex"=rep(8,length(input[,1])),
                          "gene_type"=rep("blocks",length(input[,1])))
  return(dna_segs_DGFs)

}

create_dna_segs_hits<-function(input){
  
  dna_segs_hits<-data.frame("name"=input[,3],
                            "start"=as.numeric(as.character(input[,1])),
                            "end"=as.numeric(as.character(input[,2])),
                            "strand"="+",
                            "col"=rep("grey",length(input[,1])),
                            "lty"=rep(1,length(input[,1])),
                            "pch"=rep(1,length(input[,1])),                      
                            "cex"=rep(8,length(input[,1])),
                            "gene_type"=rep("bars",length(input[,1])))
  return(dna_segs_hits)
}

#select columns of interest and get DGFs which contain functional equivalents
#this is also used later to create the comparisons vector
createComparisons<-function(y){
  
  y<-unique(y[,c(2,3,4,6,7,8)])
  temp_comparisons<-y[,c(2,3,5,6)]
  colnames(temp_comparisons)<-c("start1","end1","start2","end2")
  temp_comparisons$col<-rep("red",length(temp_comparisons[,1]))
  temp_comparisons[,c(1)]<-as.numeric(as.character(temp_comparisons[,c(1)]))
  temp_comparisons[,c(2)]<-as.numeric(as.character(temp_comparisons[,c(2)]))
  temp_comparisons[,c(3)]<-as.numeric(as.character(temp_comparisons[,c(3)]))
  temp_comparisons[,c(4)]<-as.numeric(as.character(temp_comparisons[,c(4)]))
  return(temp_comparisons)
  
}

createComparisonsDGFs<-function(temp_comparisons,colour){
  
  colnames(temp_comparisons)<-c("start1","end1","start2","end2")
  temp_comparisons$col<-rep(colour,length(temp_comparisons[,1]))
  temp_comparisons[,c(1)]<-as.numeric(as.character(temp_comparisons[,c(1)]))
  temp_comparisons[,c(2)]<-as.numeric(as.character(temp_comparisons[,c(2)]))
  temp_comparisons[,c(3)]<-as.numeric(as.character(temp_comparisons[,c(3)]))
  temp_comparisons[,c(4)]<-as.numeric(as.character(temp_comparisons[,c(4)]))
  return(temp_comparisons)  
}

createDNAsegs<-function(temp,temp_DGFs,hits){
  #takes 3 inputs a list of gene identifers for
  #select only relavant columns from the hits 
  hits<-unique(hits[,c(2,3,5)])
  
  #Create segments for gene (a), DGFs (b), functional hits (c)
  if(nrow(temp_DGFs) == 0 & nrow(hits) == 0){
    
    temp_segs_1a<-create_dna_segs_gene(temp)
    #Combine segments and make a dna_segs object
    temp_segs_1<-dna_seg(unique(temp_segs_1a))
    
  }else if(nrow(hits)==0){
    
    temp_segs_1a<-create_dna_segs_gene(temp)
    temp_segs_1b<-create_dna_segs_DGFs(temp_DGFs)
    #Combine segments and make a dna_segs object
    temp_segs_1<-dna_seg(unique(rbind(temp_segs_1a,temp_segs_1b)))
    
  }else{

    temp_segs_1a<-create_dna_segs_gene(temp)
    temp_segs_1b<-create_dna_segs_DGFs(temp_DGFs)
    temp_segs_1c<-create_dna_segs_hits(hits)
    #Combine segments and make a dna_segs object
    temp_segs_1<-dna_seg(unique(rbind(temp_segs_1a,temp_segs_1b,temp_segs_1c)))
  }
  
  
  #return the dna segment objects
  return(temp_segs_1)
}

#function to identify footprints which are conserved and occupied
getConserved<-function(gene,dgf,mapped){
  
  #gene = information on gene mapped to [ID] [start] [stop] [strand] [feature] [Chromosome].
  #dgf = all the footprints from species 2 linked to the gene of interest.
  #mapped = all the locations of footprints from species 1 mapped to species 2.
  
  #get the chromosome on which the gene is located
  #this is necessary to limit search based on stop/start position
  chr<-gene$V6
  
  #isolate only those mapped footprints on the same chromosome
  #necessary to make sure there is no ambiguity when searching on end bp.
  mapped<-subset(mapped,mapped[,1]==chr)
  
  #keep only the information about chr,start,stop,name
  dgf<-dgf[,c(1:4)]
  
  #extract all the rows which contain overlapping DGFs 
  vals<-apply(dgf,1,function(x) mapped[(x[2]<mapped[,3] & x[3]>mapped[,3] | x[2]<mapped[,2] & x[3]>mapped[,2]),])
  
  #these are retuned as a list, to convert to a dataframe it must go through a matrix,
  #the matrix is transposed then converted to a dataframe
  results<-as.data.frame(melt(vals,id=c("V1","V2","V3","V4")))
  results<-results[,c(1:4)]
  return(results)
}
############################################################################################################
#Set wotrking ditectory
setwd("/home/sjb287/documents/DNase/Files_Publication/Files_submission/")

#Gene annotation files - custom altered GFF3 format
annotation_1="Data_Figure_5_genoplot_Sb_annotations.txt"
annotation_2="Data_Figure_5_genoplot_Zm_annotations.txt"
annotation_3="Data_Figure_5_genoplot_Si_annotations.txt"
annotation_4="Data_Figure_5_genoplot_Bd_annotations.txt"

args<-commandArgs(trailingOnly = TRUE)
flag<-args[1] #variable to determine whether looking at WL or BS
gene_1<-args[2] #Sb gene ID
gene_2<-args[3] #Zm gene ID
gene_3<-args[4] #Si gene ID
gene_4<-args[5] #Bd gene ID
outfile<-args[6] #outfile

if(flag=="WL"){
  #DGF files
  DGFs_1="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Sb_WL.bed"
  DGFs_2="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Zm_WL.bed"
  DGFs_3="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Si_WL.bed"
  DGFs_4="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Bd_WL.bed"
  #Hit annotation files
  hits_1<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Sb_WL_DGF_annotated_known_and_de_novo.bed"
  hits_2<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Zm_WL_DGF_annotated_known_and_de_novo.bed"
  hits_3<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Si_WL_DGF_annotated_known_and_de_novo.bed"
  hits_4<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Bd_WL_DGF_annotated_known_and_de_novo.bed"
  #Sb to Zm
  cross_mapped_1_all_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_S.bicolor_to_Z.mays_mapped.bed"
  cross_mapped_1_all_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_Z.mays_to_S.bicolor_mapped.bed"
  cross_mapped_1_conserved_pairwise_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_S.bicolor_to_Z.mays_conserved_and_occupied_.bed"
  cross_mapped_1_conserved_pairwise_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_Z.mays_to_S.bicolor_conserved_and_occupied_.bed"
  #Zm to Si
  cross_mapped_2_all_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_Z.mays_to_S.italica_mapped.bed"
  cross_mapped_2_all_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_S.italica_to_Z.mays_mapped.bed"
  cross_mapped_2_conserved_pairwise_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_Z.mays_to_S.italica_conserved_and_occupied_.bed"
  cross_mapped_2_conserved_pairwise_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_S.italica_to_Z.mays_conserved_and_occupied_.bed"
  #Si to Bd
  cross_mapped_3_all_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_S.italica_to_B.distachyon_mapped.bed"
  cross_mapped_3_all_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_B.distachyon_to_S.italica_mapped.bed"
  cross_mapped_3_conserved_pairwise_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_S.italica_to_B.distachyon_conserved_and_occupied_.bed"
  cross_mapped_3_conserved_pairwise_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_WL_B.distachyon_to_S.italica_conserved_and_occupied_.bed"
}else{
  #DGF files
  DGFs_1="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Sb_BS.bed"
  DGFs_2="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Zm_BS.bed"
  DGFs_3="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Si_BS.bed"
  DGFs_4="/data/reads/DNase_data/2015_DNase_Hibberd_Monocot/Files_Publication/sjb287/DGF_Bd_WL.bed"
  #Hit annotation files
  hits_1<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Sb_BS_DGF_annotated_known_and_de_novo.bed"
  hits_2<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Zm_BS_DGF_annotated_known_and_de_novo.bed"
  hits_3<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Si_BS_DGF_annotated_known_and_de_novo.bed"
  hits_4<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DAP/Bd_WL_DGF_annotated_known_and_de_novo.bed"
  #Sb to Zm
  cross_mapped_1_all_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_S.bicolor_to_Z.mays_mapped.bed"
  cross_mapped_1_all_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_Z.mays_to_S.bicolor_mapped.bed"
  cross_mapped_1_conserved_pairwise_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_S.bicolor_to_Z.mays_conserved_and_occupied_.bed"
  cross_mapped_1_conserved_pairwise_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_Z.mays_to_S.bicolor_conserved_and_occupied_.bed"
  #Zm to Si
  cross_mapped_2_all_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_Z.mays_to_S.italica_mapped.bed"
  cross_mapped_2_all_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_S.italica_to_Z.mays_mapped.bed"
  cross_mapped_2_conserved_pairwise_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_Z.mays_to_S.italica_conserved_and_occupied_.bed"
  cross_mapped_2_conserved_pairwise_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_S.italica_to_Z.mays_conserved_and_occupied_.bed"
  #Si to Bd
  cross_mapped_3_all_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_S.italica_to_B.distachyon_mapped.bed"
  cross_mapped_3_all_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_B.distachyon_to_S.italica_mapped.bed"
  cross_mapped_3_conserved_pairwise_A_to_B<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_S.italica_to_B.distachyon_conserved_and_occupied_.bed"
  cross_mapped_3_conserved_pairwise_B_to_A<-"/home/sjb287/documents/DNase/Files_Publication/Files_submission/DGF_BS_B.distachyon_to_S.italica_conserved_and_occupied_.bed"
  
}
#################################################
#Producing dna_segs vector
#################################################
#Part I: Process Positional Conservation
#################################################
#read in annotation files - all gene locations
#columns: V1(GeneID), V2(start bp), V3 (end bp), V4(strand), V5(feature),V6(chromosome)
annot_1<-read.table(annotation_1,stringsAsFactors=FALSE)
annot_2<-read.table(annotation_2,stringsAsFactors=FALSE)
annot_3<-read.table(annotation_3,stringsAsFactors=FALSE)
annot_4<-read.table(annotation_4,stringsAsFactors=FALSE)
#update Si annotations
Si_Syn<-read.table("/data/genomes/setaria/Sitalica/Sitalica_312_v2.2.synonym_loci.txt")
colnames(Si_Syn)<-c("ID","V1")
Si_Syn$V1<-paste(Si_Syn$V1,".g",sep="")
annot_3<-merge(annot_3,Si_Syn,by="V1")
annot_3<-annot_3[,c("ID","V2","V3","V4","V5","V6")]
colnames(annot_3)[1]<-"V1"
#produce temporary vectors containing
#subset annotation file to contain gene information for the input ID
#columns: V1(GeneID), V2(start bp), V3 (end bp), V4(strand), V5(feature),V6(chromosome)
temp_1<-subset(annot_1,annot_1$V1==gene_1)
temp_2<-subset(annot_2,annot_2$V1==gene_2)
temp_3<-subset(annot_3,annot_3$V1==gene_3)
temp_4<-subset(annot_4,annot_4$V1==gene_4)
#processDGFs
#Read in DGF information - all DGF locations
#Columns: (chromosome), (DGF start), (DGF end), (Name), (DGF score), (strand)
annot_DGF_1<-read.table(DGFs_1,stringsAsFactors=FALSE)
annot_DGF_2<-read.table(DGFs_2,stringsAsFactors=FALSE)
annot_DGF_3<-read.table(DGFs_3,stringsAsFactors=FALSE)
annot_DGF_4<-read.table(DGFs_4,stringsAsFactors=FALSE)
#First only keep DGF sequences on the same chromosome as the gene of interest
annot_DGF_1<-subset(annot_DGF_1,annot_DGF_1$V1 == temp_1$V6)
annot_DGF_2<-subset(annot_DGF_2,annot_DGF_2$V1 == temp_2$V6)
annot_DGF_3<-subset(annot_DGF_3,annot_DGF_3$V1 == temp_3$V6)
annot_DGF_4<-subset(annot_DGF_4,annot_DGF_4$V1 == temp_4$V6)
#extract all DGFs within 3kb upstream and 1kb downstream of gene
temp_DGFs_1<-subset(annot_DGF_1,annot_DGF_1$V2 > temp_1$V2-3000 & annot_DGF_1$V2 < temp_1$V3+3000)
temp_DGFs_2<-subset(annot_DGF_2,annot_DGF_2$V2 > temp_2$V2-3000 & annot_DGF_2$V2 < temp_2$V3+3000)
temp_DGFs_3<-subset(annot_DGF_3,annot_DGF_3$V2 > temp_3$V2-3000 & annot_DGF_3$V2 < temp_3$V3+3000) 
temp_DGFs_4<-subset(annot_DGF_4,annot_DGF_4$V2 > temp_4$V2-3000 & annot_DGF_4$V2 < temp_4$V3+3000)
#processCrossMapped comparison (1) Sb/Zm
#read in the cross mapping files location of all DGFs that crossmapped
annot_cross_mapped_1a<-read.table(cross_mapped_1_all_A_to_B,stringsAsFactors=FALSE)
annot_cross_mapped_1b<-read.table(cross_mapped_1_all_B_to_A,stringsAsFactors=FALSE)
annot_cross_mapped_1c<-read.table(cross_mapped_1_conserved_pairwise_B_to_A,stringsAsFactors=FALSE)
annot_cross_mapped_1d<-read.table(cross_mapped_1_conserved_pairwise_A_to_B,stringsAsFactors=FALSE)
#processCrossMapped Zm/Si
annot_cross_mapped_2a<-read.table(cross_mapped_2_all_A_to_B,stringsAsFactors=FALSE)
annot_cross_mapped_2b<-read.table(cross_mapped_2_all_B_to_A,stringsAsFactors=FALSE)
annot_cross_mapped_2c<-read.table(cross_mapped_2_conserved_pairwise_B_to_A,stringsAsFactors=FALSE)
annot_cross_mapped_2d<-read.table(cross_mapped_2_conserved_pairwise_A_to_B,stringsAsFactors=FALSE)
#processCrossMapped Si/Bd
annot_cross_mapped_3a<-read.table(cross_mapped_3_all_A_to_B,stringsAsFactors=FALSE)
annot_cross_mapped_3b<-read.table(cross_mapped_3_all_B_to_A,stringsAsFactors=FALSE)
annot_cross_mapped_3c<-read.table(cross_mapped_3_conserved_pairwise_B_to_A,stringsAsFactors=FALSE)
annot_cross_mapped_3d<-read.table(cross_mapped_3_conserved_pairwise_A_to_B,stringsAsFactors=FALSE)
#extract all DGFs within 3kb upstream and 1kb downstream of gene, remember the gene mapping to is used as a comparison
temp_annot_cross_mapped_1a<-subset(annot_cross_mapped_1a,annot_cross_mapped_1a$V2 > temp_2$V2-3000 & annot_cross_mapped_1a$V2 < temp_2$V3+1000)
temp_annot_cross_mapped_1b<-subset(annot_cross_mapped_1b,annot_cross_mapped_1b$V2 > temp_1$V2-3000 & annot_cross_mapped_1b$V2 < temp_1$V3+1000)
temp_annot_cross_mapped_1c<-subset(annot_cross_mapped_1c,annot_cross_mapped_1c$V2 > temp_2$V2-3000 & annot_cross_mapped_1c$V2 < temp_2$V3+1000)
temp_annot_cross_mapped_1d<-subset(annot_cross_mapped_1d,annot_cross_mapped_1d$V2 > temp_2$V2-3000 & annot_cross_mapped_1d$V2 < temp_2$V3+1000)
#extract all DGFs within 3kb upstream and 1kn downstream of gene, remember the gene mapping to is used as a comparison
temp_annot_cross_mapped_2a<-subset(annot_cross_mapped_2a,annot_cross_mapped_2a$V2 > temp_3$V2-1000 & annot_cross_mapped_2a$V2 < temp_3$V3+3000)
temp_annot_cross_mapped_2b<-subset(annot_cross_mapped_2b,annot_cross_mapped_2b$V2 > temp_2$V2-3000 & annot_cross_mapped_2b$V2 < temp_2$V3+1000)
temp_annot_cross_mapped_2c<-subset(annot_cross_mapped_2c,annot_cross_mapped_2c$V2 > temp_3$V2-1000 & annot_cross_mapped_2c$V2 < temp_3$V3+3000)
temp_annot_cross_mapped_2d<-subset(annot_cross_mapped_2d,annot_cross_mapped_2d$V2 > temp_3$V2-1000 & annot_cross_mapped_2d$V2 < temp_3$V3+3000)
#extract all DGFs within 3kb upstream and 1kn downstream of gene, remember the gene mapping to is used as a comparison
temp_annot_cross_mapped_3a<-subset(annot_cross_mapped_3a,annot_cross_mapped_3a$V2 > temp_4$V2-3000 & annot_cross_mapped_3a$V2 < temp_4$V3+1000)
temp_annot_cross_mapped_3b<-subset(annot_cross_mapped_3b,annot_cross_mapped_3b$V2 > temp_3$V2-1000 & annot_cross_mapped_3b$V2 < temp_3$V3+3000)
temp_annot_cross_mapped_3c<-subset(annot_cross_mapped_3c,annot_cross_mapped_3c$V2 > temp_4$V2-3000 & annot_cross_mapped_3c$V2 < temp_4$V3+1000)
temp_annot_cross_mapped_3d<-subset(annot_cross_mapped_3d,annot_cross_mapped_3d$V2 > temp_4$V2-3000 & annot_cross_mapped_3d$V2 < temp_4$V3+1000)

#################################################
#Part II: Process Functional Annotations
#################################################
#read in all hits keeping only those annotated to the gene of interest
annot_hits_1<-read.table(hits_1,stringsAsFactors=FALSE)
annot_hits_2<-read.table(hits_2,stringsAsFactors=FALSE)
annot_hits_3<-read.table(hits_3,stringsAsFactors=FALSE)
annot_hits_4<-read.table(hits_4,stringsAsFactors=FALSE)
temp_hits_1<-subset(annot_hits_1,annot_hits_1$V5==gene_1)
temp_hits_2<-subset(annot_hits_2,annot_hits_2$V5==gene_2)
temp_hits_3<-subset(annot_hits_3,annot_hits_3$V5==gene_3)
temp_hits_4<-subset(annot_hits_4,annot_hits_4$V5==gene_4)
#find which motif hits overlap based on TF for comparison 1
x1<-unique(subset(temp_hits_1,temp_hits_1$V4 %in% temp_hits_2$V4))
x2<-unique(subset(temp_hits_2,temp_hits_2$V4 %in% temp_hits_1$V4))
#find which motif hits overlap based on TF for comparison 2
x3<-unique(subset(temp_hits_2,temp_hits_2$V4 %in% temp_hits_3$V4))
x4<-unique(subset(temp_hits_3,temp_hits_3$V4 %in% temp_hits_2$V4))
#find which motif hits overlap based on TF for comparison 3
x5<-unique(subset(temp_hits_3,temp_hits_3$V4 %in% temp_hits_4$V4))
x6<-unique(subset(temp_hits_4,temp_hits_4$V4 %in% temp_hits_3$V4))
#merge datasets to create comparison vectors
y1<-merge(x1,x2,by="V4")
y2<-merge(x3,x4,by="V4")
y3<-merge(x5,x6,by="V4")
q<-merge(y1,y2,by="V4")
temp_comp_1<-createComparisons(y1)
temp_comp_2<-createComparisons(y2)
temp_comp_3<-createComparisons(y3)

#################################################
#Identify crossmapped DGFs using names
#Comparison 1
DGF_A_to_B_1<-merge(annot_cross_mapped_1a,temp_DGFs_1, by="V4")
DGF_B_to_A_1<-merge(annot_cross_mapped_1b,temp_DGFs_2, by="V4")
#get those footprints which are conserved between species 1 and 2
temp_crossmapped_1<-getConserved(temp_2,temp_DGFs_2,annot_cross_mapped_1a)
#with the conserved DGFs we only need the ones originating from the A partner
DGF_A_to_B_1_conserved<-merge(temp_crossmapped_1,temp_DGFs_1,by="V4")
#remove any duplicate entries
DGF_A_to_B_1<-subset(DGF_A_to_B_1,!(DGF_A_to_B_1[,6] %in% DGF_A_to_B_1_conserved[,6]))

#Comparison 2
DGF_A_to_B_2<-merge(annot_cross_mapped_2a,temp_DGFs_2, by="V4")
DGF_B_to_A_2<-merge(annot_cross_mapped_2b,temp_DGFs_3, by="V4")

temp_crossmapped_2<-getConserved(temp_3,temp_DGFs_3,annot_cross_mapped_2a)
DGF_A_to_B_2_conserved<-merge(temp_crossmapped_2,temp_DGFs_2,by="V4")
#with the conserved DGFs we only need the ones originating from the A partner
#remove any duplicate entries
DGF_A_to_B_2<-subset(DGF_A_to_B_2,!(DGF_A_to_B_2[,6] %in% DGF_A_to_B_2_conserved[,6]))

#Comparison 3
DGF_A_to_B_3<-merge(annot_cross_mapped_3a,temp_DGFs_3, by="V4")
DGF_B_to_A_3<-merge(annot_cross_mapped_3b,temp_DGFs_4, by="V4")
temp_crossmapped_3<-getConserved(temp_4,temp_DGFs_4,annot_cross_mapped_3a)
DGF_A_to_B_3_conserved<-merge(temp_crossmapped_3,temp_DGFs_3,by="V4")
#with the conserved DGFs we only need the ones originating from the A partner
DGF_A_to_B_3<-subset(DGF_A_to_B_3,!(DGF_A_to_B_3[,6] %in% DGF_A_to_B_3_conserved[,6]))

#keep only columns of interest
DGF_A_to_B_1<-unique(DGF_A_to_B_1[,c(6,7,3,4)])
DGF_A_to_B_2<-unique(DGF_A_to_B_2[,c(6,7,3,4)])
DGF_A_to_B_3<-unique(DGF_A_to_B_3[,c(6,7,3,4)])
DGF_A_to_B_1_conserved<-DGF_A_to_B_1_conserved[,c(6,7,3,4)]
DGF_A_to_B_2_conserved<-DGF_A_to_B_2_conserved[,c(6,7,3,4)]
DGF_A_to_B_3_conserved<-DGF_A_to_B_3_conserved[,c(6,7,3,4)]
DGF_B_to_A_1<-unique(DGF_B_to_A_1[,c(3,4,6,7)])
DGF_B_to_A_2<-unique(DGF_B_to_A_2[,c(3,4,6,7)])
DGF_B_to_A_3<-unique(DGF_B_to_A_3[,c(3,4,6,7)])
#Create df to make comparisons object
#only functional annotations are marked in red
df_DGF_A_to_B_1<-createComparisonsDGFs(DGF_A_to_B_1,"grey")
df_DGF_A_to_B_2<-createComparisonsDGFs(DGF_A_to_B_2,"grey")
df_DGF_A_to_B_3<-createComparisonsDGFs(DGF_A_to_B_3,"grey")
df_DGF_A_to_B_1c<-createComparisonsDGFs(DGF_A_to_B_1_conserved,"red")
df_DGF_A_to_B_2c<-createComparisonsDGFs(DGF_A_to_B_2_conserved,"red")
df_DGF_A_to_B_3c<-createComparisonsDGFs(DGF_A_to_B_3_conserved,"red")
df_DGF_B_to_A_1<-createComparisonsDGFs(DGF_B_to_A_1,"grey")
df_DGF_B_to_A_2<-createComparisonsDGFs(DGF_B_to_A_2,"grey")
df_DGF_B_to_A_3<-createComparisonsDGFs(DGF_B_to_A_3,"grey")
#combine different cross mapping features into a single dataframe for each comparison
comp_1<-rbind(temp_comp_1,df_DGF_A_to_B_1,df_DGF_B_to_A_1,df_DGF_A_to_B_1c)
comp_2<-rbind(temp_comp_2,df_DGF_A_to_B_2,df_DGF_B_to_A_2,df_DGF_A_to_B_2c)
comp_3<-rbind(temp_comp_3,df_DGF_A_to_B_3,df_DGF_B_to_A_3,df_DGF_A_to_B_3c)
#convert to comparison objects
comp_1<-comparison(comp_1)
comp_2<-comparison(comp_2)
comp_3<-comparison(comp_3)
#finally create a list of comparison objects for plotting
comps<-list(comp_1,comp_2,comp_3)
#create DNA seg objects
segs_1<-createDNAsegs(temp_1,temp_DGFs_1,temp_hits_1)
segs_2<-createDNAsegs(temp_2,temp_DGFs_2,temp_hits_2)
segs_3<-createDNAsegs(temp_3,temp_DGFs_3,temp_hits_3)
segs_4<-createDNAsegs(temp_4,temp_DGFs_4,temp_hits_4)
#finally create a list of dna_seg objects for plotting
segs<-list(segs_1,segs_2,segs_3,segs_4)

#################################################

#genoplot takes two input vectors
#DNA segs - genomic locations
#comparisons - mapping locations

gene_1<-subset(segs[[1]],segs[[1]]$gene_type=="arrows")
gene_2<-subset(segs[[2]],segs[[2]]$gene_type=="arrows")
gene_3<-subset(segs[[3]],segs[[3]]$gene_type=="arrows")
gene_4<-subset(segs[[4]],segs[[4]]$gene_type=="arrows")
max_1<-gene_1$end-gene_1$start
max_2<-gene_2$end-gene_2$start
max_3<-gene_3$end-gene_3$start
max_4<-gene_4$end-gene_4$start

checkStrand<-function(gene){
  if(gene$strand=="1"){
    gene_start<-gene$start
    gene_end<-gene$end
  }else{
    gene_start<-gene$end
    gene_end<-gene$start
  }
 results<-list(gene_start,gene_end)
 return(results)
}

gene_1_se<-checkStrand(gene_1)
gene_2_se<-checkStrand(gene_2)
gene_3_se<-checkStrand(gene_3)
gene_4_se<-checkStrand(gene_4)

#calculate the maximum fragment length
max_length<-max(c(max_1,max_2,max_3,max_4))
max_length<-max_length+3000
lims <- list(c(gene_1_se[[1]]-as.numeric(gene_1$strand)*(max_length-max_1)/2, gene_1_se[[2]]+as.numeric(gene_1$strand)*(max_length-max_1)/2),
             c(gene_2_se[[1]]-as.numeric(gene_2$strand)*(max_length-max_2)/2, gene_2_se[[2]]+as.numeric(gene_2$strand)*(max_length-max_2)/2),
             c(gene_3_se[[1]]-as.numeric(gene_3$strand)*(max_length-max_3)/2, gene_3_se[[2]]+as.numeric(gene_3$strand)*(max_length-max_3)/2),
             c(gene_4_se[[1]]-as.numeric(gene_4$strand)*(max_length-max_4)/2, gene_4_se[[2]]+as.numeric(gene_4$strand)*(max_length-max_4)/2)) 

#Plot results
out_table<-paste(outfile,"csv",sep=".")
outfile<-paste(outfile,"svg",sep=".")
o_t<-rbind(comps[[1]],comps[[2]],comps[[3]])
write.table(o_t,out_table,row.names=F,quote=F)
svg(outfile,width=5,height=2)
names(segs)<-c("S. bicolor","Z. mays","S. italica", "B. distachyon")
plot_gene_map(dna_segs = segs, comparisons = comps, xlims = lims)
dev.off()