library('ChIPpeakAnno')
library("GenomicFeatures")
library(rtracklayer)
require(gdata)
##################################################
#BED annotation mapping
# A script to annotate genomic features represented
# in a bed format with their nearest gene, using a
# reference GFF3 file

##################################################
# Function to create a GRange object from a BED file
##################################################
bed2GRanges <- function(peaks) {
  #generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    ranges=IRanges(start=peaks$start, end=peaks$end, names=peaks$SiteName),
    strand="*"
  )
  return(myrange)
}

##################################################
# Function to create a GRange object from a GFF3 file
##################################################

GFF2GRanges <- function(peaks) {
  #generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$V1, #Names of all the features (Chr) peaks are mapped against
    ranges=IRanges(start=peaks$V4, end=peaks$V5, names=peaks$V9),
    strand="*",
    seqlengths=ChrSizesNamed #Lengths of all the features (Chr) peaks are mapped against
  )
  return(myrange)
}
#################################################
# Get Options
#note input GFF3 file must include exon information, use gene.exon.gff3

args<-commandArgs(trailingOnly = TRUE)
# ChrFile  = args[1]
# gff3File = args[2]
# BEDFILE = args[3]
# OUTFILE = args[4]
OUTFILE  <- args[4]

#Read in Chromsome information obtained using faSize --detailed
ChrFILE <- read.table(args[1], sep="\t", as.is=T, header=F)
ChrNames<-ChrFILE$V1
ChrSizes<-ChrFILE$V2
ChrSizesNamed<-ChrSizes
names(ChrSizesNamed)<-ChrNames

gff3<-read.table(args[2])
#Make GRanges object of GFF3
GFF3_Granges<-GFF2GRanges(gff3)

#Make a chromosome information file
chrominfo <- data.frame(chrom = ChrNames,
                        length= ChrSizes,
                        is_circular= rep(FALSE, length(ChrNames)))

## import the BED output
# Create a GRanges object
print("Importing BED data.")
bed <- read.table(args[3],sep="\t",  header=FALSE)
bed$V4<-rownames(bed)
col_names<-c("chr","start","end","SiteName")
names(bed)<-col_names
bedOutput <- bed2GRanges(bed)

#annotate peaks with GFF3 file
x<-annotatePeakInBatch(bedOutput, AnnotationData=GFF3_Granges,select="first",output="overlapping",maxgap=3000L)                       

#The program will not include peaks which are unannotated, to add them back in use the following code:
xx<-as.data.frame(unname(x))
xxx<-as.data.frame(bedOutput)

intragenic<-subset(xxx,!(row.names(xxx) %in% xx$peak))
intragenic$names<-row.names(intragenic)
tempDf<-data.frame(intragenic$seqnames,
                   intragenic$start,
                   intragenic$end,
                   intragenic$width,
                   "strand"=rep("*",length(intragenic$strand)),
                   intragenic$names,
                   "feature"=rep("intragenic",length(intragenic$names)),
                   "start_position"=rep("N/A",length(intragenic$names)),
                   "end_position"=rep("N/A",length(intragenic$names)),
                   "insideFeature"=rep("N/A",length(intragenic$names)),
                   "distancetoFeature"=rep("N/A",length(intragenic$names)),
                   "shortestDistance"=rep("N/A",length(intragenic$names)),
                   "fromOverlappingOrNearest"=rep("N/A",length(intragenic$names))
)
colnames(tempDf)<-colnames(xx)
otable<-rbind(xx,tempDf)
#write the results to file                   
write.table(otable, OUTFILE, row.names=FALSE,quote=F) 
