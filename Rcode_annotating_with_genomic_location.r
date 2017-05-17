# Script to annotate input files using ChIPpeakAnno
library('ChIPpeakAnno')
library("GenomicFeatures")
library(rtracklayer)
require(gdata)

#################################################
#Define Functions
#################################################
#Function to create a GRange object from a BED file.
bed2GRanges <- function(peaks) {
  #generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    ranges=IRanges(start=peaks$start, end=peaks$end, names=peaks$SiteName),
    strand="*"
  )
  return(myrange)
}

#Function to create GRange object from a GFF3 file.
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
# Program
#################################################
# Get options from commandline.
# ChrFile  = args[1]
# gff3File = args[2]
# BEDFILE = args[3]
# OUTFILE = args[4]

args<-commandArgs(trailingOnly = TRUE)

# Read in chromsome index information obtained using faSize (options --detailed)
print("Importing chromosome information.")
ChrFILE <- read.table(args[1], sep="\t", as.is=T, header=F)

# Read in GFF3 annotation file
# note input GFF3 file must include exon information, use gene.exon.gff3
print("Importing annotation GFF3 data.")
gff3<-read.table(args[2])

# Read in BED file for annotation.
print("Importing input BED data.")
bed <- read.table(args[3],sep="\t",  header=TRUE)

# Make a GRange object of annotations
# Create a named vector of chromosome sizes for GRange object creation.
print("Processing annotation data.")
ChrNames<-ChrFILE$V1
ChrSizes<-ChrFILE$V2
ChrSizesNamed<-ChrSizes
names(ChrSizesNamed)<-ChrNames

# Make a chromosome information file.
chrominfo <- data.frame(chrom = ChrNames,
                        length= ChrSizes,
                        is_circular= rep(FALSE, length(ChrNames)))

# Make GRanges object of GFF3.
GFF3_Granges<-GFF2GRanges(gff3)

# Make a GRanges object of input locations.
print("Processing input data.")

# Name input locations with row numbers for GRange object creation.
bed$V4<-rownames(bed)

# Name columns for GRange object creation.
col_names<-c("chr","start","end","SiteName")
names(bed)<-col_names

# Create GRanges object of input data.
bedOutput <- bed2GRanges(bed)


# Annotate data using GRanges object.
print("Annotating input data.")
# Settings
# "first" - only report first overlap
# overlapping - only annotate fragments that overlap
# maxgap - maximum distance between annotation (gene) and input.
x<-annotatePeakInBatch(bedOutput, AnnotationData=GFF3_Granges,select="first",output="overlapping",maxgap=1000L)                       

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
write.table(otable, args[4], row.names=FALSE,quote=F) 
