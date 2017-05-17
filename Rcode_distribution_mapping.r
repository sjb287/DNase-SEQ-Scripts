library('ChIPpeakAnno')
library("GenomicFeatures")
library(rtracklayer)
require(gdata)

#################################################
#A Function to create a GRange object from the MACS output file
bed2GRanges <- function(peaks) {
  #generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr, peaks$start,sep=":")),
    strand="*",
    seqlengths=ChrSizesNamed
  )
  return(myrange)
}
#################################################
# Get Options
#note input GFF3 file must include exon information, use gene.exon.gff3

args<-commandArgs(trailingOnly = TRUE)
# ChrFile  = args[1]
# gff3File = args[2]
# Species  = args[3]
# BEDFILE  = args[4]

OUTFILE<-paste(args[3],"DGF_localisation.sum", sep = "_")

##################################################
#Read in Chromsome information obtained using faSize --detailed
##################################################
ChrFILE <- read.table(args[1], sep="\t", as.is=T, header=F)
ChrNames<-ChrFILE$V1
ChrSizes<-ChrFILE$V2
ChrSizesNamed<-ChrSizes
names(ChrSizesNamed)<-ChrNames

##################################################
#Make a chromosome information file
chrominfo <- data.frame(chrom = ChrNames,
                        length= ChrSizes,
                        is_circular= rep(FALSE, length(ChrNames)))

##################################################
# Make TranscriptDb from file
print("Making TranscriptDb from GFF3 file")
Db <- makeTranscriptDbFromGFF(file = args[2] ,format = "gff3",chrominfo=chrominfo)

##################################################
## import the MACS output
# Create a GRanges object
print("Importing DGF data.")
bed <- read.table(args[4],sep="\t",  header=FALSE)
col_names<-c("chr","start","end")
names(bed)<-col_names
bedOutput <- bed2GRanges(bed)

# Check which values were out of bounds
# This returns the ends of each GRanges element, and checks it against the length of the corresponding chromosome.
# Creates a vector with the index of offending rows
#x<-which(end(macsOutput) > ChrSizesNamed[as.character(seqnames(macsOutput))])
#remove offending rows from the peak df and repeat construction of GRange object 
#macs_cleaned<-macs[-x,]
#macsOutput <- macs2GRanges(macs_cleaned)
##################################################

# By specifying precidence a peak that is found covering two regions region will be assigned to the first listed
# Therefore in this instance 'exons' refer to coding exons and promoter is upstream of 5'UTR
#The promoter region is specified as 1000bp upstream, likewise for immediately downstream regions.

print("Assigning DGFs to chromosome regions")
c <- assignChromosomeRegion(bedOutput, 
                            proximal.promoter.cutoff=3000L, immediate.downstream.cutoff=1000L,
                            nucleotideLevel=FALSE, precedence=c("fiveUTRs","threeUTRs","Exons","Introns","Promoters"), 
                            TxDb=Db)

##################################################
#Create a dataframe of the results and write to file
print("Writing results to file.")
print(c$percentage)
names(c$percentage)
rownames(c$percentage)
o <-data.frame(Feature=names(c$percentage),Percentage=c$percentage,row.names=NULL)
write.table(o, OUTFILE, sep="\t", col.names=TRUE, row.names=FALSE) 
