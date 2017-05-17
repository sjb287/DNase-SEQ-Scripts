#QC Analysis Part II

NAME=$1
OUTBED=$2
SUMMARY=$3
PEAK=$4
SAMFILE="${NAME}.mp.sam"
FILT_BAM_PREFIX="${NAME}.mp.filt.srt"
FILT_BAM_FILE_MAPSTATS="${FILT_BAM_PREFIX}.flagstat.qc" # QC file
SPP_OUT="${SAMFILE::-4}_SPP_STATS.tab"
PBC_FILE_QC="${FILT_BAM_PREFIX}.pbc.qc"
PREFIX_TAGALIGN="${FILT_BAM_PREFIX}.tagAlign"
TAGALIGN="${PREFIX_TAGALIGN}.gz"
T_BEDFILE="${FILT_BAM_PREFIX}.SPOT.bed"


##################################################
# Caluclating (Fraction of reads in peaks) FRiP score
##################################################
# See http://genome.cshlp.org/content/22/9/1813
# For point-source data sets, we calculate the fraction of all mapped reads 
# that fall into peak regions identified by a peak-calling algorithm
# We used bedtools to map the interesection of reads and peaks

# Parameters
# -f    = minimum overlap (1E-9 = 1bp)
# -wa   = write original entry A for each overlap
# -u    = write original entry A once if any overlap (i.e. just report at least one overlap
# -abam = BAM file A, in this instance it is reads
# -b    = Bed file to overlap (in this instance it is the peaks called) 

# Convert Peaks file to bed
# Take a sample of 5M reads from sample tagalign file 
echo `zcat ${TAGALIGN} | shuf | head -n 5000000 -> ${T_BEDFILE}`
SPT=$(bedtools intersect -f 1E-9 -wa -u -a ${T_BEDFILE} -b ${OUTBED} -bed | wc -l)

##################################################

NMAP=`head -1 ${FILT_BAM_FILE_MAPSTATS} | cut -d" " -f1`
PCB=`awk '{print $6}' ${PBC_FILE_QC}`
NSC=`awk '{print $9}' ${SPP_OUT}`
RSC=`awk '{print $10}' ${SPP_OUT}`
QTAG=`awk '{print $11}' ${SPP_OUT}`
PEAKSFDR=`wc -l ${OUTBED} | cut -d" " -f1`
SPOT=$(echo "scale=2; ${SPT}/5000000" | bc)
echo `printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$NMAP" "$PCB" "$NSC" "$RSC" "$QTAG" "$PEAKSFDR" "$SPOT">> ${SUMMARY}`
exit
