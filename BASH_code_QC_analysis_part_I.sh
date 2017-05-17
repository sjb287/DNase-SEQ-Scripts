#!/bin/bash
# BASH_code_QC_analysis_part_I.sh
##################################################
# Description
# Script to check quality of DNase-SEQ data based 
# on the ENCODE consortium pipeline. Calculates
# NSC: Normalized Strand Cross Correlation Coefficient.
# RSC: Relative Strand Cross Correlation Coefficient.
# PCB: PCR bottleneck coefficient.
# SPOT Score: Signal Portion of Tags (SPOT).
# Information about quality metrics, and intepretation can be found: 
# http://genome.ucsc.edu/ENCODE/qualityMetrics.html#definitions
# reference doi:  10.1101/gr.136184.111
#
# Dependencies
# SPP
# bowtie2
# samtools
# UCSC genome browser
# MACS2
##################################################

#Insert script location for SPP
RUN_SPP="{YOUR_PATH}/run_spp.R"

# Get input options
THREADS=$1 # Number of threads to use when running the script
BOWTIE=$2  # Bowtie file for the genome of interest
R1=$3      # Paired End FASTQ file mates 1
R2=$4      # Paired End FASTQ file mates 2
GENOME=$5  # FASTA file for the genome of interest
NAME=$6    # File prefix

# Define Filenames
MAPQ_THRESH=42 # Threshold for alignment quality used to filter reads (adjust accordingly)
SAMFILE="${NAME}.mp.sam" # SAM file containing aligned reads.
TEMPFILE1="${SAMFILE::-4}.bam" # BAM file containing aligned reads.
TEMPFILE2="${SAMFILE::-4}.srt" # Sorted BAM file containing aligned reads (prefix).
TEMPFILE3="$TEMPFILE2.bam" # Sorted BAM file containing aligned reads.
SPP_OUT="${SAMFILE::-4}_SPP_STATS.tab"
FILT_BAM_PREFIX="${NAME}.mp.filt.srt" # Filtered BAM file prefix.
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam" # Filtered BAM file. 
FILT_BAM_INDEX_FILE="${FILT_BAM_PREFIX}.bai" # Index of the filtered BAM file.
FILT_BAM_FILE_MAPSTATS="${FILT_BAM_PREFIX}.flagstat.qc" # QC file.
TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt" # Temporary BAM file prefix.
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam" # Temporary BAM file.
PBC_FILE_QC="${FILT_BAM_PREFIX}.pbc.qc" # Temporary file containing data relating to PCR bottleneck coefficient.
PREFIX_TAGALIGN="${FILT_BAM_PREFIX}.tagAlign" # TAGALIGN file prefix (contains read data https://genome.ucsc.edu/FAQ/FAQformat.html#format15).
TAGALIGN="${PREFIX_TAGALIGN}.gz" # TAGALIGN file.
TAGALIGNPDF="${PREFIX_TAGALIGN}.pdf"
PEAKFILE="${FILT_BAM_PREFIX}.MACS2" # Prefix for MACS2 peak file output (used to calculate SPOT score).
OUTPEAKS="${FILT_BAM_PREFIX}.MACS2_peaks.narrowPeak" # Output file from MACS2.
OUTBED="${OUTPEAKS}.bed" # BED file containing DNase-hypersensitive site locations
T_BEDFILE="${FILT_BAM_PREFIX}.SPOT.bed" # Temporary BED file used in calculation of SPOT score.
SUMMARY="${FILT_BAM_PREFIX}.sum.qc" # Summary text file containing data relating to quality metrics

##################################################
# Align reads
##################################################
echo "Step 1: Aligning reads."
echo `bowtie2 -p $THREADS --local -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 -q -x $BOWTIE -1 $R1 -2 $R2, -S $SAMFILE`

# Create sorted, indexed BAM file and clean up
echo "Step 2: Creating sorted BAM file."
echo `samtools view -@11 -bS $SAMFILE > $TEMPFILE1`
echo `samtools sort -@11 $TEMPFILE1 $TEMPFILE2`

#Remove low quality and unmapped reads
echo "Step 3: Filtering BAM file."
echo `samtools view -F 780 -h -f 2 -q ${MAPQ_THRESH} -@11 -u ${TEMPFILE3} | samtools sort -n - ${TMP_FILT_BAM_PREFIX} -@11`
echo `samtools fixmate -r -O bam ${TMP_FILT_BAM_FILE} - | samtools view -F 780 -f 2 -@11 -u - | samtools sort -@11 - ${FILT_BAM_PREFIX}`

# Index Final BAM file needs position sorted reads
echo "Step 4: Indexing sorted, filtered BAM file."
echo `samtools index ${FILT_BAM_FILE}`
echo `samtools flagstat ${FILT_BAM_FILE} > ${FILT_BAM_FILE_MAPSTATS}`
echo `rm ${TMP_FILT_BAM_FILE}`

##################################################
#Create tagAlignment file
##################################################
echo `samtools view -@11 -b ${FILT_BAM_FILE} | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $TAGALIGN`

##################################################
#Calculate NSC and RSC
##################################################
echo "Step 6: Calculating NSC and RSC values."
echo `Rscript ${RUN_SPP} -c=$TAGALIGN -savp -p=10 -out=$SPP_OUT`

##################################################
# Run MACS 
##################################################
# First Estimate genome size
echo "Step 7: Performing DHS calling using MACS2."
GENOMESIZE=$(faSize ${GENOME} | head -1 | cut -d " " -f1)
echo `macs2 callpeak -t $TAGALIGN -f BED -g $GENOMESIZE -n $PEAKFILE -B -q 0.01 --nomodel --extsize 150 --shift -75 --llocal 50000`

##################################################
# Calculate PCB - PCR bottleneck coefficient (although for DNase-Seq may have biological relevence rather than technical issue) 
##################################################
# Notes on PCB:
# A measure of library complexity, i.e. how skewed the distribution of read counts per location is towards 1 read per location.
# PBC = N1/Nd
# (where N1= number of genomic locations to which EXACTLY one unique mapping read maps, 
# and Nd = the number of genomic locations to which AT LEAST one unique mapping read maps, 
# i.e. the number of non-redundant, unique mapping reads).

# PBC is further described on the ENCODE Software Tools page. Provisionally, 0-0.5 is severe bottlenecking, 
# 0.5-0.8 is moderate bottlenecking, 0.8-0.9 is mild bottlenecking, while 0.9-1.0 is no bottlenecking. Very 
# low values can indicate a technical problem, such as PCR bias, or a biological finding, such as a very rare 
# genomic feature. Nuclease-based assays (DNase, MNase) detecting features with base-pair resolution (transcription 
# factor footprints, positioned nucleosomes) are expected to recover the same read multiple times, resulting in a 
# lower PBC score for these assays. Note that the most complex library, random DNA, would approach 1.0, thus the 
# very highest values can indicate technical problems with libraries. It is the practice for some labs outside of 
# ENCODE to remove redundant reads; after this has been done, the value for this metric is 1.0, and this metric 
# is not meaningful. 82% of TF ChIP, 89% of His ChIP, 77% of DNase, 98% of FAIRE, and 97% of control ENCODE datasets 
# have no or mild bottlenecking.
echo "Step 8: Calculating PCR bottleneck coefficient (PBC)."
echo `samtools sort -@11 -n ${FILT_BAM_FILE} temp`
echo `bedtools bamtobed -bedpe -i temp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'ChrM\ChrC'| sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}`
echo `rm temp.bam`
# Output columns
# TotalReadPairs 
# DistinctReadPairs 
# OneReadPair 
# TwoReadPairs 
# NRF=Distinct/Total 
# PBC1=OnePair/Distinct 
# PBC2=OnePair/TwoPair

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
echo "Step 9: Calculating fraction of reads in peaks (FRiP) score."
echo `awk 'BEGIN{OFS="\t"};{print $1, $2, $3}' $OUTPEAKS > ${OUTBED}`
echo `zcat ${TAGALIGN} | shuf | head -n 5000000 -> ${T_BEDFILE}`
SPT=$(bedtools intersect -f 1E-9 -wa -u -a ${T_BEDFILE} -b ${OUTBED} -bed | wc -l)

##################################################
# Create Summary File
##################################################
echo "Step 11: Creating summary file."
NMAP=`head -1 ${FILT_BAM_FILE_MAPSTATS} | cut -d" " -f1`
PCB=`awk '{print $6}' ${PBC_FILE_QC}`
NSC=`awk '{print $9}' ${SPP_OUT}`
RSC=`awk '{print $10}' ${SPP_OUT}`
QTAG=`awk '{print $11}' ${SPP_OUT}`
PEAKSFDR=`wc -l ${OUTPEAKS} | cut -d" " -f1`
SPOT=$(echo "scale=2; ${SPT}/5000000" | bc)
echo `printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$NMAP" "$PCB" "$NSC" "$RSC" "$QTAG" "$PEAKSFDR" "$SPOT"> ${SUMMARY}`
