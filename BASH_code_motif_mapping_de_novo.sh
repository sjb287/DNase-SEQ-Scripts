#script for Motif Analysis
#Usage ./BASH_code_motif_mapping.sh [GENOME] [DGFs] [Prefix] [GFF3_FILE] [CHR_FILE] [MOTIF_FILE]

##################################################
#Input Options
##################################################
#FASTA file of the genome of interest 
GENOME=$1
#BED file containing all DGF locations in the genome
BED_FILE=$2
#A Prefix to use for output files e.g. Sb for Sorghum Bicolor
PREFIX="${3}_DGF_fdr_0.01"
#GFF3 file
GFF3=$4
#Chromosome file
CHR=$5
#Motif file in MEME format
MOTIF_FILE=$6

#################################################
#File and Directory Names
##################################################
#Output FASTA file containing sequences of all DGFs
FASTA_FILE="/home/sjb287/documents/DNase/Files_Publication/Files_FASTA/${PREFIX}.fa"
#Output file containing the base pair frequencies as a background file for FIMO 
bg_FILE="/home/sjb287/documents/DNase/Files_Publication/Files_Background/${PREFIX}.meme.bg"
#Prefix for output files containing motif hits
OUTFILE="${PREFIX}_fimo_hits"
#Ouput text file for all motif hits
OUTFILE_ALL="${OUTFILE}_all.txt"
#Ouput text file renaming TF IDs
OUTFILE_NAMES="${OUTFILE}_all_names.txt"
#Output bed file for all motif hits
OUTFILE_ALL_BED="${OUTFILE}_all.bed"

##################################################
#Shell Commands
##################################################
#Step 1: Extract FASTA Sequences Using bedtools.
echo 'Extracting FASTA Sequences.'
echo `bedtools getfasta -fi ${GENOME} -bed ${BED_FILE} -fo ${FASTA_FILE}`

#Step 2: Create Background File for FIMO analysis Using fasta-get-markov from the MEME Suite.
echo 'Creating Background File for FIMO.'
echo `fasta-get-markov ${FASTA_FILE} ${bg_FILE}`

#Step 3: Split Input Motif File Containing All Motifs into Individual Motif Files.
echo 'Splitting Input Motifs File for Analysis.'
echo `csplit -f 'Motif_' -s ${MOTIF_FILE} /MEME\ version/ "{*}"`

#Step 4: Run Motif Mapping Using FIMO.
#start a counter to name output files.
echo `i=0`
echo 'Running Motif Mapping Using FIMO.'
echo `for f in ./*; do fimo --bgfile ${bg_FILE} --o "${OUTFILE}_${i}.txt" $f ${FASTA_FILE}; let i=i+1;done`
echo 'Joining Mapping Files.'
echo `find . -name \*fimo.txt -exec cat {} >> temp.txt \;`
echo `grep -v '#' temp.txt > ${OUTFILE_ALL}`

#Step 5: Annotate Hits
#Step 6: Convert Output .txt File to .bed File
echo 'Converting Mapping Files to BED Format.'
echo `awk '{OFS="\t"}{split($2,a,":");split(a[2],b,"-");print a[1],b[1],b[2],$1}' ${OUTFILE_ALL} > ${OUTFILE_ALL_BED}`

#Step 7: Cleaning Up.
echo 'Cleaning Up.'
echo `rm Motif*`
echo 'Analysis Complete.'
