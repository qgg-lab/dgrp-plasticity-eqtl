#! /bin/bash
# ===============================================
# = script to summarize RNA-seq mapping results =
# ===============================================

# this script only works with the directory
# structure output by TopHat
# including
# 1. number of input reads
# 2. number of mapped reads (%)
# 3. number of uniquely mapped reads (%)
# 4. number of reads overlapped known exons entirely (%)
# 5. number of reads mapped to known + novel junctions (%)
# 6. number of reads that overlapped with known introns (%)
# ============================================================

# command line arguments
# ============================================================
args=`getopt -o "p,b:,g:,e:,s:,o:" -l "paired,bam:,gtf:,bed:,sam:,out:" -- "$@"`
echo "command line arguments given: $args"
eval set -- "$args"

# parse arguments
# ============================================================

# default
paired=0
out=out

while true;
do
  case "$1" in
    
    -p|--paired)  paired=1; shift;;
    
    -b|--bam)     bam=$2; shift 2;;
    
    -g|--gtf)     gtf=$2; shift 2;;
    
    -e|--bed)     bed=$2; shift 2;;
    
    -s|--sam)     sam=$2; shift 2;;
    
    -o|--out)     out=$2; shift 2;;
    
    --)           shift; break;;
    
  esac
done

# check if files and programs exist
# ============================================================

# bam file and prep_reads.info
if [[ -e $bam/prep_reads.info ]]
then
  echo "use prep_read information from $bam/prep_reads.info."
else
  echo "cannot find $bam/prep_reads.info."
  exit 1
fi

if [[ -e $bam/accepted_hits.bam ]]
then
  echo "use $bam/accepted_hits.bam."
else
  echo "cannot $bam/accepted_hits.bam."
  exit 1
fi

if [[ -e $bam/unmapped.bam ]]
then
  echo "use $bam/unmapped.bam."
else
  echo "cannot $bam/unmapped.bam."
  exit 1
fi

# gtf file
if [[ -e $gtf ]]
then
  echo "use $gtf."
else
  echo "cannot find $gtf."
  exit 1
fi

# bedtools
bedcheck=`command -v $bed`
if [[ $bedcheck ]]
then
  echo "using $bedcheck."
else
  echo "cannot find bedtools."
  exit 1
fi

# samtools
samcheck=`command -v $sam`
if [[ $samcheck ]]
then
  echo "using $samcheck."
else
  echo "cannot find samtools."
  exit 1
fi

# check if $out is given
if [[ $out ]]
then
  echo "write output to $out."
else
  echo "output file not given."
  exit 1
fi

# extract number of input reads
# for paired end reads, this is the number of fragments
# for single end reads, this is simply the number of reads
# ============================================================
nReads=`grep reads_in $bam/prep_reads.info | awk -F "=" '{print $2}' | head -n 1`

if [[ $paired -eq 1 ]]
then
  nFrags=$(($nReads + $nReads))
  echo -n "total number of sequenced reads: " > $out
  printf "%'i" $nReads >> $out
  echo " x 2" >> $out
else
  nFrags=$nReads
  echo -n "total number of sequenced reads: " > $out
  printf "%'i" $nReads >> $out
  echo "" >> $out
fi

# summarize number of total reads mapped
# ============================================================

nUnMap=`$sam view $bam/unmapped.bam | awk '{ print $1":"and($2, 0x40) }' | uniq | sort | uniq | wc -l`
echo -n "number of mapped reads (%): " >> $out
nMap=$(($nFrags - $nUnMap))
printf "%'i" $nMap >> $out
perl -e '($total, $map) = @ARGV; printf " (%.2f\%)\n", $map/$total*100;' $nFrags $nMap >> $out

# uniquely mapped reads
# ============================================================
nUniq=`$sam view $bam/accepted_hits.bam | grep NH:i:1 | cut -f 1 | wc -l`
echo -n "number of uniquely mapped reads (%): " >> $out
printf "%'i" $nUniq >> $out
perl -e '($total, $uniq) = @ARGV; printf " (%.2f\%)\n", $uniq/$total*100;' $nFrags $nUniq >> $out

# reads overlapped with known exons
# ============================================================
echo -n "number of reads overlapped with known exons (%): " >> $out
nExon=`$bed intersect -abam $bam/accepted_hits.bam -b $gtf | $sam view -  | awk '{ print $1":"and($2, 0x40) }' | uniq | sort | uniq | wc -l`
printf "%'i" $nExon >> $out
perl -e '($total, $exon) = @ARGV; printf " (%.2f\%)\n", $exon/$total*100;' $nFrags $nExon >> $out

# reads mapped to splice junctions
# ============================================================
echo -n "number of junction reads (%): " >> $out
nJunc=`$sam view $bam/accepted_hits.bam | awk '$6 ~ /N/ { print $1":"and($2, 0x40) }' | sort | uniq | wc -l`
printf "%'i" $nJunc >> $out
perl -e '($total, $junc) = @ARGV; printf " (%.2f\%)\n", $junc/$total*100;' $nFrags $nJunc >> $out

# reads that mapped entirely to introns or intergenic regions
# this can be obtained by subtraction
# ============================================================
echo -n "number of reads mapped entirely within introns or intergenic regions (%): " >> $out
nIntron=$(($nMap - $nExon))
printf "%'i" $nIntron >> $out
perl -e '($total, $intron) = @ARGV; printf " (%.2f\%)\n", $intron/$total*100;' $nFrags $nIntron >> $out
