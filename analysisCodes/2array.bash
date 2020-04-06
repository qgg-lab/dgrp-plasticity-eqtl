# ====================
# = array processing =
# ====================

# 1. extract probe information from bpmap file
# ============================================================

/home/whuang9/software/R-3.1.2/bin/Rscript extractProbeInfo.R tilingArrayCels/18c/Dm_tiling2_MR_v01.bpmap tilingInt/probe.info > log/extractProbeInfo.Rout 2>&1 &

# 2. prepare sequences for mapping, only consider
#    probes on the BDGP5 assembly
# ============================================================

awk '$1 ~ /^Dm:BDGPv5;chr/' tilingInt/probe.info | perl -wne 'chomp $_; @line = split / /, $_; print "\@$line[0]:$line[1]:$line[3]\n"; print $line[6], "\n"; print "+\n"; print "I" x length($line[6]), "\n";' > tilingInt/probe.seq.fastq

/home/whuang9/software/bwa-0.7.12/bwa samse -r "@RG\tID:dmr2\tPL:affy\tLB:dmr2\tSM:dmr2" /home/whuang9/flybase/fb-r5.57/dmel.bwa0710 <(/home/whuang9/software/bwa-0.7.12/bwa aln -n 0 /home/whuang9/flybase/fb-r5.57/dmel.bwa0710 tilingInt/probe.seq.fastq 2> log/probe.aln.log) tilingInt/probe.seq.fastq > tilingInt/probe.map.sam 2> log/probe.samse.log &

# find good probes
grep -v ^@ tilingInt/probe.map.sam | awk '$5 > 0 {print $1"\t"$3"\t"$4"\t"$6"\t"length($10)"M"}' | sed 's/^Dm:BDGPv5;chr//' | sed 's/^M/dmel_mitochondrion_genome/' | sed 's/:/\t/g' | awk '$1 == $4 && $3+1 == $5 && $6 == $7' | sed 's/^dmel_mitochondrion_genome/M/' | awk '{print "Dm:BDGPv5;chr"$1"_"$2}' > tilingInt/probe.uniq.map

# 3. identify probes that overlap with polymorphisms
# ============================================================

# calculate frequency of alleles for *ALL* variants, no filter on --geno
/home/whuang9/software/plink-v1.90b3j/plink --silent --bfile /home/whuang9/dgrp/freeze2.jgil.hpp.biallelic --keep geno/exp.185line.plink.id --freq --out geno/exp.185line.freq

# calculate non-ref allele frequency
join -t $'\t' <(tail -n+2 /home/whuang9/dgrp/freeze2.jgil.hpp.biallelic.tgeno | cut -d " " -f 3,4 | sed 's/ /\t/' | sort -k1,1) <(tail -n+2 geno/exp.185line.freq.frq | awk '$6 != 0 {print $2"\t"$3"\t"$4"\t"$5}' | sort -k1,1) | awk '{ if ($2 == $4) { print $0 } else { print $1"\t"$2"\t"$3"\t"$4"\t"1-$5 } }' > geno/exp.185line.nonref.freq

# make bed file
awk '$5 >= 0.05' geno/exp.185line.nonref.freq | sed 's/_/\t/g' | awk '{print $1"\t"$2-1"\t"$2-1+length($4)"\t"$1"_"$2"_"$3}' > geno/common.nonref.bed

# intersection
awk -F " " '$1 ~ /^Dm:BDGPv5;/ {print $1"\t"$4"\t"$4 + length($7)"\t"$1"_"$2}' tilingInt/probe.info | sed 's/Dm:BDGPv5;chr//' | sed 's/^M/dmel_mitochondrion_genome/' | /home/whuang9/software/bedtools-2.22.1/bin/bedtools intersect -a stdin -b geno/common.nonref.bed -wa -loj | awk '$5 == "." {print $4}' > tilingInt/probe.novar.list

# 4. make constitutive exon files
# ============================================================

# define constitutive exons
# as bases that appear in all splice isoforms
# use bedtools multiinter
# to count coverage and output intervals
# where coverage == number of isoforms
# ============================================================

# this script runs slow
# if needed, the GTF file can be broken into several parts
# ============================================================

while read line
do
  gene=`echo $line | awk '{print $1}'`
  txs=`echo $line | awk '{print $2}' | sed 's/,/\n/g'`
  txcount=`echo $line | awk '{print $3}'`
  if [[ $txcount -eq 1 ]]
  then
    awk '$3 == "exon"' /home/whuang9/geneExpTemp/cuff/final/final.merge.gtf | awk -F "\"" '$2 == "'$gene'" && $4 == "'$txs'"' | awk '{print $1"\t"$4-1"\t"$5"\t'$gene'"}'
  else
    awk '$3 == "exon"' /home/whuang9/geneExpTemp/cuff/final/final.merge.gtf | awk -F "\"" '$2 == "'$gene'"' | awk '{print $1"\t"$4-1"\t"$5"\t'$gene'"}' | /home/whuang9/software/bedtools-2.22.1/bin/bedtools genomecov -bg -i stdin -g ~/flybase/fb-r5.57/dmel-all-chromosome-r5.57.fasta.fai | awk '$4 == "'$txcount'" {print $1"\t"$2"\t"$3"\t'$gene'"}'
  fi
done < <(awk '$3 == "exon"' /home/whuang9/geneExpTemp/cuff/final/final.merge.gtf | awk -F "\"" '{print $2"\t"$4}' | sort -k1,1 | /home/whuang9/software/bedtools-2.22.1/bin/bedtools groupby -g 1 -c 2,2 -o distinct,count_distinct) > tilingInt/gtf.const.exon.bed

# subtract regions where the exons are covered by multiple genes
sort -k1,1 -k2,2n tilingInt/gtf.const.exon.bed | /home/whuang9/software/bedtools-2.22.1/bin/bedtools genomecov -bg -i stdin -g ~/flybase/fb-r5.57/dmel-all-chromosome-r5.57.fasta.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3}' > tilingInt/constExon.multiCov.bed
/home/whuang9/software/bedtools-2.22.1/bin/bedtools subtract -a tilingInt/gtf.const.exon.bed -b tilingInt/constExon.multiCov.bed | sort -k1,1 -k2,2n > tilingInt/final.merge.uniq.constExon.bed

# 5. find probes that
#    1) fall entirely within genes
#    2) uniquely map, use tilingInt/probe.uniq.map
#    3) no overlap with variants tilingInt/probe.novar.list
# ============================================================

comm -12 <(sort tilingInt/probe.uniq.map) <(sort tilingInt/probe.novar.list) > tilingInt/probe.uniq.novar.list
awk '{print $1"_"$2"\t"$4"\t"$4+length($7)}' tilingInt/probe.info | sort -k1,1 | join -t $'\t' tilingInt/probe.uniq.novar.list - |  awk '{print $0"\t"$1}' | sed 's/^Dm:BDGPv5;chr//' | sed 's/_[0-9]\+//' | sed 's/^M/dmel_mitochondrion_genome/' | /home/whuang9/software/bedtools-2.22.1/bin/bedtools intersect -a stdin -b tilingInt/final.merge.uniq.constExon.bed -wb -f 1 | awk '{print $4"\t"$8}' | sort -k1,1 > tilingInt/probe.gene.assoc
# probes are all uniq (one probe belongs to only one gene)

# 6. extract intensities for the retained probes
# ============================================================

# find probe pmx and pmy to extract
cat <(awk '$1 ~ /^AffxCtrlBkGrGenomic:v1/ {print $1"_"$2"\t"$5"\t"$6"\tbg"}' tilingInt/probe.info) <(sed 's/ /\t/g' tilingInt/probe.info | awk '{print $1"_"$2"\t"$5"\t"$6}' | sort -k1,1 | join -t $'\t' - tilingInt/probe.gene.assoc) > tilingInt/probe.extract.index

# run R script to extract PM intensities
srun -w node1 ~/software/R-3.1.2/bin/Rscript extractArrayInt.R tilingInt/probe.extract.index tilingArrayCels/18c/ tilingInt/18c.int.RData > log/18c.int.extract.Rout 2>&1 &
srun -w node3 ~/software/R-3.1.2/bin/Rscript extractArrayInt.R tilingInt/probe.extract.index tilingArrayCels/25c/ tilingInt/25c.int.RData > log/25c.int.extract.Rout 2>&1 &

# 7. background correction
# ============================================================

srun --cpus-per-task=10 /home/whuang9/software/R-3.1.2/bin/Rscript bgCorrect.R tilingGCRMA.R tilingInt/probe.info tilingInt/18c.int.RData 10 -F tilingInt/18c.female.bg.RData > log/18c.female.bg.Rout 2>&1 &
srun --cpus-per-task=10 /home/whuang9/software/R-3.1.2/bin/Rscript bgCorrect.R tilingGCRMA.R tilingInt/probe.info tilingInt/18c.int.RData 10 -M tilingInt/18c.male.bg.RData > log/18c.male.bg.Rout 2>&1 &
srun --cpus-per-task=10 /home/whuang9/software/R-3.1.2/bin/Rscript bgCorrect.R tilingGCRMA.R tilingInt/probe.info tilingInt/25c.int.RData 10 -F tilingInt/25c.female.bg.RData > log/25c.female.bg.Rout 2>&1 &
srun --cpus-per-task=10 /home/whuang9/software/R-3.1.2/bin/Rscript bgCorrect.R tilingGCRMA.R tilingInt/probe.info tilingInt/25c.int.RData 10 -M tilingInt/25c.male.bg.RData > log/25c.male.bg.Rout 2>&1 &

# 8. preliminary normalization
# ============================================================

srun -w node1 ~/software/R-3.1.2/bin/Rscript prelimNorm.R > log/prelimNorm.Rout 2>&1 &

# 9. summarize results
# ============================================================

~/software/R-3.1.2/bin/Rscript prelimSummary.R > log/prelimSummary.Rout 2>&1 &

# extract header information for batch order
~/software/R-3.1.2/bin/Rscript -e 'load("tilingInt/18c.int.RData"); cel.headers.18c <- cel.headers; load("tilingInt/25c.int.RData"); cel.headers.25c <- cel.headers; save(cel.headers.18c, cel.headers.25c, file = "tilingInt/cel.headers.RData");' > log/extract.cel.headers.Rout 2>&1 &

# 10. final normalization
# cite the first tiling array paper and Zhou et al 2012 to 
# justify normalization within temp but for sexes separately
# temperature does not change expression that much at all
# this is also apparent in the RNA-Seq data
# after all, they are normal looking flies under 
# either cold or warm temperatures
# ============================================================

# identify files to remove
~/software/R-3.1.2/bin/Rscript -e 'load("tilingInt/prelim.norm.summary.RData"); write.table(c(colnames(female.18c.quantile)[which(female.18c.quantile[1, ] < -6 | female.18c.quantile[7, ] > 6)], colnames(female.25c.quantile)[which(female.25c.quantile[1, ] < -6 | female.25c.quantile[7, ] > 6)], colnames(male.18c.quantile)[which(male.18c.quantile[1, ] < -6 | male.18c.quantile[7, ] > 6)], colnames(male.25c.quantile)[which(male.25c.quantile[1, ] < -6 | male.25c.quantile[7, ] > 6)]), file = "tilingInt/bad.array.list", sep = " ", col.names = F, row.names = F, quote = F);'

srun ~/software/R-3.1.2/bin/Rscript normGeneExp.R > log/normGeneExp.Rout 2>&1 &

# 10. prepare file with gene information
# ============================================================

perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[2] ne "exon" || $line[1] eq "miRNA" || $line[8] =~ m/class_code \"j\"/) { next; }; if ($line[8] =~ m/gene_id \"(.*?)\"/) { $id = $1; }; $symbol = "-"; $type = "-"; if (!($id =~ /XLOC/)) { if ($line[8] =~ m/gene_name \"(.*?)\"/) { $symbol = $1; }; if ($line[8] =~ m/gene_biotype \"(.*?)\"/) { $type = $1; } } else { $type = "NTR"; }; print $id, "\t", $line[0], "\t", $line[3], "\t", $line[4], "\t", $symbol, "\t", $type, "\n"; ' /home/whuang9/geneExpTemp/cuff/final/final.merge.gtf | sort -k1,1 | ~/software/bedtools-2.22.1/bin/bedtools groupby -g 1 -c 2,3,4,5,6 -o distinct,min,max,distinct,distinct | awk '{print $1"\t"$5"\t"$6"\t"$2":"$3"-"$4}' > cuff/gene.info

# 11. the arrays have scan date, this is a natural batch variable
#   	and the scaled expression plot does show batch effect
#     perform a pca on normal quantile transformed data
#     first on two temp combined, then within temp
# ============================================================

Rscript geneExpPCA.R tilingInt/female.gene.exp.RData female.gene.exp tilingInt/cel.headers.RData tilingInt/female.pca.RData > log/female.pca.Rout 2>&1 &
Rscript geneExpPCA.R tilingInt/male.gene.exp.RData male.gene.exp tilingInt/cel.headers.RData tilingInt/male.pca.RData > log/male.pca.Rout 2>&1 &

# 12. sva adjustment
# ============================================================

Rscript sva.R tilingInt/female.gene.exp.RData female.gene.exp tilingInt/cel.headers.RData ~/dgrp/adjustData.RData tilingInt/female.sva.adjust.RData tilingInt/female.gene.exp.adjust.RData > log/female.sva.Rout 2>&1 &
Rscript sva.R tilingInt/male.gene.exp.RData male.gene.exp tilingInt/cel.headers.RData ~/dgrp/adjustData.RData tilingInt/male.sva.adjust.RData tilingInt/male.gene.exp.adjust.RData > log/male.sva.Rout 2>&1 &
