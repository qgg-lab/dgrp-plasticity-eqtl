# ================================================
# = script to merge assemblies across conditions =
# ================================================

# run in hyperion:/home/whuang9/geneExpTemp/
# on the login node, since it's light
# don't run this as a batch script, run it interactively
# except the last step cuffquant
# ============================================================

# merge assemblies from all cufflinks 
# the strategy is to work on only class_code "i", "j", and "u".

# ideally cuffmerge should serve this purpose quite well.
# However, it filters out transcripts that are present in the cufflinks assemblies
# which appears to be based on **expression level** in addition to others
# e.g. Many intronic isoforms are filtered by cuffmerge

# By investigating the cuffmerge utility, it's found that it's actually a serieas of 
# commands utilizing cufflinks and cuffcompare
# The filtering may be introduced by the cuffmerge (really cufflinks on GTF converted SAM) 
# utility with default options
# ============================================================

# 1. run cuffcompare
# ============================================================
mkdir cuff/comp
mkdir cuff/merge
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cuffcompare -o cuff/comp/comp -r /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf -s /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa cuff/female18c/transcripts.gtf cuff/female25c/transcripts.gtf cuff/male18c/transcripts.gtf cuff/male25c/transcripts.gtf > log/cuffcompare.log 2>&1 &

# 2. regardless of what cuffcompare is able to merge
#    identify transcripts whose class_code are "i", "j", and "u".
# ============================================================
cat cuff/female18c/comp.transcripts.gtf.tmap cuff/female25c/comp.transcripts.gtf.tmap cuff/male18c/comp.transcripts.gtf.tmap cuff/male25c/comp.transcripts.gtf.tmap | awk '$3 == "i" || $3 == "j" || $3 == "u" {print $5}' | sort > cuff/merge/merge.candidate.iju.id

# apply a filter for novel exons that are too short (< 8bp to be consistent with )
perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[8] =~ m/gene_id \"(.*?)\";.+transcript_id \"(.*?)\";.+oId \"(.*?)\";.+class_code \"(.*?)\";/) { print join("-", @line[(0,3,4)]), "\t", $1, "\t", $2, "\t", $3, "\n"; }' cuff/comp/comp.combined.gtf | sort -k1,1 | ~/software/bedtools-2.22.1/bin/bedtools groupby -g 1 -c 4 -o collapse | perl -wne 'chomp $_; ($exon, $id) = split /\t/, $_; @exoninfo = split /-/, $exon; $exonlen = $exoninfo[2] - $exoninfo[1] + 1; if ($exonlen < 8 && !($id =~ m/FBtr/)) { print $_, "\n"; }' | cut -f 2 | sed 's/,/\n/g' | sort | uniq | comm -13 - cuff/merge/merge.candidate.iju.id > cuff/merge/merge.candidate.id

# 3. Make SAM files and bam files
# ============================================================
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/gtf_to_sam -F cuff/female18c/transcripts.gtf cuff/merge/female18c.transcript.sam
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/gtf_to_sam -F cuff/female25c/transcripts.gtf cuff/merge/female25c.transcript.sam
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/gtf_to_sam -F cuff/male18c/transcripts.gtf cuff/merge/male18c.transcript.sam
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/gtf_to_sam -F cuff/male25c/transcripts.gtf cuff/merge/male25c.transcript.sam

cat <(/home/whuang9/software/samtools-1.2/samtools view -H map/female18c/accepted_hits.bam | grep -v ^@PG | grep -v ^@RG) <(cat cuff/merge/female18c.transcript.sam cuff/merge/female25c.transcript.sam cuff/merge/male18c.transcript.sam cuff/merge/male25c.transcript.sam | sort -k1,1 | join -t $'\t' cuff/merge/merge.candidate.id -) | /home/whuang9/software/samtools-1.2/samtools view -b - | /home/whuang9/software/samtools-1.2/samtools sort - cuff/merge/merge.candidate

# 4. Run cufflinks for the merged transcripts
#    these should not overlap with any reference transcripts
#    so no need to supply a reference GTF
# ============================================================
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cufflinks -o cuff/merge/ --label merge --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --library-type transfrags --total-hits-norm --max-bundle-frags 1000000 --min-isoform-fraction 0 --pre-mrna-fraction 0 --max-intron-length 25000 --junc-alpha 1 --small-anchor-fraction 0 --min-frags-per-transfrag 0 --overhang-tolerance 8 --max-bundle-length 400000 --min-intron-length 20 --trim-3-avgcov-thresh 0 --trim-3-dropoff-frac 0 --max-multiread-fraction 1 --overlap-radius 50 --3-overhang-tolerance 200 --intron-overhang-tolerance 50 --no-update-check --no-5-extend cuff/merge/merge.candidate.bam > log/merge.candidate.log 2>&1

# 5. Run cuffcompare
# ============================================================
/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cuffcompare -o cuff/merge/comp -r /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf -s /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa cuff/merge/transcripts.gtf > log/compare.merge.log 2>&1

# 6. Extract candidate transcripts
# ============================================================
awk '$3 == "exon" {print $0" class_code \"=\";"}' /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf > cuff/merge/fb-r5.57.plus.candidate.gtf
perl extractGTF.pl cuff/merge/comp.combined.gtf /home/whuang9/flybase/fb-r5.57/gene.symbol >> cuff/merge/fb-r5.57.plus.candidate.gtf

# 7. Run cuffquant to re-estimate gene expression
#    cuffquant can only be run on one sample at a time
# ============================================================
mkdir cuff/quant

mkdir cuff/quant/female18c
mkdir cuff/quant/female25c
mkdir cuff/quant/male18c
mkdir cuff/quant/male25c

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/quant/female18c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/merge/fb-r5.57.plus.candidate.gtf map/female18c/accepted_hits.bam > log/cuff.female18c.quant.log 2>&1 &

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/quant/female25c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/merge/fb-r5.57.plus.candidate.gtf map/female25c/accepted_hits.bam > log/cuff.female25c.quant.log 2>&1 &

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/quant/male18c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/merge/fb-r5.57.plus.candidate.gtf map/male18c/accepted_hits.bam > log/cuff.male18c.quant.log 2>&1 &

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/quant/male25c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/merge/fb-r5.57.plus.candidate.gtf map/male25c/accepted_hits.bam > log/cuff.male25c.quant.log 2>&1 &

# 8. normalize gene expression
# ============================================================
srun --cpus-per-task=4 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffnorm --output-dir cuff/quant/ --labels f18c,f25c,m18c,m25c --num-threads 4 --library-type fr-unstranded --library-norm-method geometric --output-format simple-table --verbose --no-update-check cuff/merge/fb-r5.57.plus.candidate.gtf cuff/quant/female18c/abundances.cxb cuff/quant/female25c/abundances.cxb cuff/quant/male18c/abundances.cxb cuff/quant/male25c/abundances.cxb > log/cuff.quant.log 2>&1 &

# 9. estimate isoform fraction
#    only retain isoforms whose expression >= 0.05 of 
#    total gene expression in at least one of the four samples
# ============================================================
awk -F "\"" '{print $4"\t"$2}' cuff/merge/fb-r5.57.plus.candidate.gtf | sort -k1,1 | uniq | join -t $'\t' - <(tail -n+2 cuff/quant/isoforms.fpkm_table | sort -k1,1) | sort -k2,2 | ~/software/bedtools-2.22.1/bin/bedtools groupby -g 2 -c 1,3,4,5,6 -o collapse,collapse,collapse,collapse,collapse | perl -wne 'use List::Util qw(sum); chomp $_; @line = split /\t/, $_; @txs = split /,/, $line[1]; for (my $i = 2; $i <= 5; $i++) { @fpkm = split /,/, $line[$i]; $sum = sum(@fpkm) + 1e-20; @ratio = map { $_/$sum } @fpkm; for (my $j = 0; $j <= $#txs; $j++) { $txs[$j] .= ",".$ratio[$j]; } }; for (my $k = 0; $k <= $#txs; $k++) { @txratio = split /,/, $txs[$k]; print $line[0], "\t", join("\t", @txratio), "\n"; }' | awk '$2 ~ /TCONS/ && ($3 >= 0.05 || $4 >= 0.05 || $5 >= 0.05 || $6 >= 0.05) { print $2 }' > cuff/merge/novel.candidate.id

# 10. obtain final assembly
# ============================================================
mkdir cuff/final
awk '$3 == "exon" {print $0" class_code \"=\";"}' /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf > cuff/final/final.merge.gtf
cat cuff/merge/fb-r5.57.plus.candidate.gtf | perl -we '%txs = (); open TX, "<cuff/merge/novel.candidate.id"; while (<TX>) { chomp $_; $txs{$_} = 1; } close TX; while (<>) { if ($_ =~ m/transcript_id \"(.*?)\";/) { if (defined($txs{$1})) { print $_; } } }' >> cuff/final/final.merge.gtf

# 11. one last expression + normalization
# ============================================================
mkdir cuff/final/female18c
mkdir cuff/final/female25c
mkdir cuff/final/male18c
mkdir cuff/final/male25c

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/final/female18c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/final/final.merge.gtf map/female18c/accepted_hits.bam > log/cuff.final.female18c.quant.log 2>&1 &

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/final/female25c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/final/final.merge.gtf map/female25c/accepted_hits.bam > log/cuff.final.female25c.quant.log 2>&1 &

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/final/male18c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/final/final.merge.gtf map/male18c/accepted_hits.bam > log/cuff.final.male18c.quant.log 2>&1 &

srun --cpus-per-task=8 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffquant --output-dir cuff/final/male25c/ --num-threads 8 --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --max-bundle-frags 1000000 --no-update-check --verbose cuff/final/final.merge.gtf map/male25c/accepted_hits.bam > log/cuff.final.male25c.quant.log 2>&1 &

srun --cpus-per-task=4 ~/software/cufflinks-2.2.1.Linux_x86_64/cuffnorm --output-dir cuff/final/ --labels f18c,f25c,m18c,m25c --num-threads 4 --library-type fr-unstranded --library-norm-method geometric --output-format simple-table --verbose --no-update-check cuff/final/final.merge.gtf cuff/final/female18c/abundances.cxb cuff/final/female25c/abundances.cxb cuff/final/male18c/abundances.cxb cuff/final/male25c/abundances.cxb > log/cuff.final.quant.log 2>&1 &
