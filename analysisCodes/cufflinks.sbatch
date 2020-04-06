#! /bin/bash
# run in hyperion:/home/whuang9/geneExpTemp/
# ============================================================

# run cufflinks for each sample separately
# this serves as the initial cufflinks to obtain
# GTF files, which will be subject to merge later
# mapped RNA-seq bam files are located
# map/female18c/ etc.
# ============================================================

# request resource
# ============================================================

#SBATCH --cpus-per-task=8
#SBATCH --output=/home/whuang9/geneExpTemp/log/cufflinks.log

# run cufflinks
# ============================================================

/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cufflinks --output-dir cuff/female18c/ --label f18c --num-threads 8 --GTF-guide /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --total-hits-norm --max-bundle-frags 1000000 --min-isoform-fraction 0.1 --pre-mrna-fraction 0.05 --max-intron-length 25000 --junc-alpha 0.001 --small-anchor-fraction 0.08 --min-frags-per-transfrag 10 --overhang-tolerance 8 --max-bundle-length 400000 --min-intron-length 20 --trim-3-avgcov-thresh 10 --trim-3-dropoff-frac 0.1 --max-multiread-fraction 0.50 --overlap-radius 50 --3-overhang-tolerance 200 --intron-overhang-tolerance 50 --no-update-check map/female18c/accepted_hits.bam

/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cufflinks --output-dir cuff/male18c/ --label m18c --num-threads 8 --GTF-guide /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --total-hits-norm --max-bundle-frags 1000000 --min-isoform-fraction 0.1 --pre-mrna-fraction 0.05 --max-intron-length 25000 --junc-alpha 0.001 --small-anchor-fraction 0.08 --min-frags-per-transfrag 10 --overhang-tolerance 8 --max-bundle-length 400000 --min-intron-length 20 --trim-3-avgcov-thresh 10 --trim-3-dropoff-frac 0.1 --max-multiread-fraction 0.50 --overlap-radius 50 --3-overhang-tolerance 200 --intron-overhang-tolerance 50 --no-update-check map/male18c/accepted_hits.bam

/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cufflinks --output-dir cuff/female25c/ --label f25c --num-threads 8 --GTF-guide /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --total-hits-norm --max-bundle-frags 1000000 --min-isoform-fraction 0.1 --pre-mrna-fraction 0.05 --max-intron-length 25000 --junc-alpha 0.001 --small-anchor-fraction 0.08 --min-frags-per-transfrag 10 --overhang-tolerance 8 --max-bundle-length 400000 --min-intron-length 20 --trim-3-avgcov-thresh 10 --trim-3-dropoff-frac 0.1 --max-multiread-fraction 0.50 --overlap-radius 50 --3-overhang-tolerance 200 --intron-overhang-tolerance 50 --no-update-check map/female25c/accepted_hits.bam

/home/whuang9/software/cufflinks-2.2.1.Linux_x86_64/cufflinks --output-dir cuff/male25c/ --label m25c --num-threads 8 --GTF-guide /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf --mask-file /home/whuang9/flybase/fb-r5.57/fb-r5.57.rRNA.mito.gtf --frag-bias-correct /home/whuang9/flybase/fb-r5.57/dmel.rmsk.fa --multi-read-correct --library-type fr-unstranded --total-hits-norm --max-bundle-frags 1000000 --min-isoform-fraction 0.1 --pre-mrna-fraction 0.05 --max-intron-length 25000 --junc-alpha 0.001 --small-anchor-fraction 0.08 --min-frags-per-transfrag 10 --overhang-tolerance 8 --max-bundle-length 400000 --min-intron-length 20 --trim-3-avgcov-thresh 10 --trim-3-dropoff-frac 0.1 --max-multiread-fraction 0.50 --overlap-radius 50 --3-overhang-tolerance 200 --intron-overhang-tolerance 50 --no-update-check map/male25c/accepted_hits.bam
