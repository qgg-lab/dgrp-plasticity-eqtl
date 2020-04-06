#! /bin/bash
# run in hyperion:/home/whuang9/geneExpTemp/
# ============================================================

# summarize run results
# ============================================================

# request resource
# ============================================================

#SBATCH --cpus-per-task=16
#SBATCH --output=/home/whuang9/geneExpTemp/log/summarizeMap.log

# ============================================================

# summarize mapping statistics
# ============================================================

bash summarizeMap.bash --paired --bed /home/whuang9/software/bedtools-2.22.1/bin/bedtools --bam map/female18c/ --out map/female18c.summary --sam /home/whuang9/software/samtools-1.2/samtools --gtf /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf > log/female18c.summary.log
bash summarizeMap.bash --paired --bed /home/whuang9/software/bedtools-2.22.1/bin/bedtools --bam map/female25c/ --out map/female25c.summary --sam /home/whuang9/software/samtools-1.2/samtools --gtf /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf > log/female25c.summary.log
bash summarizeMap.bash --paired --bed /home/whuang9/software/bedtools-2.22.1/bin/bedtools --bam map/male18c/ --out map/male18c.summary --sam /home/whuang9/software/samtools-1.2/samtools --gtf /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf > log/male18c.summary.log
bash summarizeMap.bash --paired --bed /home/whuang9/software/bedtools-2.22.1/bin/bedtools --bam map/male25c/ --out map/male25c.summary --sam /home/whuang9/software/samtools-1.2/samtools --gtf /home/whuang9/flybase/fb-r5.57/fb-r5.57.clean.gtf > log/male25c.summary.log
