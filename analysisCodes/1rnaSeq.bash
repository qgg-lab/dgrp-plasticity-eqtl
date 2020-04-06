# ==================================
# = rna-seq to improve annotations =
# ==================================

# 1. tophat
# ============================================================

sbatch rnaseqMap18C.sbatch
sbatch rnaseqMap25C.sbatch

# 2. cufflinks, in sex/temperature separately
# ============================================================

sbatch cufflinks.sbatch

# 3. merge GTF assemblies
# ============================================================

# run commands in the script mergeAssembly.bash
# which produces the final GTF cuff/final/final.merge.gtf
# and abundance estimates

# 4. summarize RNA-seq
# ============================================================

sbatch runSummarize.sbatch
