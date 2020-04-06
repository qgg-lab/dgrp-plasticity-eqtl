#! /bin/bash
# run in hyperion:/home/whuang9/geneExpTemp/
# ============================================================

# job array to run eQTL mapping
# ============================================================

# request resource and sbatch options
# ============================================================

#SBATCH --job-name=eqtl.male.18c
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/whuang9/geneExpTemp/log/eqtl.male.18c.JobArray%a.log
#SBATCH --error=/home/whuang9/geneExpTemp/log/eqtl.male.18c.JobArray%a.log

/home/whuang9/software/plink-v1.90b3w/plink --silent --bfile geno/dgrp.common --pheno eqtl/permPheno/male.18c.perm"${SLURM_ARRAY_TASK_ID}".pheno --all-pheno --assoc --pfilter 0.00001 --out eqtl/perm/male.18c.perm"${SLURM_ARRAY_TASK_ID}"/eqtl.perm
