#! /bin/bash
# run in hyperion:/home/whuang9/geneExpTemp/
# ============================================================

# job array to simulate expression and correlation
# ============================================================

# request resource and sbatch options
# ============================================================

#SBATCH --job-name=male.sim.corr.18c
#SBATCH --array=1-1000
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/whuang9/geneExpTemp/log/male.sim.corr.18c.JobArray%a.log
#SBATCH --error=/home/whuang9/geneExpTemp/log/male.sim.corr.18c.JobArray%a.log

/home/whuang9/software/R-3.2.2/bin/Rscript simCorr18C.R qg/male.sig.gene.exp.RData "${SLURM_ARRAY_TASK_ID}" simCorr/male.sim.corr.18c.perm"${SLURM_ARRAY_TASK_ID}".RData > /home/whuang9/geneExpTemp/log/male.sim.corr.18c.JobArray"${SLURM_ARRAY_TASK_ID}".log 2>&1
