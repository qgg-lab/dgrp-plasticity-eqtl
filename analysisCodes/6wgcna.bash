# =====================================
# = co-expression network using wgcna =
# =====================================

# 1. get expression data
# ============================================================

~/software/R-3.2.2/bin/Rscript -e 'load("qg/female.adj.qgSingleTemp.RData"); sig.18c.gene <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.25c.gene <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); sig.gene <- sort(unique(c(sig.18c.gene, sig.25c.gene))); geno.line.order <- read.table("geno/dgrp.common.fam", as.is = TRUE)[, 1]; load("qg/female.adj.qgBLUP.RData"); blup.18c.data <- blup.18c[sig.gene, match(geno.line.order, line.order)]; blup.25c.data <- blup.25c[sig.gene, match(geno.line.order, line.order)]; sig.18c.gene <- gene.name[sig.18c.gene]; sig.25c.gene <- gene.name[sig.25c.gene]; sig.gene.name <- gene.name[sig.gene]; save(blup.18c.data, blup.25c.data, sig.18c.gene, sig.25c.gene, sig.gene.name, file = "qg/female.sig.gene.exp.RData")' > log/female.sig.gene.exp.Rout 2>&1 &

~/software/R-3.2.2/bin/Rscript -e 'load("qg/male.adj.qgSingleTemp.RData"); sig.18c.gene <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.25c.gene <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); sig.gene <- sort(unique(c(sig.18c.gene, sig.25c.gene))); geno.line.order <- read.table("geno/dgrp.common.fam", as.is = TRUE)[, 1]; load("qg/male.adj.qgBLUP.RData"); blup.18c.data <- blup.18c[sig.gene, match(geno.line.order, line.order)]; blup.25c.data <- blup.25c[sig.gene, match(geno.line.order, line.order)]; sig.18c.gene <- gene.name[sig.18c.gene]; sig.25c.gene <- gene.name[sig.25c.gene]; sig.gene.name <- gene.name[sig.gene]; save(blup.18c.data, blup.25c.data, sig.18c.gene, sig.25c.gene, sig.gene.name, file = "qg/male.sig.gene.exp.RData")' > log/male.sig.gene.exp.Rout 2>&1 &

# 2. simulation
# ============================================================

sbatch -x node2 femaleSimCorr18CJobArray.sbatch
sbatch -x node2 maleSimCorr18CJobArray.sbatch 

# 3. summarize simulation results
# ============================================================

srun -w node3 ~/software/R-3.2.2/bin/Rscript summarizeSimCorr18C.R qg/female.sig.gene.exp.RData qg/female.adj.qgGxE.RData 10 simCorr/female.sim.corr.18c.perm 1000 simCorr/female.simCorr18C.RData > log/female.simCorr18C.Rout 2>&1 &
srun -w node4 ~/software/R-3.2.2/bin/Rscript summarizeSimCorr18C.R qg/male.sig.gene.exp.RData qg/male.adj.qgGxE.RData 10 simCorr/male.sim.corr.18c.perm 1000 simCorr/male.simCorr18C.RData > log/male.simCorr18C.Rout 2>&1 &
