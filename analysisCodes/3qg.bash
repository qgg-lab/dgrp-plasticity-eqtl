# 1. QG and blup
# ============================================================

srun -w node3 ~/software/R-3.2.2/bin/Rscript qgSingleTemp.R tilingInt/female.gene.exp.adjust.RData 16 qg/female.adj.qgSingleTemp.RData > log/female.adj.qgSingleTemp.Rout 2>&1 &
srun -w node4 ~/software/R-3.2.2/bin/Rscript qgSingleTemp.R tilingInt/male.gene.exp.adjust.RData 16 qg/male.adj.qgSingleTemp.RData > log/male.adj.qgSingleTemp.Rout 2>&1 &

srun -w node3 ~/software/R-3.2.2/bin/Rscript qgVarHet.R tilingInt/female.gene.exp.adjust.RData 16 qg/female.adj.qgVarHet.RData > log/female.adj.qgVarHet.Rout 2>&1 &
srun -w node4 ~/software/R-3.2.2/bin/Rscript qgVarHet.R tilingInt/male.gene.exp.adjust.RData 16 qg/male.adj.qgVarHet.RData > log/male.adj.qgVarHet.Rout 2>&1 &

srun -w node3 ~/software/R-3.2.2/bin/Rscript qgGxE.R tilingInt/female.gene.exp.adjust.RData 16 qg/female.adj.qgGxE.RData > log/female.adj.qgGxE.Rout 2>&1 &
srun -w node4 ~/software/R-3.2.2/bin/Rscript qgGxE.R tilingInt/male.gene.exp.adjust.RData 16 qg/male.adj.qgGxE.RData > log/male.adj.qgGxE.Rout 2>&1 &

srun -w node3 ~/software/R-3.2.2/bin/Rscript qgBLUP.R tilingInt/female.gene.exp.adjust.RData ~/dgrp/adjustData.RData 16 qg/female.adj.qgBLUP.RData > log/female.adj.qgBLUP.Rout 2>&1 &
srun -w node4 ~/software/R-3.2.2/bin/Rscript qgBLUP.R tilingInt/male.gene.exp.adjust.RData ~/dgrp/adjustData.RData 16 qg/male.adj.qgBLUP.RData > log/male.adj.qgBLUP.Rout 2>&1 &

# examples of GxE small and large (find some Hsps)
# ============================================================

srun ~/software/R-3.2.2/bin/Rscript qgGxEexample.R qg/female.adj.qgBLUP.RData qg/female.adj.qgGxE.RData cuff/gene.info qg/female.adj.qgGxEexample.RData > log/female.adj.qgGxEexample.Rout 2>&1 &
