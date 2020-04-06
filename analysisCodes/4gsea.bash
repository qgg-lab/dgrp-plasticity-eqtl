# 1. run GSEA for canalization/decanalization score
# ============================================================

~/software/R-3.2.2/bin/Rscript -e 'load("qg/female.adj.qgVarHet.RData"); load("qg/female.adj.qgSingleTemp.RData"); sig.h2 <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05 | p.adjust(single.temp[, 12], method = "BH") < 0.05); h2.18c <- as.numeric(var.het[, 8])[sig.h2]; h2.25c <- as.numeric(var.het[, 7])[sig.h2]; score <- (h2.18c - h2.25c)/(h2.18c + h2.25c); write.table(na.omit(cbind(gene.name[sig.h2], score)), file = "gsea/female.decan.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/software/R-3.2.2/bin/Rscript -e 'load("qg/male.adj.qgVarHet.RData"); load("qg/male.adj.qgSingleTemp.RData"); sig.h2 <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05 | p.adjust(single.temp[, 12], method = "BH") < 0.05); h2.18c <- as.numeric(var.het[, 8])[sig.h2]; h2.25c <- as.numeric(var.het[, 7])[sig.h2]; score <- (h2.18c - h2.25c)/(h2.18c + h2.25c); write.table(na.omit(cbind(gene.name[sig.h2], score)), file = "gsea/male.decan.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/software/R-3.2.2/bin/Rscript prepGSEA.R gsea/female.decan.score ~/flybase/fb-r5.57/gene_go.table 20 gsea/female.decan.gsea.data.RData > log/female.decan.prepGSEA.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript prepGSEA.R gsea/male.decan.score ~/flybase/fb-r5.57/gene_go.table 20 gsea/male.decan.gsea.data.RData > log/male.decan.prepGSEA.Rout 2>&1 &

~/software/R-3.2.2/bin/Rscript gsea.R gsea/female.decan.gsea.data.RData gsea/female.decan.score 1000 gsea/female.decan.gsea.RData > log/female.decan.gsea.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript gsea.R gsea/male.decan.gsea.data.RData gsea/male.decan.score 1000 gsea/male.decan.gsea.RData > log/male.decan.gsea.Rout 2>&1 &

~/software/R-3.2.2/bin/Rscript gseaFDR.R gsea/female.decan.gsea.RData gsea/female.decan.gsea.fdr.RData > log/female.decan.gsea.fdr.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript gseaFDR.R gsea/male.decan.gsea.RData gsea/male.decan.gsea.fdr.RData > log/male.decan.gsea.fdr.Rout 2>&1 &

# extract examples
# ============================================================

~/software/R-3.2.2/bin/Rscript gseaExample.R gsea/female.decan.gsea.data.RData gsea/female.decan.score "MF>GO:0005214" gsea/female.decan.gsea.example.RData > log/female.decan.gsea.example.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript gseaExample.R gsea/gsea.RData gsea/male.decan.score "MF>GO:0005549,MF>GO:0003700" gsea/male.decan.gsea.example.RData > log/male.decan.gsea.example.Rout 2>&1 &

# 3. prepare data for GEI GSEA
# ============================================================

~/software/R-3.2.2/bin/Rscript -e 'load("qg/female.adj.qgGxE.RData"); gei <- as.numeric(gxe[, 6])/(as.numeric(gxe[, 6]) + as.numeric(gxe[, 7])); int.fdr <- p.adjust(as.numeric(gxe[, 10]), method = "BH"); line.fdr <- p.adjust(as.numeric(gxe[, 9]), method = "BH"); write.table(na.omit(cbind(gene.name[int.fdr < 0.05 | line.fdr < 0.05], 2*gei[int.fdr < 0.05 | line.fdr < 0.05] - 1)), file = "gsea/female.gei.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/software/R-3.2.2/bin/Rscript -e 'load("qg/male.adj.qgGxE.RData"); gei <- as.numeric(gxe[, 6])/(as.numeric(gxe[, 6]) + as.numeric(gxe[, 7])); int.fdr <- p.adjust(as.numeric(gxe[, 10]), method = "BH"); line.fdr <- p.adjust(as.numeric(gxe[, 9]), method = "BH"); write.table(na.omit(cbind(gene.name[int.fdr < 0.05 | line.fdr < 0.05], 2*gei[int.fdr < 0.05 | line.fdr < 0.05] - 1)), file = "gsea/male.gei.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/software/R-3.2.2/bin/Rscript prepGSEA.R gsea/female.gei.score ~/flybase/fb-r5.57/gene_go.table 20 gsea/female.gei.gsea.data.RData > log/female.gei.prepGSEA.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript prepGSEA.R gsea/male.gei.score ~/flybase/fb-r5.57/gene_go.table 20 gsea/male.gei.gsea.data.RData > log/male.gei.prepGSEA.Rout 2>&1 &

~/software/R-3.2.2/bin/Rscript gsea.R gsea/female.gei.gsea.data.RData gsea/male.gei.score 1000 gsea/male.gei.gsea.RData > log/male.gei.gsea.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript gsea.R gsea/male.gei.gsea.data.RData gsea/female.gei.score 1000 gsea/female.gei.gsea.RData > log/female.gei.gsea.Rout 2>&1 &

~/software/R-3.2.2/bin/Rscript gseaFDR.R gsea/female.gei.gsea.RData gsea/female.gei.gsea.fdr.RData > log/female.gei.gsea.fdr.Rout 2>&1 &
~/software/R-3.2.2/bin/Rscript gseaFDR.R gsea/male.gei.gsea.RData gsea/male.gei.gsea.fdr.RData > log/male.gei.gsea.fdr.Rout 2>&1 &

