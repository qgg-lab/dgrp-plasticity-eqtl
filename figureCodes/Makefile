# ===============================================
# = make file for generating tables and figures =
# ===============================================

all: report/Table_MappingSummary.xlsx report/Figure_RNASeqExpressionScatterPlot.pdf report/Figure_PreliminaryNormalizationOutlier.pdf report/Figure_PCA.pdf report/Figure_Decanalization.pdf report/Figure_NoWolbaDecanalization.pdf report/Table_VarComp.xlsx report/Table_DecanGSEA.xlsx report/Figure_GxEHist.pdf report/Figure_DecanGSEAExample.pdf reportData/female.wgcna.RData reportData/male.wgcna.RData report/Figure_Network.pdf reportData/female.module1.go.RData report/Figure_CorrDiff.pdf report/Table_GEIGSEA.xlsx report/Table_eqtl.xlsx report/Figure_eQTLsummary.pdf report/Table_female.module1.go.xlsx report/Table_NetRep.xlsx report/Table_GeneComparison.xlsx report/Figure_Connectivity.pdf
	
FONT="Myriad Pro"

# mapping statistics table in excel format
# ============================================================

report/Table_MappingSummary.xlsx:	reportCode/tableMapSummary.R reportData/female18c.summary reportData/male18c.summary reportData/female25c.summary reportData/male25c.summary
	Rscript reportCode/tableMapsummary.R reportData/female18c.summary reportData/male18c.summary reportData/female25c.summary reportData/male25c.summary report/Table_MappingSummary.xlsx
	cp report/Table_MappingSummary.xlsx paper/TableS1.xlsx

# figure for temperature effect based on RNA-Seq
# ============================================================

report/Figure_RNASeqExpressionScatterPlot.pdf: reportCode/figureRNASeqExpressionScatterPlot.R reportData/genes.fpkm_table reportData/gene.symbol
	Rscript reportCode/figureRNASeqExpressionScatterPlot.R reportData/genes.fpkm_table reportData/gene.symbol report/Figure_RNASeqExpressionScatterPlot.pdf $(FONT)
	cp report/Figure_RNASeqExpressionScatterPlot.pdf paper/figures/FigureS1.pdf

# figure for outlier detection for arrays
# ============================================================

report/Figure_PreliminaryNormalizationOutlier.pdf: reportCode/figurePrelimNormOutlier.R reportData/prelim.norm.summary.RData reportData/cel.headers.RData
	Rscript reportCode/figurePrelimNormOutlier.R reportData/prelim.norm.summary.RData reportData/cel.headers.RData report/Figure_PreliminaryNormalizationOutlier.pdf $(FONT)
	cp report/Figure_PreliminaryNormalizationOutlier.pdf paper/figures/FigureS2.pdf

# figure for batch correction
# ============================================================

report/Figure_PCA.pdf: reportCode/figurePCA.R reportData/female.pca.RData reportData/female.sva.adjust.RData reportData/male.pca.RData reportData/male.sva.adjust.RData
	Rscript reportCode/figurePCA.R reportData/female.pca.RData reportData/female.sva.adjust.RData reportData/male.pca.RData reportData/male.sva.adjust.RData report/Figure_PCA.pdf $(FONT)
	cp report/Figure_PCA.pdf paper/figures/FigureS3.pdf

# figure for variance components (adjusting for everything)
# ============================================================

report/Figure_Decanalization.pdf report/Figure_NoWolbaDecanalization.pdf: reportCode/figureDecanalization.R reportData/female.adj.qgSingleTemp.RData reportData/male.adj.qgSingleTemp.RData reportData/female.adj.qgVarHet.RData reportData/male.adj.qgVarHet.RData reportData/gene.info
	Rscript reportCode/figureDecanalizationR1.R reportData/female.adj.qgSingleTemp.RData reportData/male.adj.qgSingleTemp.RData reportData/female.adj.qgVarHet.RData reportData/male.adj.qgVarHet.RData reportData/gene.info report/Figure_Decanalization.pdf report/Figure_NoWolbaDecanalization.pdf $(FONT)
	cp report/Figure_Decanalization.pdf paper/figures/Figure2.pdf
	cp report/Figure_NoWolbaDecanalization.pdf paper/figures/FigureS4.pdf

# table for variance components
# ============================================================

report/Table_VarComp.xlsx: reportCode/tableVarComponents.R reportData/female.adj.qgSingleTemp.RData reportData/male.adj.qgSingleTemp.RData reportData/female.adj.qgVarHet.RData reportData/male.adj.qgVarHet.RData reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData reportData/gene.info
	Rscript reportCode/tableVarComponents.R reportData/female.adj.qgSingleTemp.RData reportData/male.adj.qgSingleTemp.RData reportData/female.adj.qgVarHet.RData reportData/male.adj.qgVarHet.RData reportData/female.adj.qgGxE.RData reportData/male.qgGxE.RData reportData/gene.info report/Table_VarComp.xlsx

# table for decanlization GSEA
# ============================================================

report/Table_DecanGSEA.xlsx: reportCode/tableDecanGSEA.R reportData/female.decan.gsea.fdr.RData reportData/male.decan.gsea.fdr.RData
	Rscript reportCode/tableDecanGSEA.R reportData/female.decan.gsea.fdr.RData reportData/male.decan.gsea.fdr.RData report/Table_DecanGSEA.xlsx
	
# figure for example and histograms
# ============================================================

report/Figure_GxEHist.pdf: reportCode/figureGxEsummary.R reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData reportData/female.adj.qgGxEexample.RData
	Rscript reportCode/figureGxEsummary.R reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData reportData/female.adj.qgGxEexample.RData 18 sosie 1 Cp15 report/Figure_GxEHist.pdf $(FONT)
	cp report/Figure_GxEHist.pdf paper/figures/Figure3.pdf

# figure for decanalization GSEA
# ============================================================

report/Figure_DecanGSEAExample.pdf: reportCode/figureDecanGSEAExample.R reportData/female.decan.gsea.example.RData reportData/male.decan.gsea.example.RData
	Rscript reportCode/figureDecanGSEAExample.R reportData/female.decan.gsea.example.RData reportData/male.decan.gsea.example.RData report/Figure_DecanGSEAExample.pdf $(FONT)
	cp report/Figure_DecanGSEAExample.pdf paper/figures/FigureS5.pdf

# run code for network analysis
# ============================================================

reportData/female.wgcna.RData: reportCode/runNetwork.R reportData/female.sig.gene.exp.RData
	Rscript reportCode/runNetwork.R reportData/female.sig.gene.exp.RData reportData/female.wgcna.RData
reportData/male.wgcna.RData: reportCode/runNetwork.R reportData/male.sig.gene.exp.RData
	Rscript reportCode/runNetwork.R reportData/female.module1.RData reportData/male.wgcna.RData
	
# figure for modules
# ============================================================

report/Figure_Network.pdf: reportCode/figureNetwork.R reportData/female.wgcna.RData reportData/male.wgcna.RData
	Rscript reportCode/figureNetwork.R reportData/female.wgcna.RData reportData/male.wgcna.RData reportData/female.simCorr18C.RData reportData/male.simCorr18C.RData report/Figure_Network.pdf $(FONT)
	cp report/Figure_Network.pdf paper/figures/Figure5.pdf
	
# run module enrichment
# ============================================================

reportData/female.module1.go.RData: reportCode/runModuleGo.R reportData/female.wgcna.RData reportData/gene_go.table
	Rscript -e 'load("reportData/female.wgcna.RData"); all.gene <- grep("FBgn", colnames(blup.25c.cor), value = T); module.gene <- grep("FBgn", colnames(blup.25c.cor)[which(blup.25c.tree == 1)], value = TRUE); save(all.gene, module.gene, file = "reportData/female.module1.RData");'
	Rscript reportCode/runModuleGo.R reportData/female.module1.RData reportData/gene_go.table 20 reportData/female.module1.go.RData

# produce table for female module 1 enrichment
# ============================================================

report/Table_female.module1.go.xlsx: reportCode/tableFemaleModule1Go.R reportData/female.module1.go.RData
	Rscript reportCode/tableFemaleModule1Go.R reportData/female.module1.go.RData report/Table_female.module1.go.xlsx
	
# table for GEI gsea
# ============================================================

report/Table_GEIGSEA.xlsx: reportCode/tableGEIGSEA.R reportData/female.gei.gsea.fdr.RData reportData/male.gei.gsea.fdr.RData
	Rscript reportCode/tableGEIGSEA.R reportData/female.gei.gsea.fdr.RData reportData/male.gei.gsea.fdr.RData report/Table_GEIGSEA.xlsx

# table for eqtl
# ============================================================

report/Table_eqtl.xlsx:	reportCode/tableModelSelectTable.R reportData/female.18c.fdr05.eqtl.table.txt reportData/female.25c.fdr05.eqtl.table.txt reportData/male.18c.fdr05.eqtl.table.txt reportData/male.25c.fdr05.eqtl.table.txt reportData/gene.info reportData/female.18c.genvar.id reportData/female.25c.genvar.id reportData/male.18c.genvar.id reportData/male.25c.genvar.id
	Rscript reportCode/tableModelSelectTable.R reportData/female.18c.fdr05.eqtl.table.txt reportData/female.25c.fdr05.eqtl.table.txt reportData/male.18c.fdr05.eqtl.table.txt reportData/male.25c.fdr05.eqtl.table.txt reportData/gene.info reportData/female.18c.genvar.id reportData/female.25c.genvar.id reportData/male.18c.genvar.id reportData/male.25c.genvar.id reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData report/Table_eqtl.xlsx

# figure for eqtl
# ============================================================

report/Figure_eQTLsummary.pdf: reportCode/figureSummaryeQTL.R reportData/female.18c.fdr05.eqtl.table.txt reportData/female.25c.fdr05.eqtl.table.txt reportData/male.18c.fdr05.eqtl.table.txt reportData/male.25c.fdr05.eqtl.table.txt reportData/gene.info reportData/female.18c.genvar.id reportData/female.25c.genvar.id reportData/male.18c.genvar.id reportData/male.25c.genvar.id reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData
	Rscript reportCode/figureSummaryeQTL.R reportData/female.18c.fdr05.eqtl.table.txt reportData/female.25c.fdr05.eqtl.table.txt reportData/male.18c.fdr05.eqtl.table.txt reportData/male.25c.fdr05.eqtl.table.txt reportData/gene.info reportData/female.18c.genvar.id reportData/female.25c.genvar.id reportData/male.18c.genvar.id reportData/male.25c.genvar.id reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData reportData/female.eqtl.eff.RData reportData/male.eqtl.eff.RData reportData/female.eqtl.pred.RData reportData/male.eqtl.pred.RData reportData/female.tf.summary.RData reportData/male.tf.summary.RData report/Figure_eQTLsummary.pdf $(FONT)
	cp report/Figure_eQTLsummary.pdf paper/figures/Figure4.pdf

# table for network preservation analysis
# ============================================================

report/Table_NetRep.xlsx: reportCode/tableNetRep.R reportData/female.wgcna.RData reportData/male.wgcna.RData
	Rscript reportCode/tableNetRep.R reportData/female.wgcna.RData reportData/male.wgcna.RData report/Table_NetRep.xlsx

# connectivity versus stabilizing selection
# ============================================================

report/Figure_Connectivity.pdf: reportCode/figureConnectivity.R reportData/female.wgcna.RData reportData/male.wgcna.RData reportData/sigma2.RData
	Rscript reportCode/figureConnectivity.R reportData/female.wgcna.RData reportData/male.wgcna.RData reportData/sigma2.RData report/Figure_Connectivity.pdf $(FONT)
	cp report/Figure_Connectivity.pdf paper/figures/Figure6.pdf

# gene comparison
# ============================================================

report/Table_GeneComparison.xlsx: reportCode/tableComparison.R reportData/cg.fbgn.txt reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData reportData/gene.info reportData/Hutter2008CGs.txt reportData/Levine2011CGs.txt reportData/Zhao2015dmel21cFBgns.txt reportData/Zhao2015dmel29cFBgns.txt
	Rscript reportCode/tableComparison.R reportData/cg.fbgn.txt reportData/female.adj.qgGxE.RData reportData/male.adj.qgGxE.RData reportData/gene.info reportData/Hutter2008CGs.txt reportData/Levine2011CGs.txt reportData/Zhao2015dmel21cFBgns.txt reportData/Zhao2015dmel29cFBgns.txt report/Table_GeneComparison.xlsx
	
