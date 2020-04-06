# ==================================
# = co-expression network analysis =
# ==================================

library("WGCNA")
library("NetRep")
args <- commandArgs(TRUE)
load(args[1])

# 1. run wgcna
# ============================================================

allowWGCNAThreads(nThreads = 8)
blup.18c <- t(blup.18c.data[match(intersect(sig.18c.gene, sig.25c.gene), sig.gene.name), ])
colnames(blup.18c) <- intersect(sig.18c.gene, sig.25c.gene)
blup.25c <- t(blup.25c.data[match(intersect(sig.18c.gene, sig.25c.gene), sig.gene.name), ])
colnames(blup.25c) <- intersect(sig.18c.gene, sig.25c.gene)

rownames(blup.18c) <- 1:nrow(blup.18c)
rownames(blup.25c) <- 1:nrow(blup.25c)

power.scan = 1:20

blup.25c.soft.threshold = pickSoftThreshold(blup.25c, dataIsExpr = TRUE, powerVector = power.scan, corFnc = cor, corOptions = list(use = 'pairwise', method = "pearson"), networkType = "unsigned")

blup.25c.adj <- adjacency(blup.25c, power = blup.25c.soft.threshold$powerEstimate, corOptions = list(use = 'pairwise', method = "pearson"))
blup.25c.cor <- cor(blup.25c, use = "pairwise", method = "pearson")
blup.25c.tom <- TOMdist(blup.25c.adj)
blup.25c.clust <- hclust(as.dist(blup.25c.tom), method = "average")

blup.25c.tree <- cutreeDynamic(dendro = blup.25c.clust, method = "tree", minClusterSize = 20)
blup.25c.col <- labels2colors(blup.25c.tree)

colnames(blup.25c.tom) <- colnames(blup.25c)
rownames(blup.25c.tom) <- colnames(blup.25c)
names(blup.25c.tree) <- colnames(blup.25c)

# blup.18c.soft.threshold = pickSoftThreshold(blup.18c, dataIsExpr = TRUE, powerVector = power.scan, corFnc = cor, corOptions = list(use = 'pairwise', method = "pearson"), networkType = "unsigned")
# note to use the 25c threshold
blup.18c.adj <- adjacency(blup.18c, power = blup.25c.soft.threshold$powerEstimate, corOptions = list(use = 'pairwise', method = "pearson"))
blup.18c.cor <- cor(blup.18c, use = "pairwise", method = "pearson")
blup.18c.tom <- TOMdist(blup.18c.adj)
blup.18c.clust <- hclust(as.dist(blup.18c.tom), method = "average")

blup.18c.tree <- cutreeDynamic(dendro = blup.18c.clust, method = "tree", minClusterSize = 20)
blup.18c.col <- labels2colors(blup.18c.tree)

colnames(blup.18c.tom) <- colnames(blup.18c)
rownames(blup.18c.tom) <- colnames(blup.18c)

set.seed(1)
data_list <- list(cohort1 = blup.25c, cohort2 = blup.18c)
correlation_list <- list(cohort1 = cor(blup.25c), cohort2 = cor(blup.18c))
network_list <- list(cohort1 = blup.25c.tom, cohort2 = blup.18c.tom)
mod.pres <- modulePreservation(network = network_list, data = data_list, correlation = correlation_list, moduleAssignments = blup.25c.tree, discovery = "cohort1", test = "cohort2", nPerm = 10000)

save(mod.pres, blup.25c.soft.threshold, blup.25c.cor, blup.25c.tom, blup.25c.clust, blup.25c.tree, blup.18c.cor, blup.18c.tom, blup.18c.clust, blup.18c.tree, file = args[2])
