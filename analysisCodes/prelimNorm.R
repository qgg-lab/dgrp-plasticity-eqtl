# ============================================================
# = preliminary normalization of arrays to identify outliers =
# ============================================================

library("limma")

load("tilingInt/18c.female.bg.RData")
female.18c <- cel.int.bg
colnames(female.18c) <- cel.headers[, 1]

load("tilingInt/18c.male.bg.RData")
male.18c <- cel.int.bg
colnames(male.18c) <- cel.headers[, 1]

load("tilingInt/25c.female.bg.RData")
female.25c <- cel.int.bg
colnames(female.25c) <- cel.headers[, 1]

load("tilingInt/25c.male.bg.RData")
male.25c <- cel.int.bg
colnames(male.25c) <- cel.headers[, 1]

# normalization
# ============================================================

female.18c.norm <- normalizeQuantiles(female.18c)
male.18c.norm <- normalizeQuantiles(male.18c)
female.25c.norm <- normalizeQuantiles(female.25c)
male.25c.norm <- normalizeQuantiles(male.25c)

# median polish
# ============================================================

female.18c.gene.exp <- t(sapply(split(female.18c.norm, probe.gene.assoc[, 2]), function(x) { suppressWarnings(y <- medpolish(matrix(log2(x), ncol = ncol(female.18c.norm)), trace.iter = F)); return(y$overall + y$col) } ))
colnames(female.18c.gene.exp) <- colnames(female.18c.norm)
cat("female.18c.gene.exp\n")

male.18c.gene.exp <- t(sapply(split(male.18c.norm, probe.gene.assoc[, 2]), function(x) { suppressWarnings(y <- medpolish(matrix(log2(x), ncol = ncol(male.18c.norm)), trace.iter = F)); return(y$overall + y$col) } ))
colnames(male.18c.gene.exp) <- colnames(male.18c.norm)
cat("male.18c.gene.exp\n")

female.25c.gene.exp <- t(sapply(split(female.25c.norm, probe.gene.assoc[, 2]), function(x) { suppressWarnings(y <- medpolish(matrix(log2(x), ncol = ncol(female.25c.norm)), trace.iter = F)); return(y$overall + y$col) } ))
colnames(female.25c.gene.exp) <- colnames(female.25c.norm)
cat("male.18c.gene.exp\n")

male.25c.gene.exp <- t(sapply(split(male.25c.norm, probe.gene.assoc[, 2]), function(x) { suppressWarnings(y <- medpolish(matrix(log2(x), ncol = ncol(male.25c.norm)), trace.iter = F)); return(y$overall + y$col) } ))
colnames(male.25c.gene.exp) <- colnames(male.25c.norm)
cat("male.25c.gene.exp\n")

save(female.18c.gene.exp, male.18c.gene.exp, female.25c.gene.exp, male.25c.gene.exp, file = "tilingInt/prelim.norm.gene.exp.RData")

sessionInfo()
