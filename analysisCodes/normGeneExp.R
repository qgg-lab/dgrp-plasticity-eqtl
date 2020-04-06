# =======================================================
# = R script to normalize and summarize gene expression =
# =======================================================

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

bad.array <- scan("tilingInt/bad.array.list", what = "", quiet = TRUE)
exp.line.id <- scan("geno/exp.185line.id", what = "", quiet = TRUE)

female.18c <- female.18c[, -which(colnames(female.18c) %in% bad.array)]
female.25c <- female.25c[, -which(colnames(female.25c) %in% bad.array)]
male.18c <- male.18c[, -which(colnames(male.18c) %in% bad.array)]
male.25c <- male.25c[, -which(colnames(male.25c) %in% bad.array)]

# extract line number
female.18c <- female.18c[, which(paste("line_", gsub("-F.*$", "", colnames(female.18c)), sep = "") %in% exp.line.id)]
female.25c <- female.25c[, which(paste("line_", gsub("^T", "", gsub("-F.*$", "", colnames(female.25c))), sep = "") %in% exp.line.id)]
male.18c <- male.18c[, which(paste("line_", gsub("-M.*$", "", colnames(male.18c)), sep = "") %in% exp.line.id)]
male.25c <- male.25c[, which(paste("line_", gsub("^T", "", gsub("-M.*$", "", colnames(male.25c))), sep = "") %in% exp.line.id)]

# quantile normalization
# ============================================================

female.norm <- normalizeQuantiles(cbind(female.18c, female.25c))
male.norm <- normalizeQuantiles(cbind(male.18c, male.25c))

# median polish
# ============================================================

female.gene.exp <- t(sapply(split(female.norm, probe.gene.assoc[, 2]), function(x) { suppressWarnings(y <- medpolish(matrix(log2(x), ncol = ncol(female.norm)), trace.iter = F)); return(y$overall + y$col) } ))
colnames(female.gene.exp) <- colnames(female.norm)

save(female.gene.exp, probe.gene.assoc, file = "tilingInt/female.gene.exp.RData")

male.gene.exp <- t(sapply(split(male.norm, probe.gene.assoc[, 2]), function(x) { suppressWarnings(y <- medpolish(matrix(log2(x), ncol = ncol(male.norm)), trace.iter = F)); return(y$overall + y$col) } ))
colnames(male.gene.exp) <- colnames(male.norm)

save(male.gene.exp, probe.gene.assoc, file = "tilingInt/male.gene.exp.RData")

sessionInfo()
