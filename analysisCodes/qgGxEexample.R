# ================================================================
# = blup for Hsp proteins only, without adjusting for covariates =
# ================================================================

# use normal quantiles to transform gene expression
# within each temperature
# ============================================================

args <- commandArgs(TRUE) # args <- c("qg/female.adj.qgBLUP.RData", "qg/female.adj.qgGxE.RData", "cuff/gene.info")

# load data
# ============================================================

load(args[1])
load(args[2])
gene.info <- read.table(args[3], header = FALSE, as.is = TRUE, quote = "\"")

# find hsp genes
# ============================================================

hs.gene <- match(gene.info[grepl("Hs", gene.info[, 2]), 1], gene.name)

# calculate rank correlation
# ============================================================

gene.rank.cor <- numeric(length(gene.name))

for (i in 1:length(gene.rank.cor)) {
  
  gene.rank.cor[i] <- cor(blup.18c[i, ], blup.25c[i, ], method = "spearman")
  
}

# calculate sample variance ratio
# ============================================================

var.ratio <- numeric(length(gene.name))

for (i in 1:length(var.ratio)) {
  
  var.ratio[i] <- var(blup.18c[i, ])/var(blup.25c[i, ])
  
}

# calculate gei
# ============================================================

gei <- gxe[, 6]/(gxe[, 7] + gxe[, 6])

# range
# ============================================================

range.18c <- apply(blup.18c, 1, max) - apply(blup.18c, 1, min)
range.25c <- apply(blup.25c, 1, max) - apply(blup.25c, 1, min)


# find candidate genes
# ============================================================

candidate.gene <- unique(c(which( (gene.rank.cor > 0.75 & var.ratio < 1.5 & var.ratio > 0.5 & gei < 0.1 & range.18c > 0.5 & range.25c > 0.5) |
                                  (gene.rank.cor < 0.2 & (var.ratio > 2 | var.ratio < 0.5) & gei > 0.8 & range.18c > 0.5 & range.25c > 0.5) ),
                           hs.gene))
gene.info <- gene.info[match(gene.name[candidate.gene], gene.info[, 1]), ]
blup.25c <- blup.25c[candidate.gene, ]
blup.18c <- blup.18c[candidate.gene, ]
gxe <- gxe[candidate.gene, ]

save(gene.info, blup.25c, blup.18c, gxe, file = args[4])

# session info
# ============================================================

sessionInfo()
