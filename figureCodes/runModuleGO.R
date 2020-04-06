# ========================================
# = test for enrichemnt of GO in modules =
# ========================================

args <- commandArgs(TRUE) # args <- c("reportData/female.module1.RData", "reportData/gene_go.table", 20, "reportData/female.module1.go.RData")
load(args[1])
annot.file <- args[2]
min.go <- as.numeric(args[3])
out.file <- args[4]
library("GO.db")

annot <- read.table(annot.file, header = TRUE, as.is = TRUE, sep = " ")

mf.annot <- annot[annot[, 1] %in% all.gene & annot[, 2] != "", c(1, 2)]
bp.annot <- annot[annot[, 1] %in% all.gene & annot[, 3] != "", c(1, 3)]
cc.annot <- annot[annot[, 1] %in% all.gene & annot[, 4] != "", c(1, 4)]



# function to test for significance
go.stat <- function(annot.table, list.gene, type, min.go) {
  
  # identify go terms to be tested
  go.terms <- table(unlist(strsplit(annot.table[, 2], split = ",")))
  go.terms <- go.terms[go.terms >= min.go]
  n.total <- nrow(annot.table)
  # filter gene
  list.gene <- intersect(list.gene, annot.table[, 1])
  
  res <- NULL
  k <- length(list.gene)
  
  for (go in names(go.terms)) {
    m <- sum(grepl(go, annot.table[, 2]))
    n <- n.total - m
    q = length(intersect(list.gene, annot.table[grepl(go, annot.table[, 2]), 1]))
    res <- rbind(res, c(type, go, Term(go), m, n.total, k, q, phyper(q, m, n, k, lower.tail = F)))
  }
  res <- cbind(res, p.adjust(as.numeric(res[, 8]), method = "BH"))
  return (res)
}

mf.out <- go.stat(mf.annot, module.gene, "MF", min.go)
bp.out <- go.stat(bp.annot, module.gene, "BP", min.go)
cc.out <- go.stat(cc.annot, module.gene, "CC", min.go)

save(mf.out, bp.out, cc.out, file = out.file)
