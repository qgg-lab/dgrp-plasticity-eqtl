# ================================================
# = prepare annotation data sets for GSEA to use =
# ================================================

args <- commandArgs(TRUE)

gene.id <- read.table(args[1], header = FALSE, as.is = TRUE)[, 1]
annot <- read.table(args[2], header = TRUE, as.is = TRUE, sep = " ") # annotation table
min.go <- as.numeric(args[3])

# subset annotations
# ============================================================

mf.annot <- annot[annot[, 1] %in% gene.id & annot[, 2] != "", c(1, 2)]
bp.annot <- annot[annot[, 1] %in% gene.id & annot[, 3] != "", c(1, 3)]
cc.annot <- annot[annot[, 1] %in% gene.id & annot[, 4] != "", c(1, 4)]

mf.gene.id <- mf.annot[, 1]
bp.gene.id <- bp.annot[, 1]
cc.gene.id <- cc.annot[, 1]

unique.mf <- unique(unlist(strsplit(mf.annot[, 2], split = ",")))
unique.bp <- unique(unlist(strsplit(bp.annot[, 2], split = ",")))
unique.cc <- unique(unlist(strsplit(cc.annot[, 2], split = ",")))

# prepare list of GO genes
# ============================================================

mf.gene.set <- list()
i = 1
for (go in unique.mf) {
  go.genes <- mf.annot[grepl(go, mf.annot[, 2]), 1]
  if (length(go.genes) >= min.go) {
    mf.gene.set[[i]] <- go.genes
    names(mf.gene.set)[i] <- go
    i = i + 1
  }
}


bp.gene.set <- list()
i = 1
for (go in unique.bp) {
  go.genes <- bp.annot[grepl(go, bp.annot[, 2]), 1]
  if (length(go.genes) >= min.go) {
    bp.gene.set[[i]] <- go.genes
    names(bp.gene.set)[i] <- go
    i = i + 1
  }
}


cc.gene.set <- list()
i = 1
for (go in unique.cc) {
  go.genes <- cc.annot[grepl(go, cc.annot[, 2]), 1]
  if (length(go.genes) >= min.go) {
    cc.gene.set[[i]] <- go.genes
    names(cc.gene.set)[i] <- go
    i = i + 1
  }
}

# save data set
save(mf.gene.set, mf.gene.id, bp.gene.set, bp.gene.id, cc.gene.set, cc.gene.id, file = args[4])
