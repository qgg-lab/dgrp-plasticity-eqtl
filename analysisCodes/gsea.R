# =======================
# = GSEA main R program =
# =======================

args <- commandArgs(TRUE) # args <- c("gsea/female.gei.gsea.data.RData", "gsea/female.gei.score", 1000, "gsea/female.gei.gsea.RData")

load(args[1]) # load GSEA annotation data
score.list <- read.table(args[2], header = FALSE, as.is = TRUE)

# rank the score list
# ============================================================

score.list <- score.list[order(score.list[, 2], decreasing = TRUE), ]
rownames(score.list) <- score.list[, 1]

# function to calculate ES
# takes 1) gene.set (set gene IDs), which has been processed
# to contain only gene IDs present in the score list
# this is not exactly necessary but will improve efficiency
# if pre-computed
# 2) score list
# ============================================================

gsea <- function(gene.set, score) {
  
  N <- nrow(score)
  
  p <- rep(-1, N)
  S <- sort(match(gene.set, score[, 1]))
  p[S] <- abs(score[S, 2])/sum(abs(score[S, 2]))
  p[-S] <- p[-S]/(N - length(S))
  
  p.cumsum <- cumsum(p)
  return(c(
    max(c(max(p.cumsum), 0)),
    min(c(0, min(p.cumsum)))
  ))
  
}

# loop the main program
# ============================================================

bp.es <- matrix(NA, nrow = length(bp.gene.set), ncol = 2)
mf.es <- matrix(NA, nrow = length(mf.gene.set), ncol = 2)
cc.es <- matrix(NA, nrow = length(cc.gene.set), ncol = 2)

bp.score <- score.list[score.list[, 1] %in% bp.gene.id, ]
mf.score <- score.list[score.list[, 1] %in% mf.gene.id, ]
cc.score <- score.list[score.list[, 1] %in% cc.gene.id, ]

rownames(bp.es) <- names(bp.gene.set)
rownames(mf.es) <- names(mf.gene.set)
rownames(cc.es) <- names(cc.gene.set)

for (i in 1:length(bp.gene.set)) {
  
  bp.es[i, ] <- gsea(bp.gene.set[[i]], bp.score)
  
}


for (i in 1:length(mf.gene.set)) {
  
  mf.es[i, ] <- gsea(mf.gene.set[[i]], mf.score)
  
}


for (i in 1:length(cc.gene.set)) {
  
  cc.es[i, ] <- gsea(cc.gene.set[[i]], cc.score)
  
}

# permutation results
# ============================================================

bp.es.perm <- matrix(NA, nrow = length(bp.gene.set), ncol = 2*as.numeric(args[3]))
mf.es.perm <- matrix(NA, nrow = length(mf.gene.set), ncol = 2*as.numeric(args[3]))
cc.es.perm <- matrix(NA, nrow = length(cc.gene.set), ncol = 2*as.numeric(args[3]))

set.seed(as.numeric(args[3]))

for (i in 1:as.numeric(args[3])) {
  
  score.permute <- data.frame(gene = sample(score.list[, 1]), score = score.list[, 2])
  
  bp.score.permute <- score.permute[score.permute[, 1] %in% bp.gene.id, ]
  mf.score.permute <- score.permute[score.permute[, 1] %in% mf.gene.id, ]
  cc.score.permute <- score.permute[score.permute[, 1] %in% cc.gene.id, ]
  
  for (j in 1:length(bp.gene.set)) {
    bp.es.perm[j, c(2*i - 1, 2*i)] <- gsea(bp.gene.set[[j]], bp.score.permute)
  }
  
  for (j in 1:length(mf.gene.set)) {
    mf.es.perm[j, c(2*i - 1, 2*i)] <- gsea(mf.gene.set[[j]], mf.score.permute)
  }
  
  for (j in 1:length(cc.gene.set)) {
    cc.es.perm[j, c(2*i - 1, 2*i)] <- gsea(cc.gene.set[[j]], cc.score.permute)
  }
  
  cat(i, "\n")
  
}

# save results
# ============================================================

save(bp.es, mf.es, cc.es, bp.es.perm, mf.es.perm, cc.es.perm, file = args[4])
