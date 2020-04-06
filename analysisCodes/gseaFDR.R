# ==========================
# = calculate FDR for GSEA =
# ==========================

args <- commandArgs(TRUE) # args <- c("gsea/female.gei.gsea.RData")
load(args[1])

# 1. normalize es score using means of the permuted data sets
# ============================================================

bp.mean.es.pos <- rowMeans(bp.es.perm[, seq(1, ncol(bp.es.perm), 2)])
bp.mean.es.neg <- rowMeans(bp.es.perm[, seq(2, ncol(bp.es.perm), 2)])

bp.es.pos <- bp.es[, 1]/bp.mean.es.pos
bp.es.neg <- bp.es[, 2]/bp.mean.es.neg

bp.mean.es.pos.perm <- apply(bp.es.perm[, seq(1, ncol(bp.es.perm), 2)], 2, function(x) { return(x/bp.mean.es.pos) })
bp.mean.es.neg.perm <- apply(bp.es.perm[, seq(2, ncol(bp.es.perm), 2)], 2, function(x) { return(x/bp.mean.es.neg) })

mf.mean.es.pos <- rowMeans(mf.es.perm[, seq(1, ncol(mf.es.perm), 2)])
mf.mean.es.neg <- rowMeans(mf.es.perm[, seq(2, ncol(mf.es.perm), 2)])

mf.es.pos <- mf.es[, 1]/mf.mean.es.pos
mf.es.neg <- mf.es[, 2]/mf.mean.es.neg

mf.mean.es.pos.perm <- apply(mf.es.perm[, seq(1, ncol(mf.es.perm), 2)], 2, function(x) { return(x/mf.mean.es.pos) })
mf.mean.es.neg.perm <- apply(mf.es.perm[, seq(2, ncol(mf.es.perm), 2)], 2, function(x) { return(x/mf.mean.es.neg) })

cc.mean.es.pos <- rowMeans(cc.es.perm[, seq(1, ncol(cc.es.perm), 2)])
cc.mean.es.neg <- rowMeans(cc.es.perm[, seq(2, ncol(cc.es.perm), 2)])

cc.es.pos <- cc.es[, 1]/cc.mean.es.pos
cc.es.neg <- cc.es[, 2]/cc.mean.es.neg

cc.mean.es.pos.perm <- apply(cc.es.perm[, seq(1, ncol(cc.es.perm), 2)], 2, function(x) { return(x/cc.mean.es.pos) })
cc.mean.es.neg.perm <- apply(cc.es.perm[, seq(2, ncol(cc.es.perm), 2)], 2, function(x) { return(x/cc.mean.es.neg) })

# 2. calculate FDR
# ============================================================

bp.es.pos.sorted <- sort(bp.es.pos, decreasing = TRUE)
bp.es.neg.sorted <- sort(bp.es.neg, decreasing = TRUE)

bp.pos.fdr <- data.frame(go = names(bp.es.pos.sorted), nes = bp.es.pos.sorted, fdr = NA)
bp.neg.fdr <- data.frame(go = names(bp.es.neg.sorted), nes = bp.es.neg.sorted, fdr = NA)

for (i in 1:length(bp.es.pos.sorted)) {
  
  bp.pos.fdr[i, 3] <- sum(bp.mean.es.pos.perm >= bp.es.pos.sorted[i])/(ncol(bp.es.perm)/2)/i
  
}

for (i in 1:length(bp.es.neg.sorted)) {
  
  bp.neg.fdr[i, 3] <- sum(bp.mean.es.neg.perm >= bp.es.neg.sorted[i])/(ncol(bp.es.perm)/2)/i
  
}

bp.pos.fdr[, 3] <- ifelse(cummax(bp.pos.fdr[, 3]) > 1, 1, cummax(bp.pos.fdr[, 3]))
bp.neg.fdr[, 3] <- ifelse(cummax(bp.neg.fdr[, 3]) > 1, 1, cummax(bp.neg.fdr[, 3]))


mf.es.pos.sorted <- sort(mf.es.pos, decreasing = TRUE)
mf.es.neg.sorted <- sort(mf.es.neg, decreasing = TRUE)

mf.pos.fdr <- data.frame(go = names(mf.es.pos.sorted), nes = mf.es.pos.sorted, fdr = NA)
mf.neg.fdr <- data.frame(go = names(mf.es.neg.sorted), nes = mf.es.neg.sorted, fdr = NA)

for (i in 1:length(mf.es.pos.sorted)) {
  
  mf.pos.fdr[i, 3] <- sum(mf.mean.es.pos.perm >= mf.es.pos.sorted[i])/(ncol(mf.es.perm)/2)/i
  
}

for (i in 1:length(mf.es.neg.sorted)) {
  
  mf.neg.fdr[i, 3] <- sum(mf.mean.es.neg.perm >= mf.es.neg.sorted[i])/(ncol(mf.es.perm)/2)/i
  
}

mf.pos.fdr[, 3] <- ifelse(cummax(mf.pos.fdr[, 3]) > 1, 1, cummax(mf.pos.fdr[, 3]))
mf.neg.fdr[, 3] <- ifelse(cummax(mf.neg.fdr[, 3]) > 1, 1, cummax(mf.neg.fdr[, 3]))


cc.es.pos.sorted <- sort(cc.es.pos, decreasing = TRUE)
cc.es.neg.sorted <- sort(cc.es.neg, decreasing = TRUE)

cc.pos.fdr <- data.frame(go = names(cc.es.pos.sorted), nes = cc.es.pos.sorted, fdr = NA)
cc.neg.fdr <- data.frame(go = names(cc.es.neg.sorted), nes = cc.es.neg.sorted, fdr = NA)

for (i in 1:length(cc.es.pos.sorted)) {
  
  cc.pos.fdr[i, 3] <- sum(cc.mean.es.pos.perm >= cc.es.pos.sorted[i])/(ncol(cc.es.perm)/2)/i
  
}

for (i in 1:length(cc.es.neg.sorted)) {
  
  cc.neg.fdr[i, 3] <- sum(cc.mean.es.neg.perm >= cc.es.neg.sorted[i])/(ncol(cc.es.perm)/2)/i
  
}

cc.pos.fdr[, 3] <- ifelse(cummax(cc.pos.fdr[, 3]) > 1, 1, cummax(cc.pos.fdr[, 3]))
cc.neg.fdr[, 3] <- ifelse(cummax(cc.neg.fdr[, 3]) > 1, 1, cummax(cc.neg.fdr[, 3]))

save(bp.pos.fdr, bp.neg.fdr, mf.pos.fdr, mf.neg.fdr, cc.pos.fdr, cc.neg.fdr, file = args[2])
