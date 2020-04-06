# ========================================
# = extract example GSEA result for plot =
# ========================================

args <- commandArgs(TRUE)

load(args[1]) # load GSEA annotation data
score.list <- read.table(args[2], header = FALSE, as.is = TRUE)
go.list <- matrix(unlist(strsplit(unlist(strsplit(args[3], split = ",")), split = ">")), ncol = 2, byrow = TRUE)

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

gsea <- function(gene.set, score, go) {
  
  N <- nrow(score)
  
  p <- rep(-1, N)
  S <- sort(match(gene.set, score[, 1]))
  p[S] <- abs(score[S, 2])/sum(abs(score[S, 2]))
  p[-S] <- p[-S]/(N - length(S))
  
  return(list(p = p, S = S, go = go, gene.set = gene.set, score = score))
  
}

# loop the main program
# ============================================================

bp.score <- score.list[score.list[, 1] %in% bp.gene.id, ]
mf.score <- score.list[score.list[, 1] %in% mf.gene.id, ]
cc.score <- score.list[score.list[, 1] %in% cc.gene.id, ]

gsea.example <- list()


for (i in 1:nrow(go.list)) {
  
  if (go.list[i, 1] == "MF") {
    
    gsea.example[[i]] <- gsea(mf.gene.set[[which(names(mf.gene.set) == go.list[i, 2])]], mf.score, go.list[i, 2])
    
  }
  
  if (go.list[i, 1] == "BP") {
    
    gsea.example[[i]] <- gsea(bp.gene.set[[which(names(bp.gene.set) == go.list[i, 2])]], bp.score, go.list[i, 2])
    
  }
  
  if (go.list[i, 1] == "CC") {
    
    gsea.example[[i]] <- gsea(cc.gene.set[[which(names(cc.gene.set) == go.list[i, 2])]], cc.score, go.list[i, 2])
    
  }
  
}

# save results
# ============================================================

save(gsea.example, file = args[4])
