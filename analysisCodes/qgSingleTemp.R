# ============================================
# = quantitative genetics of gene expression =
# ============================================

# use normal quantiles to transform gene expression
# within each temperature
# ============================================================

args <- commandArgs(TRUE) # args <- c("tilingInt/female.gene.exp.adjust.RData", 8)
rdata.file <- args[1]
n.cpu <- as.numeric(args[2])

# load libraries
# ============================================================

library("lmerTest")
library("car")
library("doMC")
registerDoMC(n.cpu)

# load data
# ============================================================

load(rdata.file)

# re-order wolbachia, inversions, etc.
# ============================================================

line.25c <- line.id[temp == "25C"]
line.18c <- line.id[temp == "18C"]
wolba.25c <- wolba[temp == "25C"]
wolba.18c <- wolba[temp == "18C"]

# perform QG analysis
# ============================================================

gene.name <- rownames(gene.exp.adj)
# group.id <- factor(c(rep(1, length(levels(line.18c))), rep(2, length(levels(line.25c)))))
# residual.group.id <- factor(c(rep(1, length(line.18c)), rep(2, length(line.25c))))

single.temp <- foreach (i = 1:nrow(gene.exp.adj), .combine = rbind) %dopar% {
  
  cat(gene.name[i], "\n")
  res <- rep(NA, 12)
  
  y <- unlist(gene.exp.adj[i, ])
  
  # fit model in temperature separately
  # ============================================================
  
  y.25c <- y[temp == "25C"]
  y.18c <- y[temp == "18C"]
  
  lmer.fit.18c <- lmer(y.18c ~ (1|line.18c)) # rand(lmer.fit.18c)$rand.table[1, 3] is p value for random effect
  # as.data.frame(VarCorr(lmer.fit.18c))[, 4] is variance components
  lmer.fit.25c <- lmer(y.25c ~ (1|line.25c)) 
  res[1:3] <- c(as.data.frame(VarCorr(lmer.fit.18c))[, 4], rand(lmer.fit.18c)$rand.table[1, 3])
  res[4:6] <- c(as.data.frame(VarCorr(lmer.fit.25c))[, 4], rand(lmer.fit.25c)$rand.table[1, 3])
  
  # fit model in the presence of wolbachia
  # ============================================================
  
  lmer.fit.18c.wolba <- lmer(y.18c ~ wolba.18c + (1|line.18c))
  lmer.fit.25c.wolba <- lmer(y.25c ~ wolba.25c + (1|line.25c))
  res[7:9] <- c(as.data.frame(VarCorr(lmer.fit.18c.wolba))[, 4], rand(lmer.fit.18c.wolba)$rand.table[1, 3])
  res[10:12] <- c(as.data.frame(VarCorr(lmer.fit.25c.wolba))[, 4], rand(lmer.fit.25c.wolba)$rand.table[1, 3])

  # return result
  # ============================================================
  
  res
  
}

colnames(single.temp) <- c("var.line.18c", "var.e.18c", "p.line.18c", "var.line.25c", "var.e.25c", "p.line.25c", "wolba.var.line.18c", "wolba.var.e.18c", "wolba.p.line.18c", "wolba.var.line.25c", "wolba.var.e.25c", "wolba.p.line.25c")
save(single.temp, gene.name, file = args[3])

# session info
# ============================================================

sessionInfo()
