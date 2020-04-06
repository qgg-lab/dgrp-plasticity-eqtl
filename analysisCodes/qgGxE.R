# ============================================
# = quantitative genetics of gene expression =
# ============================================

# use normal quantiles to transform gene expression
# within each temperature
# ============================================================

args <- commandArgs(TRUE) # args <- c("tilingInt/female.gene.exp.adjust.RData", 8)
rdata.file <- args[1]
n.cpu <- as.numeric(args[2])

library("lmerTest")
library("doMC")
registerDoMC(n.cpu)

# load data
# ============================================================

load(rdata.file)

# transformation function
# ============================================================

nqt <- function(x) {
  
  return(qnorm(rank(x)/(length(x) + 1), sd = mad(x)) + median(x))
  
}

# find line id and temperature
# ============================================================

unique.id <- factor(paste(line.id, temp, sep = "-"))

# perform QG analysis
# need
# 1. fit a linear mixed model: temp + (1|line) + (1|line:temp)
# 2. test for significance of line by temp interaction
# 3. get blup from the full model
# 4. get variance components from single temperature, not blup
# ============================================================

gene.name <- rownames(gene.exp.adj)

gxe <- foreach (i = 1:nrow(gene.exp.adj), .combine = rbind) %dopar% {
  
  cat(gene.name[i], "\n")
  
  y <- unlist(gene.exp.adj[i, ])
  
  # result vector
  res <- rep(NA, 10)
  
  # no wolbachia
  # ============================================================
  
  lmer.pool <- lmer(y ~ temp + (1|line.id) + (1|line.id:temp))
  res[1:3] <- as.data.frame(VarCorr(lmer.pool))[, 4] # var line x temp, line, res
  res[4:5] <- rand(lmer.pool)$rand.table$p.value # p value for line, line x temp
  
  # in the presence of wolbachia
  # ============================================================
  
  lmer.pool.wolba <- lmer(y ~ temp*wolba + (1|line.id) + (1|line.id:temp))
  res[6:8] <- as.data.frame(VarCorr(lmer.pool.wolba))[, 4]
  res[9:10] <- rand(lmer.pool.wolba)$rand.table$p.value
  
  # return result
  # ============================================================
  
  res
  
}

colnames(gxe) <- c("var.gxe", "var.g", "var.e", "p.g", "p.gxe", "wolba.var.gxe", "wolba.var.g", "wolba.var.e", "wolba.p.g", "wolba.p.gxe")
save(gxe, gene.name, file = args[3])

# session info
# ============================================================

sessionInfo()
