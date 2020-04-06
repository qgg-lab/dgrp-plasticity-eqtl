# ============================================
# = quantitative genetics of gene expression =
# ============================================

# use normal quantiles to transform gene expression
# within each temperature
# ============================================================

args <- commandArgs(TRUE) # args <- c("tilingInt/female.gene.exp.adjust.RData", 8)
rdata.file <- args[1]
n.cpu <- as.numeric(args[2])

library("nlme")
library("doMC")
registerDoMC(n.cpu)
lme.ctrl <- lmeControl(opt = "optim")

# load data
# ============================================================

load(rdata.file)

# find line id and temperature
# ============================================================

unique.id <- factor(paste(line.id, temp, sep = "-"))

# perform QG analysis
# ============================================================

gene.name <- rownames(gene.exp.adj)

var.het <- foreach(i = 1:nrow(gene.exp.adj), .combine = rbind) %dopar% {
  
  y <- unlist(gene.exp.adj[i, ])
  res <- rep(NA, 12)
  
  # fit pooled model, test for significance
  # ============================================================
  
  lme.pool.fit0 <- try(lme(y ~ temp, random = ~ 0 + temp | unique.id, weights = varIdent(form = ~ 1 | temp), control = lme.ctrl))
  
  # single residual variance
  lme.pool.fit1 <- try(lme(y ~ temp, random = ~ 0 + temp | unique.id, control = lme.ctrl))
  
  # single genetic variance
  lme.pool.fit2 <- try(lme(y ~ temp, random = ~ 1 | unique.id, weights = varIdent(form = ~ 1 | temp), control = lme.ctrl))
  
  # also do the full model with wolbachia
  # ============================================================
  
  lme.pool.fit0.wolba <- try(lme(y ~ temp*wolba, random = ~ 0 + temp | unique.id, weights = varIdent(form = ~ 1 | temp), control = lme.ctrl))
  
  # single residual variance
  lme.pool.fit1.wolba <- try(lme(y ~ temp*wolba, random = ~ 0 + temp | unique.id, control = lme.ctrl))
  
  # single genetic variance
  lme.pool.fit2.wolba <- try(lme(y ~ temp*wolba, random = ~ 1 | unique.id, weights = varIdent(form = ~ 1 | temp), control = lme.ctrl))
  
  # save results
  # ============================================================

  if (attr(lme.pool.fit0, "class") == "lme" & attr(lme.pool.fit1, "class") == "lme" & attr(lme.pool.fit2, "class") == "lme") {
    
    res[1:6] <- c(VarCorr(lme.pool.fit0)[, 1],
                   as.numeric(VarCorr(lme.pool.fit0)[3, 1])*exp(summary(lme.pool.fit0)$modelStruct$varStruct[[1]])^2,
                   pchisq(2*(logLik(lme.pool.fit0) - logLik(lme.pool.fit1)), df = 1, lower.tail = FALSE),
                   pchisq(2*(logLik(lme.pool.fit0) - logLik(lme.pool.fit2)), df = 1, lower.tail = FALSE))
  }
    
  if (attr(lme.pool.fit0.wolba, "class") == "lme" & attr(lme.pool.fit1.wolba, "class") == "lme" & attr(lme.pool.fit2.wolba, "class") == "lme") {
    
    res[7:12] <- c(VarCorr(lme.pool.fit0.wolba)[, 1],
                   as.numeric(VarCorr(lme.pool.fit0.wolba)[3, 1])*exp(summary(lme.pool.fit0.wolba)$modelStruct$varStruct[[1]])^2,
                   pchisq(2*(logLik(lme.pool.fit0.wolba) - logLik(lme.pool.fit1.wolba)), df = 1, lower.tail = FALSE),
                   pchisq(2*(logLik(lme.pool.fit0.wolba) - logLik(lme.pool.fit2.wolba)), df = 1, lower.tail = FALSE))
  }
  
  cat(gene.name[i], "\n")
  
  res
    
} 

colnames(var.het) <- c("var.25c", "var.18c", "var.e.18c", "var.e.25c", "p.single.e", "p.single.g", "wolba.var.25c", "wolba.var.18c", "wolba.var.e.18c", "wolba.var.e.25c", "wolba.p.single.e", "wolba.p.single.g")
save(var.het, gene.name, file = args[3])

# session info
# ============================================================

sessionInfo()
