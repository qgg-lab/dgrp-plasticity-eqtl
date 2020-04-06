# ============================================
# = summarize correlation simulation results =
# ============================================

args <- commandArgs(TRUE) # args <- c("qg/female.sig.gene.exp.RData", "qg/female.adj.qgGxE.RData", 10, "simCorr/female.sim.corr.18c.perm", 1000)

load(args[1])
load(args[2])
library("doMC")
registerDoMC(as.numeric(args[3]))
sim.file.prefix <- args[4]
out.file <- args[6]
n.sim <- as.numeric(args[5])

# read data
# ============================================================

exp.18c <- t(blup.18c.data[match(intersect(sig.18c.gene, sig.25c.gene), sig.gene.name), ])
colnames(exp.18c) <- intersect(sig.18c.gene, sig.25c.gene)
exp.25c <- t(blup.25c.data[match(intersect(sig.18c.gene, sig.25c.gene), sig.gene.name), ])
colnames(exp.25c) <- intersect(sig.18c.gene, sig.25c.gene)

# gxe criteria
# ============================================================

int.fdr <- p.adjust(as.numeric(gxe[, 10]), method = "BH")
int.gene <- na.omit(match(gene.name[which(int.fdr < 0.05)], colnames(exp.25c)))

# observed correlation and difference
# ============================================================

corr.25c <- cor(exp.25c)
corr.25c <- corr.25c[c(setdiff(1:ncol(exp.25c), int.gene), int.gene), c(setdiff(1:ncol(exp.25c), int.gene), int.gene)]
corr.25c.gxe <- corr.25c[(ncol(exp.25c) - length(int.gene) + 1):ncol(exp.25c), (ncol(exp.25c) - length(int.gene) + 1):ncol(exp.25c)]
corr.25c.g <- corr.25c[1:(ncol(exp.25c) - length(int.gene)), 1:(ncol(exp.25c) - length(int.gene))]
corr.25c.bt <- corr.25c[1:(ncol(exp.25c) - length(int.gene)), (ncol(exp.25c) - length(int.gene) + 1):ncol(exp.25c)]

corr.18c <- cor(exp.18c)
corr.18c <- corr.18c[c(setdiff(1:ncol(exp.18c), int.gene), int.gene), c(setdiff(1:ncol(exp.18c), int.gene), int.gene)]
corr.18c.gxe <- corr.18c[(ncol(exp.18c) - length(int.gene) + 1):ncol(exp.18c), (ncol(exp.18c) - length(int.gene) + 1):ncol(exp.18c)]
corr.18c.g <- corr.18c[1:(ncol(exp.18c) - length(int.gene)), 1:(ncol(exp.18c) - length(int.gene))]
corr.18c.bt <- corr.18c[1:(ncol(exp.18c) - length(int.gene)), (ncol(exp.18c) - length(int.gene) + 1):ncol(exp.18c)]

# correlation differece by class
# ============================================================

corr.diff.gxe.dens <- density(corr.18c.gxe[upper.tri(corr.18c.gxe, diag = FALSE)] - corr.25c.gxe[upper.tri(corr.25c.gxe, diag = FALSE)], from = -2, to = 2, n = 512)
corr.diff.gxe.dens$mean <- mean(corr.18c.gxe[upper.tri(corr.18c.gxe, diag = FALSE)] - corr.25c.gxe[upper.tri(corr.25c.gxe, diag = FALSE)])
corr.diff.gxe.dens$sd <- sd(corr.18c.gxe[upper.tri(corr.18c.gxe, diag = FALSE)] - corr.25c.gxe[upper.tri(corr.25c.gxe, diag = FALSE)])

corr.diff.g.dens <- density(corr.18c.g[upper.tri(corr.18c.g, diag = FALSE)] - corr.25c.g[upper.tri(corr.25c.g, diag = FALSE)], from = -2, to = 2, n = 512)
corr.diff.g.dens$mean <- mean(corr.18c.g[upper.tri(corr.18c.g, diag = FALSE)] - corr.25c.g[upper.tri(corr.25c.g, diag = FALSE)])
corr.diff.g.dens$sd <- sd(corr.18c.g[upper.tri(corr.18c.g, diag = FALSE)] - corr.25c.g[upper.tri(corr.25c.g, diag = FALSE)])

corr.diff.bt.dens <- density(corr.18c.bt - corr.25c.bt, from = -2, to = 2, n = 512)
corr.diff.bt.dens$mean <- mean(corr.18c.bt - corr.25c.bt)
corr.diff.bt.dens$sd <- sd(corr.18c.bt - corr.25c.bt)

corr.diff.dens <- density(corr.18c[upper.tri(corr.18c, diag = FALSE)] - corr.25c[upper.tri(corr.25c, diag = FALSE)], from = -2, to = 2, n = 512)
corr.diff.dens$mean <- mean(corr.18c[upper.tri(corr.18c, diag = FALSE)] - corr.25c[upper.tri(corr.25c, diag = FALSE)])
corr.diff.dens$sd <- sd(corr.18c[upper.tri(corr.18c, diag = FALSE)] - corr.25c[upper.tri(corr.25c, diag = FALSE)])

# simulated correlation difference
# ============================================================

sim.corr.diff.gxe.dens <- foreach(i = 1:n.sim) %dopar% {
  
  load(paste(sim.file.prefix, i, ".RData", sep = ""))
  sim.corr <- sim.corr[c(setdiff(1:ncol(exp.25c), int.gene), int.gene), c(setdiff(1:ncol(exp.25c), int.gene), int.gene)]
  sim.corr.gxe <- sim.corr[(ncol(exp.18c) - length(int.gene) + 1):ncol(exp.18c), (ncol(exp.18c) - length(int.gene) + 1):ncol(exp.18c)]
  sim.corr.diff.gxe.dens.sim <- density(sim.corr.gxe[upper.tri(sim.corr.gxe, diag = FALSE)] - corr.25c.gxe[upper.tri(corr.25c.gxe, diag = FALSE)], from = -2, to = 2, n = 512)
  sim.corr.diff.gxe.dens.sim$mean <- mean(sim.corr.gxe[upper.tri(sim.corr.gxe, diag = FALSE)] - corr.25c.gxe[upper.tri(corr.25c.gxe, diag = FALSE)])
	sim.corr.diff.gxe.dens.sim$sd <- sd(sim.corr.gxe[upper.tri(sim.corr.gxe, diag = FALSE)] - corr.25c.gxe[upper.tri(corr.25c.gxe, diag = FALSE)])
  cat(i, "\n")
	sim.corr.diff.gxe.dens.sim
  
}

sim.corr.diff.g.dens <- foreach(i = 1:n.sim) %dopar% {
  
  load(paste(sim.file.prefix, i, ".RData", sep = ""))
  sim.corr <- sim.corr[c(setdiff(1:ncol(exp.25c), int.gene), int.gene), c(setdiff(1:ncol(exp.25c), int.gene), int.gene)]
  sim.corr.g <- sim.corr[1:(ncol(exp.18c) - length(int.gene)), 1:(ncol(exp.18c) - length(int.gene))]
  sim.corr.diff.g.dens.sim <- density(sim.corr.g[upper.tri(sim.corr.g, diag = FALSE)] - corr.25c.g[upper.tri(corr.25c.g, diag = FALSE)], from = -2, to = 2, n = 512)
	sim.corr.diff.g.dens.sim$mean <- mean(sim.corr.g[upper.tri(sim.corr.g, diag = FALSE)] - corr.25c.g[upper.tri(corr.25c.g, diag = FALSE)])
	sim.corr.diff.g.dens.sim$sd <- sd(sim.corr.g[upper.tri(sim.corr.g, diag = FALSE)] - corr.25c.g[upper.tri(corr.25c.g, diag = FALSE)])
  cat(i, "\n")
  sim.corr.diff.g.dens.sim
  
}

sim.corr.diff.bt.dens <- foreach(i = 1:n.sim) %dopar% {
  
  load(paste(sim.file.prefix, i, ".RData", sep = ""))
  sim.corr <- sim.corr[c(setdiff(1:ncol(exp.25c), int.gene), int.gene), c(setdiff(1:ncol(exp.25c), int.gene), int.gene)]
  sim.corr.bt <- sim.corr[1:(ncol(exp.18c) - length(int.gene)), (ncol(exp.18c) - length(int.gene) + 1):ncol(exp.18c)]
  sim.corr.diff.bt.dens.sim <- density(sim.corr.bt - corr.25c.bt, from = -2, to = 2, n = 512)
	sim.corr.diff.bt.dens.sim$mean <- mean(sim.corr.bt - corr.25c.bt)
	sim.corr.diff.bt.dens.sim$sd <- sd(sim.corr.bt - corr.25c.bt)
  cat(i, "\n")
  sim.corr.diff.bt.dens.sim
  
}

sim.corr.diff.dens <- foreach(i = 1:n.sim) %dopar% {
  
  load(paste(sim.file.prefix, i, ".RData", sep = ""))
  sim.corr <- sim.corr[c(setdiff(1:ncol(exp.25c), int.gene), int.gene), c(setdiff(1:ncol(exp.25c), int.gene), int.gene)]
  sim.corr.diff.dens.sim <- density(sim.corr - corr.25c, from = -2, to = 2, n = 512)
	sim.corr.diff.dens.sim$mean <- mean(sim.corr - corr.25c)
	sim.corr.diff.dens.sim$sd <- sd(sim.corr - corr.25c)
  cat(i, "\n")
  sim.corr.diff.dens.sim
  
}


save(corr.diff.gxe.dens, corr.diff.g.dens, corr.diff.bt.dens, corr.diff.dens, sim.corr.diff.gxe.dens, sim.corr.diff.g.dens, sim.corr.diff.bt.dens, sim.corr.diff.dens, file = out.file)
