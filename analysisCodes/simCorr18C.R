# ==============================================================
# = simulate expression at 18C and caculate correlation matrix =
# ==============================================================

args <- commandArgs(TRUE)
load(args[1])

# read data
# ============================================================

exp.18c <- t(blup.18c.data[match(intersect(sig.18c.gene, sig.25c.gene), sig.gene.name), ])
colnames(exp.18c) <- intersect(sig.18c.gene, sig.25c.gene)
exp.25c <- t(blup.25c.data[match(intersect(sig.18c.gene, sig.25c.gene), sig.gene.name), ])
colnames(exp.25c) <- intersect(sig.18c.gene, sig.25c.gene)

# function to simulate new expression
# ============================================================

sim.exp <- function(exp1, exp2) {
  
  exp.cor <- cor(exp1, exp2)
  return(sd(exp2)*(exp.cor*scale(exp1) + sqrt(1-exp.cor^2)*rnorm(length(exp2))) + mean(exp2))
  
}

# main program
# ============================================================

set.seed(as.numeric(args[2]))

for (i in 1:ncol(exp.18c)) {
  
  exp.18c[, i] <- sim.exp(exp.25c[, i], exp.18c[, i])
  
}

sim.corr <- cor(exp.18c)
sim.corr[lower.tri(sim.corr)] <- 0

save(sim.corr, file = args[3])
