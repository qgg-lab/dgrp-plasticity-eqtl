# ===============================================
# = summarize preliminary normalization results =
# ===============================================

# load data
# ============================================================

load("tilingInt/prelim.norm.gene.exp.RData")

# scale data
# ============================================================

female.18c.gene.exp <- t(scale(t(female.18c.gene.exp)))
female.25c.gene.exp <- t(scale(t(female.25c.gene.exp)))

male.18c.gene.exp <- t(scale(t(male.18c.gene.exp)))
male.25c.gene.exp <- t(scale(t(male.25c.gene.exp)))

# summarize quantiles
# ============================================================

female.18c.quantile <- apply(female.18c.gene.exp, 2, function(x) { quantile(x, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) })
female.25c.quantile <- apply(female.25c.gene.exp, 2, function(x) { quantile(x, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) })
male.18c.quantile <- apply(male.18c.gene.exp, 2, function(x) { quantile(x, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) })
male.25c.quantile <- apply(male.25c.gene.exp, 2, function(x) { quantile(x, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) })

# save data
# ============================================================

save(female.18c.quantile, male.18c.quantile, female.25c.quantile, male.25c.quantile, file = "tilingInt/prelim.norm.summary.RData")
sessionInfo()
