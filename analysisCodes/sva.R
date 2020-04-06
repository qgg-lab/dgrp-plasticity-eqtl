# ==================
# = sva correction =
# ==================

args <- commandArgs(TRUE) # args <- c("~/Downloads/female.gene.exp.RData", "female.gene.exp", "~/work/projects/geneExpTempEff/reportData/cel.headers.RData", "~/dgrp/adjustData.RData")
rdata.file <- args[1]
data.set <- args[2]
library("lubridate")
library("sva")

# load data
# ============================================================

load(rdata.file)
eval(parse(text = paste("gene.exp = ", data.set, sep = "")))
load(args[3])
load(args[4])

nqt <- function(x) {
  
  return(qnorm(rank(x)/(length(x) + 1), sd = mad(x)) + median(x))
  
}

# date
# ============================================================

rownames(cel.headers.18c) <- cel.headers.18c[, 1]
rownames(cel.headers.25c) <- cel.headers.25c[, 1]

cel.date.18c <- parse_date_time(cel.headers.18c[, 2], "mdy HMS")
cel.date.25c <- parse_date_time(cel.headers.25c[, 2], "mdy HMS")

names(cel.date.18c) <- cel.headers.18c[, 1]
names(cel.date.25c) <- cel.headers.25c[, 1]

# find line id and temperature
# ============================================================

temp <- rep("25C", ncol(gene.exp))
temp[grep("18C", colnames(gene.exp))] <- "18C"
temp <- factor(temp, levels = c("25C", "18C"))
line.id <- factor(paste("line_", gsub("T", "", gsub("-.*", "", colnames(gene.exp))), sep = ""))
scan.date <- date(c(cel.date.18c, cel.date.25c)[colnames(gene.exp)])
wolba <- factor(wolba[as.character(line.id), 1])

# estimate sv number
# ============================================================
gene.exp.nqt <- gene.exp
for (i in 1:nrow(gene.exp)) { gene.exp.nqt[i, ] <- nqt(gene.exp[i, ]) }

mod = model.matrix( ~ temp + wolba)
mod0 = matrix(1, nrow = ncol(gene.exp.nqt), ncol = 1)
colnames(mod0) <- "Intercept"
nsv <- num.sv(gene.exp.nqt, mod = mod, method = "be", B = 20, seed = 1)

sva.adj <- sva(gene.exp.nqt, mod, mod0, n.sv = nsv)

gene.exp.adj <- gene.exp.nqt
for (i in 1:nrow(gene.exp.nqt)) {
  this.fit <- lm(gene.exp.nqt[i, ] ~ sva.adj$sv)
  gene.exp.adj[i, ] <- coefficients(this.fit)[1] + residuals(this.fit)
}

save(sva.adj, scan.date, file = args[5])
save(gene.exp.adj, wolba, temp, line.id, file = args[6])

# session info
# ============================================================

sessionInfo()
