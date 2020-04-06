# =========================
# = eqtl cross prediction =
# =========================

args <- commandArgs(TRUE) # args <- c("eqtl/female.18c.model.select.RData", "eqtl/female.18c.genvar.id", "eqtl/female.25c.model.select.RData", "eqtl/female.25c.genvar.id", "eqtl/eqtl.snp.geno.tped", "eqtl/eqtl.snp.geno.tfam", "eqtl/female.18c.pheno", "eqtl/female.25c.pheno")

# load data
# ============================================================

load(args[1])
eqtl.18c.pair <- mc.select.eqtl
rownames(eqtl.18c.pair) <- eqtl.18c.pair[, 1]
eqtl.18c.genvar <- scan(args[2], what = "")

load(args[3])
eqtl.25c.pair <- mc.select.eqtl
rownames(eqtl.25c.pair) <- eqtl.25c.pair[, 1]
eqtl.25c.genvar <- scan(args[4], what = "")

common.genvar <- intersect(eqtl.18c.genvar, eqtl.25c.genvar)


tfam <- read.table(args[6], header = FALSE, as.is = TRUE)
n.line <- nrow(tfam)
tped <- matrix(scan(args[5], what = ""), ncol = 4 + 2*n.line, byrow = TRUE)

geno.code <- matrix(4 - (as.numeric(tped[, seq(from = 5, length = n.line, by = 2)]) + as.numeric(tped[, seq(from = 6, length = n.line, by = 2)])), ncol = n.line)
geno.code[geno.code == 4] <- NA
rownames(geno.code) <- tped[, 2]

# read phenotype
pheno.18c <- read.table(args[7], header = TRUE, as.is = TRUE, row.names = 1)
pheno.18c <- pheno.18c[tfam[, 1], ]

pheno.25c <- read.table(args[8], header = TRUE, as.is = TRUE, row.names = 1)
pheno.25c <- pheno.25c[tfam[, 1], ]

eqtl.pred.18c <- matrix(nrow = length(common.genvar), ncol = 4)

for (i in 1:length(common.genvar)) {
  gene <- common.genvar[i]
  eqtl.pred.18c[i, 1] <- gene
  eqtl.pred.18c[i, 2] <- cor(pheno.18c[, gene], pheno.25c[, gene])
  if (gene %in% eqtl.18c.pair[, 1]) {
    snps <- gsub("VAR_", "", unlist(strsplit(eqtl.18c.pair[gene, 2], split = ",")))
    lm.fit.18c <- lm(pheno.18c[, gene] ~ t(matrix(geno.code[snps, ], nrow = length(snps))))
    fitted <- (cbind(1, t(matrix(geno.code[snps, ], nrow = length(snps)))) %*% coefficients(lm.fit.18c))[, 1]
    eqtl.pred.18c[i, 3] <- cor(unlist(pheno.18c[, gene]), fitted, use = "pairwise")
    eqtl.pred.18c[i, 4] <- cor(unlist(pheno.25c[, gene]), fitted, use = "pairwise")
  }
  
}


eqtl.pred.25c <- matrix(nrow = length(common.genvar), ncol = 4)

for (i in 1:length(common.genvar)) {
  gene <- common.genvar[i]
  eqtl.pred.25c[i, 1] <- gene
  eqtl.pred.25c[i, 2] <- cor(pheno.25c[, gene], pheno.18c[, gene])
  if (gene %in% eqtl.25c.pair[, 1]) {
    snps <- gsub("VAR_", "", unlist(strsplit(eqtl.25c.pair[gene, 2], split = ",")))
    lm.fit.25c <- lm(pheno.25c[, gene] ~ t(matrix(geno.code[snps, ], nrow = length(snps))))
    fitted <- (cbind(1, t(matrix(geno.code[snps, ], nrow = length(snps)))) %*% coefficients(lm.fit.25c))[, 1]
    eqtl.pred.25c[i, 3] <- cor(unlist(pheno.25c[, gene]), fitted, use = "pairwise")
    eqtl.pred.25c[i, 4] <- cor(unlist(pheno.18c[, gene]), fitted, use = "pairwise")
  }
  
}


save(eqtl.pred.18c, eqtl.pred.25c, file = args[9])
