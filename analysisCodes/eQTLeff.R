# =========================
# = eqtl effect estimates =
# =========================

args <- commandArgs(TRUE) # args <- c("eqtl/female.18c.model.select.RData", "eqtl/female.18c.eqtl.fdr.out", "eqtl/female.18c.genvar.id", "eqtl/female.25c.model.select.RData", "eqtl/female.25c.eqtl.fdr.out", "eqtl/female.25c.genvar.id", "eqtl/eqtl.snp.geno.tped", "eqtl/eqtl.snp.geno.tfam", "eqtl/female.18c.pheno", "eqtl/female.25c.pheno")

# load data
# ============================================================

load(args[1])
eqtl.18c.pair <- mc.select.eqtl
rownames(eqtl.18c.pair) <- eqtl.18c.pair[, 1]
eqtl.18c.fdr <- read.table(args[2], header = FALSE, as.is = TRUE)
rownames(eqtl.18c.fdr) <- eqtl.18c.fdr[, 1]
eqtl.18c.genvar <- scan(args[3], what = "")

load(args[4])
eqtl.25c.pair <- mc.select.eqtl
rownames(eqtl.25c.pair) <- eqtl.25c.pair[, 1]
eqtl.25c.fdr <- read.table(args[5], header = FALSE, as.is = TRUE)
rownames(eqtl.25c.fdr) <- eqtl.25c.fdr[, 1]
eqtl.25c.genvar <- scan(args[6], what = "")

common.genvar <- intersect(eqtl.18c.genvar, eqtl.25c.genvar)

eqtl.pair <- NULL

for (i in 1:length(common.genvar)) {
  
  gene <- common.genvar[i]
  
  # find unique pairs
  if (gene %in% eqtl.18c.pair[, 1]) {
    gene.eqtl.18c <- gsub("VAR_", "", unlist(strsplit(eqtl.18c.pair[gene, 2], split = ",")))
    gene.fdr.18c <- unlist(strsplit(eqtl.18c.fdr[gene, 3], split = ","))
    if (gene %in% eqtl.25c.pair[, 1]) {
      gene.eqtl.25c <- gsub("VAR_", "", unlist(strsplit(eqtl.25c.pair[gene, 2], split = ",")))
      gene.fdr.25c <- unlist(strsplit(eqtl.25c.fdr[gene, 3], split = ","))
      uniq.eqtl <- unique(c(gene.eqtl.18c, gene.eqtl.25c))
      for (eqtl in uniq.eqtl) {
        if (eqtl %in% gene.fdr.18c) {
          if (eqtl %in% gene.fdr.25c) {
            eqtl.pair <- rbind(eqtl.pair, c(gene, eqtl, "both"))
          } else {
            eqtl.pair <- rbind(eqtl.pair, c(gene, eqtl, "18c"))
          }
        } else {
          if (eqtl %in% gene.fdr.25c) {
            eqtl.pair <- rbind(eqtl.pair, c(gene, eqtl, "25c"))
          }
        }
      }
    } else {
      for (eqtl in gene.eqtl.18c) {
        eqtl.pair <- rbind(eqtl.pair, c(gene, eqtl, "18c"))
      }
    }
  } else {
    if (gene %in% eqtl.25c.pair[, 1]) {
      gene.eqtl.25c <- gsub("VAR_", "", unlist(strsplit(eqtl.25c.pair[gene, 2], split = ",")))
      gene.fdr.25c <- unlist(strsplit(eqtl.25c.fdr[gene, 3], split = ","))
      for (eqtl in gene.eqtl.25c) {
        eqtl.pair <- rbind(eqtl.pair, c(gene, eqtl, "25c"))
      }
    }
  }
}       

tfam <- read.table(args[8], header = FALSE, as.is = TRUE)
n.line <- nrow(tfam)
tped <- matrix(scan(args[7], what = ""), ncol = 4 + 2*n.line, byrow = TRUE)

geno.code <- matrix(4 - (as.numeric(tped[, seq(from = 5, length = n.line, by = 2)]) + as.numeric(tped[, seq(from = 6, length = n.line, by = 2)])), ncol = n.line)
geno.code[geno.code == 4] <- NA
rownames(geno.code) <- tped[, 2]

# read phenotype
pheno.18c <- read.table(args[9], header = TRUE, as.is = TRUE, row.names = 1)
pheno.18c <- pheno.18c[tfam[, 1], ]

pheno.25c <- read.table(args[10], header = TRUE, as.is = TRUE, row.names = 1)
pheno.25c <- pheno.25c[tfam[, 1], ]

# estimate effects
# ============================================================

eqtl.pair <- cbind(eqtl.pair, NA, NA, NA, NA, NA, NA)
col.names.18c <- colnames(pheno.18c)
col.names.25c <- colnames(pheno.25c)

for (i in 1:nrow(eqtl.pair)) {
  # some genes don't exist in one temperature because they are not genetically variable, need to check these
  if (sum(col.names.18c == eqtl.pair[i, 1]) > 0) {
    eqtl.pair[i, 4:5] <- summary(lm(pheno.18c[, eqtl.pair[i, 1]] ~ geno.code[eqtl.pair[i, 2], ]))$coefficients[2, 1:2]
    eqtl.pair[i, 6] <- var(pheno.18c[, eqtl.pair[i, 1]])
  }
  if (sum(col.names.25c == eqtl.pair[i, 1]) > 0) {
    eqtl.pair[i, 7:8] <- summary(lm(pheno.25c[, eqtl.pair[i, 1]] ~ geno.code[eqtl.pair[i, 2], ]))$coefficients[2, 1:2]
    eqtl.pair[i, 9] <- var(pheno.25c[, eqtl.pair[i, 1]])
  }
}

save(eqtl.pair, file = args[11])

