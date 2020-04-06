# ==========================
# = estimate empirical fdr =
# ==========================

args <- commandArgs(TRUE) # args <- c("eqtl/female.18c.pheno", "eqtl/female.18c/", "eqtl/perm/female.18c.perm", 0.20, 100, 5, "test.out")

pheno.file <- args[1]
obs.dir <- args[2]
perm.dir <- args[3]
fdr <- as.numeric(args[4])
n.perm <- as.numeric(args[5])
n.cpu <- as.numeric(args[6])
out.file <- args[7]
library(doMC)
registerDoMC(n.cpu)

# phenotype file
# ============================================================

gene.list <- colnames(read.table(pheno.file, header = TRUE, as.is = TRUE, row.names = 1)[, -1])

fdr.out <- foreach (i = 1:length(gene.list), .combine = rbind) %dopar% {
  
  # first read the observed p value
  
  obs.pval <- read.table(paste(obs.dir, "/eqtl.", gene.list[i], ".qassoc", sep = ""),
                         header = TRUE, as.is = TRUE)
  res <- c(gene.list[i], fdr, NA, NA, NA)
  
  if (nrow(obs.pval) > 0) {
    obs.pval <- obs.pval[order(obs.pval[, 9]), c(2, 9, 5)]
    # read permuted data
    perm.pval <- data.frame(matrix(scan(con <- pipe(paste("cat ", perm.dir, "*/eqtl.perm.", gene.list[i], ".qassoc | awk \'$1 != \"CHR\" {print $2\"\\t\"$9}\'", sep = "")), 
                                        what = ""), ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
    perm.pval[, 2] <- as.numeric(perm.pval[, 2])
    obs.pval <- cbind(obs.pval, rep(0, nrow(obs.pval)))
    for (j in 1:nrow(obs.pval)) {
      obs.pval[j, 4] <- sum(perm.pval[, 2] <= obs.pval[j, 2])/n.perm/j
    }
    # identify the cutoff point
    thres <- which(cummax(obs.pval[, 4]) < fdr)
    if (length(thres) > 0) {
      res <- c(gene.list[i], obs.pval[max(thres), 4], paste(obs.pval[1:max(thres), 1], collapse = ","), paste(obs.pval[1:max(thres), 2], collapse = ","),  paste(obs.pval[1:max(thres), 3], collapse = ","))
    } 
  }
  
  cat(gene.list[i], "\n")
  res
}

write.table(fdr.out, file = out.file, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

sessionInfo()
