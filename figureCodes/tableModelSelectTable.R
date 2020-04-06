# ======================================
# = table for eqtl and model selection =
# ======================================

args <- commandArgs(TRUE) # args <- c("reportData/female.18c.fdr05.eqtl.table.txt", "reportData/female.25c.fdr05.eqtl.table.txt", "reportData/male.18c.fdr05.eqtl.table.txt", "reportData/male.25c.fdr05.eqtl.table.txt", "reportData/gene.info", "reportData/female.18c.genvar.id", "reportData/female.25c.genvar.id", "reportData/male.18c.genvar.id", "reportData/male.25c.genvar.id", "reportData/female.adj.qgGxE.RData", "reportData/male.adj.qgGxE.RData")
library(openxlsx)

# read data
# ============================================================

female.18c <- read.table(args[1], header = FALSE, as.is = TRUE)
female.25c <- read.table(args[2], header = FALSE, as.is = TRUE)
male.18c <- read.table(args[3], header = FALSE, as.is = TRUE)
male.25c <- read.table(args[4], header = FALSE, as.is = TRUE)
rownames(female.18c) <- female.18c[, 1]
rownames(female.25c) <- female.25c[, 1]
rownames(male.18c) <- male.18c[, 1]
rownames(male.25c) <- male.25c[, 1]

gene.info <- read.table(args[5], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

female.18c.genvar <- scan(args[6], what = "", quiet = TRUE)
female.25c.genvar <- scan(args[7], what = "", quiet = TRUE)
male.18c.genvar <- scan(args[8], what = "", quiet = TRUE)
male.25c.genvar <- scan(args[9], what = "", quiet = TRUE)
load(args[10])
female.gxe <- gxe
load(args[11])
male.gxe <- gxe


# combine data
# ============================================================

final.table <- cbind(c(rep("female.18c", nrow(female.18c)), rep("female.25c", nrow(female.25c)),
											 rep("male.18c", nrow(male.18c)), rep("male.25c", nrow(male.25c))),
											 rbind(female.18c[, -3], female.25c[, -3], male.18c[, -3], male.25c[, -3]))
final.table <- cbind(final.table[, 1:2], gene.info[final.table[, 2], ], final.table[, -(1:2)])	

# output table
# ============================================================

cat("number of selected eqtls:\n")
cat(table(sapply(strsplit(final.table[, 6], split = ","), length)), "\n")

cat("number of eqtls mapped among shared gen var genes in females:\n")
cat("18C:", length(intersect(intersect(female.18c.genvar, female.25c.genvar), female.18c[, 1])), "\n")
cat("25C:", length(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1])), "\n")
cat("bothC:", length(intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1])), "\n")
cat("share eqtl:", sum(apply(cbind(female.18c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)], female.25c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0), "\n")

female.gene.share.eqtl <- intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1])[apply(cbind(female.18c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)], female.25c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0]

female.int.fdr <- p.adjust(as.numeric(female.gxe[, 10]), method = "BH")
cat("share eqtl and GxE sig:", length(intersect(gene.name[female.int.fdr < 0.05], female.gene.share.eqtl)), "\n")



cat("number of eqtls mapped among shared gen var genes in males:\n")
cat("18C:", length(intersect(intersect(male.18c.genvar, male.25c.genvar), male.18c[, 1])), "\n")
cat("25C:", length(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1])), "\n")
cat("bothC:", length(intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1])), "\n")
cat("share eqtl:", sum(apply(cbind(male.18c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)], male.25c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0), "\n")

male.gene.share.eqtl <- intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1])[apply(cbind(male.18c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)], male.25c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0]

male.int.fdr <- p.adjust(as.numeric(male.gxe[, 10]), method = "BH")
cat("share eqtl and GxE sig:", length(intersect(gene.name[male.int.fdr < 0.05], male.gene.share.eqtl)), "\n")





colnames(final.table) <- c("condition", "FBgn", "symbol", "type", "pos", "selected.eqtl", "eqtl", "pval", "eff")
                   

# write excel file
# ============================================================

write.xlsx(final.table, file = args[12], col.names = TRUE, row.names = FALSE)
