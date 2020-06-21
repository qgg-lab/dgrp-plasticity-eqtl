# ================================================
# = overlap between GxE genes with other studies =
# ================================================

args <- commandArgs(TRUE) # args <- c("reportData/cg.fbgn.txt", "reportData/female.adj.qgGxE.RData", "reportData/male.adj.qgGxE.RData", "reportData/gene.info", "reportData/Hutter2008CGs.txt", "reportData/Levine2011CGs.txt", "reportData/Zhao2015dmel21cFBgns.txt", "reportData/Zhao2015dmel29cFBgns.txt", "report/Table_GeneComparison.xlsx")
library(xlsx)

# read data
# ============================================================
cg.fbgn <- read.table(args[1], sep = "\t", as.is = TRUE, header = FALSE)
load(args[2])
female.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)
female.gene.name <- gene.name
load(args[3])
male.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)
male.gene.name <- gene.name
gene.info <- read.table(args[4], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)
hutter2008 <- scan(args[5], what = "") # male expression
levine2011 <- scan(args[6], what = "") # male expression
zhao2015c21 <- scan(args[7], what = "")
zhao2015c29 <- scan(args[8], what = "")
# gxe gene
# ============================================================

female.gei <- female.gxe[, 6]/(female.gxe[, 6] + female.gxe[, 7])
female.int.fdr <- p.adjust(as.numeric(female.gxe[, 10]), method = "BH")
female.gei.gene <- female.gene.name[female.int.fdr < 0.05]

male.gei <- male.gxe[, 6]/(male.gxe[, 6] + male.gxe[, 7])
male.int.fdr <- p.adjust(as.numeric(male.gxe[, 10]), method = "BH")
male.gei.gene <- male.gene.name[which(male.int.fdr < 0.05)]

# comparison with hutter genes, overlap / with annotation / original hutter gene count
# Hutter et al. 2008 additional file 4
# male expression
# ============================================================

cat("Hutter: ", length(intersect(cg.fbgn[match(hutter2008, cg.fbgn[, 1]), 2], male.gei.gene)), "/", length(intersect(cg.fbgn[match(hutter2008, cg.fbgn[, 1]), 2], male.gene.name)), "/", length(hutter2008), "\n")

# comparison with levine genes, Levine et al 2011 Table s1c, GEI genes
# male expression
# ============================================================

cat("Levine: ", length(intersect(cg.fbgn[match(levine2011, cg.fbgn[, 1]), 2], male.gei.gene)), "/", length(intersect(cg.fbgn[match(levine2011, cg.fbgn[, 1]), 2], male.gene.name)), "/", length(levine2011), "\n")

# comparison with Zhao et al 2015 gene, table s2 (sheet3) at 21 degrees
# ============================================================

cat("Zhao 21c: ", length(intersect(zhao2015c21, male.gei.gene)), "/", length(intersect(zhao2015c21, male.gene.name)), "/", length(zhao2015c21), "\n")

# comparison with Zhao et al 2015 gene, table s2 (sheet2) at 29 degrees
# ============================================================

cat("Zhao 29c: ", length(intersect(zhao2015c29, male.gei.gene)), "/", length(intersect(zhao2015c29, male.gene.name)), "/", length(zhao2015c29), "\n")

# all genes
# ============================================================

cat("Combined: ", length(intersect(male.gei.gene, unique(c(cg.fbgn[match(hutter2008, cg.fbgn[, 1]), 2], cg.fbgn[match(levine2011, cg.fbgn[, 1]), 2], zhao2015c21, zhao2015c29)))), "\n")

# generate data for table
all.gene.table <- NULL
# 1. unique set of all genes
all.gene <- unique(c(intersect(cg.fbgn[match(hutter2008, cg.fbgn[, 1]), 2], male.gene.name), intersect(cg.fbgn[match(levine2011, cg.fbgn[, 1]), 2], male.gene.name), intersect(zhao2015c21, male.gene.name), intersect(zhao2015c29, male.gene.name)))
all.gene.table <- cbind(all.gene, gene.info[all.gene, ])

# 2. gei annotation
all.gene.table <- cbind(all.gene.table, ifelse(all.gene %in% male.gei.gene, "x", ""),
												ifelse(all.gene %in% cg.fbgn[match(hutter2008, cg.fbgn[, 1]), 2], "x", ""),
												ifelse(all.gene %in% cg.fbgn[match(levine2011, cg.fbgn[, 1]), 2], "x", ""),
												ifelse(all.gene %in% zhao2015c21, "x", ""),
												ifelse(all.gene %in% zhao2015c29, "x", ""))
all.gene.table <- all.gene.table[order(all.gene.table[, 1]), ]

# write excel file
# ============================================================
write.xlsx(all.gene.table, file = args[9], col.names = TRUE, row.names = FALSE)
