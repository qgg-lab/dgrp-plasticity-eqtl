# =====================================================
# = make a table for the GEI GSEA analysis =
# =====================================================

args <- commandArgs(TRUE) # args <- c("reportData/female.gei.gsea.fdr.RData", "reportData/male.gei.gsea.fdr.RData")
library("openxlsx")
library("GO.db")

# read and process female data
# ============================================================
load(args[1])

female.bp <- unname(cbind("BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                          as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
female.mf <- unname(cbind("MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                          as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
load(args[2])

male.bp <- unname(cbind("BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                        as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
male.mf <- unname(cbind("MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                        as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        



# combine data
# ============================================================

final.table <- rbind(cbind("Female", rbind(as.matrix(female.bp), as.matrix(female.mf))),
                     cbind("Male", rbind(as.matrix(male.bp), as.matrix(male.mf))))

# write excel file
# ============================================================

write.xlsx(final.table, file = args[3], col.names = TRUE, row.names = FALSE)
