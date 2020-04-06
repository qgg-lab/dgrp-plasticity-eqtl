# ==================================================
# = make a table for female module 1 GO enrichment =
# ==================================================

args <- commandArgs(TRUE) # args <- c("reportData/female.module1.go.RData")

library("openxlsx")
library("GO.db")

# read and process data
# ============================================================
load(args[1])

final.table <- rbind(bp.out, mf.out)

# write excel file
# ============================================================

write.xlsx(final.table, file = args[2], col.names = FALSE, row.names = FALSE)
