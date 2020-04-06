# ==================================
# = table for network preservation =
# ==================================

args <- commandArgs(TRUE) # args <- c("reportData/female.wgcna.RData", "reportData/male.wgcna.RData")

library("openxlsx")

# read and process data
# ============================================================

load(args[1])
female.mod.pres <- mod.pres

load(args[2])
male.mod.pres <- mod.pres

# final.table
# ============================================================

final.table <- rbind(cbind("female", 1:nrow(female.mod.pres$observed), female.mod.pres$observed, female.mod.pres$p.values),
										 cbind("male", 1:nrow(male.mod.pres$observed), male.mod.pres$observed, male.mod.pres$p.values))

# write excel file
# ============================================================

write.xlsx(final.table, file = args[3], col.names = TRUE, row.names = FALSE)
