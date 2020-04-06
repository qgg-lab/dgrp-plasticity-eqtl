# ==============================================
# = make summary tables for mapping statistics =
# ==============================================
args <- commandArgs(TRUE)
library(xlsx)

# read data
# ============================================================
female.18c <- read.table(args[1], header = FALSE, as.is = TRUE, sep = ":")
male.18c <- read.table(args[2], header = FALSE, as.is = TRUE, sep = ":")
female.25c <- read.table(args[3], header = FALSE, as.is = TRUE, sep = ":")
male.25c <- read.table(args[4], header = FALSE, as.is = TRUE, sep = ":")

# combine data
# ============================================================
final.table <- cbind(female.18c, male.18c[, 2], female.25c[, 2], male.25c[, 2])
for (i in 1:nrow(final.table)) {
  final.table[, 1] <- paste(toupper(substring(final.table[, 1], 1, 1)), substring(final.table[, 1], 2), sep = "")
}
colnames(final.table) <- c("", "Female (18C)", "Male (18C)", "Female (25C)", "Male (25C)")

# write excel file
# ============================================================
write.xlsx(final.table, file = args[5], col.names = TRUE, row.names = FALSE)
