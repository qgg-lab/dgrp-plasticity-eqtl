# ==================================================
# = figure for eQTL mapping summary, effects, etc. =
# ==================================================

args <- commandArgs(TRUE) # args <- c("reportData/female.18c.fdr05.eqtl.table.txt", "reportData/female.25c.fdr05.eqtl.table.txt", "reportData/male.18c.fdr05.eqtl.table.txt", "reportData/male.25c.fdr05.eqtl.table.txt", "reportData/gene.info", "reportData/female.18c.genvar.id", "reportData/female.25c.genvar.id", "reportData/male.18c.genvar.id", "reportData/male.25c.genvar.id", "reportData/female.adj.qgGxE.RData", "reportData/male.adj.qgGxE.RData", "reportData/female.eqtl.eff.RData", "reportData/male.eqtl.eff.RData", "reportData/female.eqtl.pred.RData", "reportData/male.eqtl.pred.RData", "reportData/female.tf.summary.RData", "reportData/male.tf.summary.RData")
library("RColorBrewer")

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
rownames(female.gxe) <- gene.name
load(args[11])
male.gxe <- gxe
rownames(male.gxe) <- gene.name

load(args[12])
female.eqtl.pair <- eqtl.pair
load(args[13])
male.eqtl.pair <- eqtl.pair

load(args[14])
female.eqtl.pred.25c <- eqtl.pred.25c
female.eqtl.pred.18c <- eqtl.pred.18c

load(args[15])
male.eqtl.pred.25c <- eqtl.pred.25c
male.eqtl.pred.18c <- eqtl.pred.18c

load(args[16])

female.tf.fold.enrich <- tf.fold.enrich[tf.list != "mE1_TFBS_HSA", ]
female.tf.fold.enrich.p <- tf.fold.enrich.p[tf.list != "mE1_TFBS_HSA", ]

load(args[17])
male.tf.fold.enrich <- tf.fold.enrich[tf.list != "mE1_TFBS_HSA", ]
male.tf.fold.enrich.p <- tf.fold.enrich.p[tf.list != "mE1_TFBS_HSA", ]

tf.list <- tf.list[tf.list != "mE1_TFBS_HSA"]

# set up file
# ============================================================

file.width = 183
cairo_pdf(file = args[18], width = file.width/25.4, height = file.width*0.45/25.4, family = args[19])
par(las = 1, tcl = -0.2, ps = 7, lwd = 0.5, mfrow = c(3, 2), xpd = TRUE)
layout(mat = matrix(c(1, 2, 3, 4, 1, 2, 5, 6, 7, 7, 5, 6), ncol = 4, byrow = T), height = c(1, 0.3, 0.7))

# ribbons in females
# ============================================================

left.count <- numeric(6)
right.count <- numeric(6)

# 25c, 
left.count[1] <- length(intersect(female.25c[, 1], female.18c[, 1]))
left.count[2] <- length(setdiff(intersect(intersect(female.25c.genvar, female.18c.genvar), female.25c[, 1]), female.18c[, 1]))
left.count[3] <- length(setdiff(intersect(intersect(female.25c.genvar, female.18c.genvar), female.18c[, 1]), female.25c[, 1]))
left.count[4] <- length(setdiff(setdiff(intersect(female.25c.genvar, female.18c.genvar), female.18c[, 1]), female.25c[, 1]))
left.count[5] <- length(intersect(female.25c[, 1], setdiff(female.25c.genvar, female.18c.genvar)))
left.count[6] <- length(setdiff(setdiff(female.25c.genvar, female.18c.genvar), female.25c[, 1]))

# 18c
right.count[1] <- left.count[1]
right.count[2] <- left.count[3]
right.count[3] <- left.count[2]
right.count[4] <- left.count[4]
right.count[5] <- length(intersect(female.18c[, 1], setdiff(female.18c.genvar, female.25c.genvar)))
right.count[6] <- length(setdiff(setdiff(female.18c.genvar, female.25c.genvar), female.18c[, 1]))

# make boxes
# ============================================================
par(mar = c(0, 4.5, 1.2, 1.5))

plot(c(0, 4), c(-max(c(length(female.18c.genvar), length(female.25c.genvar), length(male.18c.genvar), length(male.25c.genvar)))*1.2, 0), type = "n", xlab = "", ylab = "", axes = FALSE)
rect(0, 0 - sum(left.count[1:2]), 1, 0, col = brewer.pal(9, "Purples")[9], border = NA)
rect(0, 0 - sum(left.count[1:4]), 1, 0 - sum(left.count[1:2]), col = brewer.pal(9, "Purples")[3], border = NA)
rect(0, 0 - sum(left.count[1:5]), 1, 0 - sum(left.count[1:4]), col = brewer.pal(9, "Reds")[9], border = NA)
rect(0, 0 - sum(left.count[1:6]), 1, 0 - sum(left.count[1:5]), col = brewer.pal(9, "Reds")[3], border = NA)

rect(3, 0 - sum(right.count[1:2]), 4, 0, col = brewer.pal(9, "Purples")[9], border = NA)
rect(3, 0 - sum(right.count[1:4]), 4, 0 - sum(right.count[1:2]), col = brewer.pal(9, "Purples")[3], border = NA)
rect(3, 0 - sum(right.count[1:5]), 4, 0 - sum(right.count[1:4]), col = brewer.pal(9, "Blues")[9], border = NA)
rect(3, 0 - sum(right.count[1:6]), 4, 0 - sum(right.count[1:5]), col = brewer.pal(9, "Blues")[3], border = NA)


xspline(c(1, 3, 3, 1), c(-left.count[1], -right.count[1], 0, 0), shape = c(0, 0, 0, 0), open = FALSE, col = rgb(col2rgb(brewer.pal(9, "Purples")[9])[1, 1], col2rgb(brewer.pal(9, "Purples")[9])[2, 1], col2rgb(brewer.pal(9, "Purples")[9])[3, 1], 75, max = 255), border = NA)

xspline(c(1, 1.5, 2, 2.5, 3, 3, 2.5, 2, 1.5, 1), c(-sum(left.count[1:2]), -sum(left.count[1:2]), (-sum(left.count[1:2]) - sum(right.count[1:3]))/2, -sum(right.count[1:3]), -sum(right.count[1:3]), -sum(right.count[1:2]), -sum(right.count[1:2]), (-sum(left.count[1]) - sum(right.count[1:2]))/2, -left.count[1], -left.count[1]), shape = c(0, 1, 0, 1, 0), open = FALSE, col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 75, max = 255), border = NA)

xspline(c(1, 1.5, 2, 2.5, 3, 3, 2.5, 2, 1.5, 1), c(-sum(left.count[1:3]), -sum(left.count[1:3]), (-sum(left.count[1:3]) - sum(right.count[1:2]))/2, -sum(right.count[1:2]), -sum(right.count[1:2]), -sum(right.count[1]), -sum(right.count[1]), (-sum(left.count[1:2]) - sum(right.count[1]))/2, -sum(left.count[1:2]), -sum(left.count[1:2])), shape = c(0, 1, 0, 1, 0), open = FALSE, col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 75, max = 255), border = NA)

rect(1, -sum(left.count[1:4]), 3, -sum(left.count[1:3]), col = rgb(col2rgb(brewer.pal(9, "Purples")[3])[1, 1], col2rgb(brewer.pal(9, "Purples")[3])[2, 1], col2rgb(brewer.pal(9, "Purples")[3])[3, 1], 160, max = 255), border = NA)

text(2, 190, expression(paste("Genes with ", italic(H)^2, " >0")), cex = 7/par("ps")/par("cex"), pos = 3)
text(-1, -3000, "\u2640", family = "Arial", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(0.5, -150, expression(paste("25 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"), pos = 3)
text(3.5, -150, expression(paste("18 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"), pos = 3)

segments(c(-0.15, -0.15), c(0, -sum(left.count[1:4])), c(0, 0), c(0, -sum(left.count[1:4])), lwd = 0.5)
segments(c(-0.15), c(0), c(-0.15), c(-sum(left.count[1:4])), lwd = 0.5)

text(-1.05, -sum(left.count[1:4])/2 + 300, expression(paste(italic(H)^2, " > 0")), cex = 7/par("ps")/par("cex"))
text(-1.05, -sum(left.count[1:4])/2 - 300, paste("at both\nn = ", formatC(sum(left.count[1:4]), big.mark = ",")), cex = 7/par("ps")/par("cex"))
text(2, -sum(left.count[1])/2, formatC(left.count[1], big.mark = ","), cex = 7/par("ps")/par("cex"))

arrows(1.3, -sum(left.count[1:3]) - 200, 1, -sum(left.count[1:2]) - left.count[3]/2, length = 0.05)
text(1.3, -sum(left.count[1:3]) - 100, formatC(left.count[3], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

arrows(2.7, -sum(right.count[1:3]) - 200, 3, -sum(right.count[1:2]) - right.count[3]/2, length = 0.05)
text(2.7, -sum(right.count[1:3]) - 100, formatC(right.count[3], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

arrows(1.3, -sum(left.count[1:5]) - 200, 1, -sum(left.count[1:4]) - left.count[5]/2, length = 0.05)
text(1.3, -sum(left.count[1:5]) - 100, formatC(left.count[5], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

arrows(2.7, -sum(right.count[1:5]) - 200, 3, -sum(right.count[1:4]) - right.count[5]/2, length = 0.05)
text(2.7, -sum(right.count[1:5]) - 100, formatC(right.count[5], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

# legend
rect(c(-1, -1, -1, -1), -sum(left.count[1:6]) - 1800 - (1:4)*350, c(-0.5, -0.5, -0.5, -0.5), -sum(left.count[1:6]) - 1800 - (1:4)*350 + 160, border = NA, col = c(brewer.pal(9, "Purples")[9], brewer.pal(9, "Purples")[3], brewer.pal(9, "Reds")[9], brewer.pal(9, "Reds")[3]))

rect(c(1.5, 1.5, 1.5, 1.5), -sum(left.count[1:6]) - 1800 - (1:4)*350, c(2,2,2,2), -sum(left.count[1:6]) - 1800 - (1:4)*350 + 160, border = NA, col = c(brewer.pal(9, "Purples")[9], brewer.pal(9, "Purples")[3], brewer.pal(9, "Blues")[9], brewer.pal(9, "Blues")[3]))

text(0.5, -sum(left.count[1:6]) - 1800 - (1:4)*350 + 80, c("eQTL > 0", "eQTL = 0", "eQTL > 0", "eQTL = 0"), cex = 7/par("ps")/par("cex"))

text(-0.75, -sum(left.count[1:6]) - 1720, expression(paste("25 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"))
text(1.75, -sum(left.count[1:6]) - 1720, expression(paste("18 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# ribbons in males
# ============================================================

left.count <- numeric(6)
right.count <- numeric(6)

# 25c, 
left.count[1] <- length(intersect(male.25c[, 1], male.18c[, 1]))
left.count[2] <- length(setdiff(intersect(intersect(male.25c.genvar, male.18c.genvar), male.25c[, 1]), male.18c[, 1]))
left.count[3] <- length(setdiff(intersect(intersect(male.25c.genvar, male.18c.genvar), male.18c[, 1]), male.25c[, 1]))
left.count[4] <- length(setdiff(setdiff(intersect(male.25c.genvar, male.18c.genvar), male.18c[, 1]), male.25c[, 1]))
left.count[5] <- length(intersect(male.25c[, 1], setdiff(male.25c.genvar, male.18c.genvar)))
left.count[6] <- length(setdiff(setdiff(male.25c.genvar, male.18c.genvar), male.25c[, 1]))

# 18c
right.count[1] <- left.count[1]
right.count[2] <- left.count[3]
right.count[3] <- left.count[2]
right.count[4] <- left.count[4]
right.count[5] <- length(intersect(male.18c[, 1], setdiff(male.18c.genvar, male.25c.genvar)))
right.count[6] <- length(setdiff(setdiff(male.18c.genvar, male.25c.genvar), male.18c[, 1]))

# make boxes
# ============================================================

plot(c(0, 4), c(-max(c(length(female.18c.genvar), length(female.25c.genvar), length(male.18c.genvar), length(male.25c.genvar)))*1.2, 0), type = "n", xlab = "", ylab = "", axes = FALSE)
rect(0, 0 - sum(left.count[1:2]), 1, 0, col = brewer.pal(9, "Purples")[9], border = NA)
rect(0, 0 - sum(left.count[1:4]), 1, 0 - sum(left.count[1:2]), col = brewer.pal(9, "Purples")[3], border = NA)
rect(0, 0 - sum(left.count[1:5]), 1, 0 - sum(left.count[1:4]), col = brewer.pal(9, "Reds")[9], border = NA)
rect(0, 0 - sum(left.count[1:6]), 1, 0 - sum(left.count[1:5]), col = brewer.pal(9, "Reds")[3], border = NA)

rect(3, 0 - sum(right.count[1:2]), 4, 0, col = brewer.pal(9, "Purples")[9], border = NA)
rect(3, 0 - sum(right.count[1:4]), 4, 0 - sum(right.count[1:2]), col = brewer.pal(9, "Purples")[3], border = NA)
rect(3, 0 - sum(right.count[1:5]), 4, 0 - sum(right.count[1:4]), col = brewer.pal(9, "Blues")[9], border = NA)
rect(3, 0 - sum(right.count[1:6]), 4, 0 - sum(right.count[1:5]), col = brewer.pal(9, "Blues")[3], border = NA)


xspline(c(1, 3, 3, 1), c(-left.count[1], -right.count[1], 0, 0), shape = c(0, 0, 0, 0), open = FALSE, col = rgb(col2rgb(brewer.pal(9, "Purples")[9])[1, 1], col2rgb(brewer.pal(9, "Purples")[9])[2, 1], col2rgb(brewer.pal(9, "Purples")[9])[3, 1], 75, max = 255), border = NA)

xspline(c(1, 1.5, 2, 2.5, 3, 3, 2.5, 2, 1.5, 1), c(-sum(left.count[1:2]), -sum(left.count[1:2]), (-sum(left.count[1:2]) - sum(right.count[1:3]))/2, -sum(right.count[1:3]), -sum(right.count[1:3]), -sum(right.count[1:2]), -sum(right.count[1:2]), (-sum(left.count[1]) - sum(right.count[1:2]))/2, -left.count[1], -left.count[1]), shape = c(0, 1, 0, 1, 0), open = FALSE, col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 75, max = 255), border = NA)

xspline(c(1, 1.5, 2, 2.5, 3, 3, 2.5, 2, 1.5, 1), c(-sum(left.count[1:3]), -sum(left.count[1:3]), (-sum(left.count[1:3]) - sum(right.count[1:2]))/2, -sum(right.count[1:2]), -sum(right.count[1:2]), -sum(right.count[1]), -sum(right.count[1]), (-sum(left.count[1:2]) - sum(right.count[1]))/2, -sum(left.count[1:2]), -sum(left.count[1:2])), shape = c(0, 1, 0, 1, 0), open = FALSE, col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 75, max = 255), border = NA)

rect(1, -sum(left.count[1:4]), 3, -sum(left.count[1:3]), col = rgb(col2rgb(brewer.pal(9, "Purples")[3])[1, 1], col2rgb(brewer.pal(9, "Purples")[3])[2, 1], col2rgb(brewer.pal(9, "Purples")[3])[3, 1], 160, max = 255), border = NA)

text(2, 190, expression(paste("Genes with ", italic(H)^2, " >0")), cex = 7/par("ps")/par("cex"), pos = 3)
text(-1, -3000, "\u2642", family = "Arial", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0.5, -150, expression(paste("25 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"), pos = 3)
text(3.5, -150, expression(paste("18 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"), pos = 3)

segments(c(-0.15, -0.15), c(0, -sum(left.count[1:4])), c(0, 0), c(0, -sum(left.count[1:4])), lwd = 0.5)
segments(c(-0.15), c(0), c(-0.15), c(-sum(left.count[1:4])), lwd = 0.5)

text(-1.05, -sum(left.count[1:4])/2 + 300, expression(paste(italic(H)^2, " > 0")), cex = 7/par("ps")/par("cex"))
text(-1.05, -sum(left.count[1:4])/2 - 300, paste("at both\nn = ", formatC(sum(left.count[1:4]), big.mark = ",")), cex = 7/par("ps")/par("cex"))
text(2, -sum(left.count[1])/2, formatC(left.count[1], big.mark = ","), cex = 7/par("ps")/par("cex"))

arrows(1.3, -sum(left.count[1:3]) - 200, 1, -sum(left.count[1:2]) - left.count[3]/2, length = 0.05)
text(1.3, -sum(left.count[1:3]) - 100, formatC(left.count[3], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

arrows(2.7, -sum(right.count[1:3]) - 200, 3, -sum(right.count[1:2]) - right.count[3]/2, length = 0.05)
text(2.7, -sum(right.count[1:3]) - 100, formatC(right.count[3], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

arrows(1.3, -sum(left.count[1:5]) - 200, 1, -sum(left.count[1:4]) - left.count[5]/2, length = 0.05)
text(1.3, -sum(left.count[1:5]) - 100, formatC(left.count[5], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

arrows(2.7, -sum(right.count[1:5]) - 200, 3, -sum(right.count[1:4]) - right.count[5]/2, length = 0.05)
text(2.7, -sum(right.count[1:5]) - 100, formatC(right.count[5], big.mark = ","), pos = 1, cex = 7/par("ps")/par("cex"))

# legend
rect(c(-1, -1, -1, -1), -sum(left.count[1:6]) - 500 - (1:4)*350, c(-0.5, -0.5, -0.5, -0.5), -sum(left.count[1:6]) - 500 - (1:4)*350 + 160, border = NA, col = c(brewer.pal(9, "Purples")[9], brewer.pal(9, "Purples")[3], brewer.pal(9, "Reds")[9], brewer.pal(9, "Reds")[3]))

rect(c(1.5, 1.5, 1.5, 1.5), -sum(left.count[1:6]) - 500 - (1:4)*350, c(2,2,2,2), -sum(left.count[1:6]) - 500 - (1:4)*350 + 160, border = NA, col = c(brewer.pal(9, "Purples")[9], brewer.pal(9, "Purples")[3], brewer.pal(9, "Blues")[9], brewer.pal(9, "Blues")[3]))

text(0.5, -sum(left.count[1:6]) - 500 - (1:4)*350 + 80, c("eQTL > 0", "eQTL = 0", "eQTL > 0", "eQTL = 0"), cex = 7/par("ps")/par("cex"))

text(-0.75, -sum(left.count[1:6]) - 420, expression(paste("25 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"))
text(1.75, -sum(left.count[1:6]) - 420, expression(paste("18 ", {}*degree, "C")), cex = 7/par("ps")/par("cex"))
text(grconvertX(0.05 + file.width/4/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot effect size comparison
# ============================================================
par(mar = c(2.5, 3, 1, 0.3))

female.col <- rep(brewer.pal(9, "Purples")[9], nrow(female.eqtl.pair))
female.col[female.eqtl.pair[, 3] == "18c"] <- brewer.pal(9, "Blues")[9]
female.col[female.eqtl.pair[, 3] == "25c"] <- brewer.pal(9, "Reds")[9]

cat("range of effects (standardized):", range(as.numeric(female.eqtl.pair[, 4])/sqrt(as.numeric(female.eqtl.pair[, 6]))), range(as.numeric(female.eqtl.pair[, 7])/sqrt(as.numeric(female.eqtl.pair[, 9]))), "\n")

cat("there are a total of", nrow(female.eqtl.pair), "eQTL-gene pairs.\n")
print(table(female.eqtl.pair[,3]))

plot(as.numeric(female.eqtl.pair[, 7])/sqrt(as.numeric(female.eqtl.pair[, 9])), as.numeric(female.eqtl.pair[, 4])/sqrt(as.numeric(female.eqtl.pair[, 6])), col = female.col, cex = 0.2, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", axes = FALSE, lwd = 0.2)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = -0.5, cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))

title(xlab = expression(paste("Standardized effect at 25 ", {}*degree, "C")), mgp = c(0.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Standardized effect at 18 ", {}*degree, "C")), mgp = c(1.35, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

legend("topleft", col = c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9], brewer.pal(9, "Purples")[9]), pch = 1, legend = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), "Both"), y.intersp = 0.6, x.intersp = 0.5, pt.lwd = 1, cex = 6/par("ps")/par("cex"), bty = "n")

text(1, -0.5, "\u2640", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 1, family = "Arial")
text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# male
# ============================================================

male.col <- rep(brewer.pal(9, "Purples")[9], nrow(male.eqtl.pair))
male.col[male.eqtl.pair[, 3] == "18c"] <- brewer.pal(9, "Blues")[9]
male.col[male.eqtl.pair[, 3] == "25c"] <- brewer.pal(9, "Reds")[9]

cat("range of effects (standardized):", range(as.numeric(male.eqtl.pair[, 4])/sqrt(as.numeric(male.eqtl.pair[, 6]))), range(as.numeric(male.eqtl.pair[, 7])/sqrt(as.numeric(male.eqtl.pair[, 9]))), "\n")

cat("there are a total of", nrow(male.eqtl.pair), "eQTL-gene pairs.\n")
print(table(male.eqtl.pair[,3]))

plot(as.numeric(male.eqtl.pair[, 7])/sqrt(as.numeric(male.eqtl.pair[, 9])), as.numeric(male.eqtl.pair[, 4])/sqrt(as.numeric(male.eqtl.pair[, 6])), col = male.col, cex = 0.2, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", axes = FALSE, lwd = 0.2)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = -0.5, cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))

title(xlab = expression(paste("Standardized effect at 25 ", {}*degree, "C")), mgp = c(0.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Standardized effect at 18 ", {}*degree, "C")), mgp = c(1.35, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

legend("topleft", col = c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9], brewer.pal(9, "Purples")[9]), pch = 1, legend = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), "Both"), y.intersp = 0.6, x.intersp = 0.5, pt.lwd = 1, cex = 6/par("ps")/par("cex"), bty = "n")

text(1, -0.5, "\u2642", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 1, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# female prediction
# ============================================================

rownames(female.eqtl.pred.25c) <- female.eqtl.pred.25c[, 1]
female.eqtl.gene <- na.omit(female.eqtl.pred.25c)[, 1]
female.eqtl.col <- rep("grey80", length(female.eqtl.gene))
female.eqtl.col[female.eqtl.gene %in% female.eqtl.pair[female.eqtl.pair[, 3] == "both", 1]] <- brewer.pal(9, "Reds")[9]
female.eqtl.gene <- female.eqtl.gene[order(female.eqtl.col, decreasing = T)]
female.eqtl.col <- female.eqtl.col[order(female.eqtl.col, decreasing = T)]

plot(female.gxe[female.eqtl.gene, 6]/(female.gxe[female.eqtl.gene, 6] + female.gxe[female.eqtl.gene, 7]), as.numeric(female.eqtl.pred.25c[female.eqtl.gene, 4]), col = female.eqtl.col, pch = 16, cex = 0.4, axes = FALSE, xlim = c(0, 1), ylim = c(-0.1, 1))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(xlab = "GxE", mgp = c(0.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Prediction correlation", mgp = c(1.35, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

legend("topright", col = c("grey80", brewer.pal(9, "Reds")[9]), pch = 16, legend = c("No shared eQTLs", "Shared eQTLs"), y.intersp = 0.6, x.intersp = 0.5, pt.lwd = 1, cex = 6/par("ps")/par("cex"), bty = "n")

text(0.9, 0.3, "\u2640", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 1, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# male prediction
# ============================================================

rownames(male.eqtl.pred.25c) <- male.eqtl.pred.25c[, 1]
male.eqtl.gene <- na.omit(male.eqtl.pred.25c)[, 1]
male.eqtl.col <- rep("grey80", length(male.eqtl.gene))
male.eqtl.col[male.eqtl.gene %in% male.eqtl.pair[male.eqtl.pair[, 3] == "both", 1]] <- brewer.pal(9, "Blues")[9]
male.eqtl.gene <- male.eqtl.gene[order(male.eqtl.col, decreasing = T)]
male.eqtl.col <- male.eqtl.col[order(male.eqtl.col, decreasing = T)]

plot(male.gxe[male.eqtl.gene, 6]/(male.gxe[male.eqtl.gene, 6] + male.gxe[male.eqtl.gene, 7]), as.numeric(male.eqtl.pred.25c[male.eqtl.gene, 4]), col = male.eqtl.col, pch = 16, cex = 0.4, axes = FALSE, xlim = c(0, 1), ylim = c(-0.1, 1))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(xlab = "GxE", mgp = c(0.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Prediction correlation", mgp = c(1.35, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

legend("topright", col = c("grey80", brewer.pal(9, "Blues")[9]), pch = 16, legend = c("No shared eQTLs", "Shared eQTLs"), y.intersp = 0.6, x.intersp = 0.5, pt.lwd = 1, cex = 6/par("ps")/par("cex"), bty = "n")

text(0.9, 0.3, "\u2642", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 1, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("f")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# remove hot spot analysis TF
# ============================================================

par(mar = c(0, 1, 0.5, 0.5))

tf.name <- gsub(".*_TFBS_", "", tf.list)

plot(c(-10, length(tf.name) + 1), c(0, 14), type = "n", axes = FALSE, xlab = "", ylab = "")

image(1:length(tf.name), 1:8, cbind(male.tf.fold.enrich[, 3], female.tf.fold.enrich[, 3], NA, male.tf.fold.enrich[, 2], female.tf.fold.enrich[, 2], NA, male.tf.fold.enrich[, 1], female.tf.fold.enrich[, 1]), add = T, breaks = c(seq(0, 1, 0.2), 2:6), col = c(brewer.pal(9, "Blues")[c(8, 6, 4, 2, 1)], brewer.pal(9, "Reds")[c(1, 2, 4, 6, 8)]))

segments(0.5, seq(0.5, 8.5, 1), length(tf.name) + 0.5, seq(0.5, 8.5, 1))
segments(seq(0.5, length(tf.name) + 0.5, 1), 0.5, seq(0.5, length(tf.name) + 0.5, 1), 2.5)
segments(seq(0.5, length(tf.name) + 0.5, 1), 3.5, seq(0.5, length(tf.name) + 0.5, 1), 5.5)
segments(seq(0.5, length(tf.name) + 0.5, 1), 6.5, seq(0.5, length(tf.name) + 0.5, 1), 8.5)

for (i in 1:length(tf.name)) {
	
	this.tf <- tf.name[i]
	text(i - 0.6, 9, parse(text = paste("italic(\"", this.tf, "\")", sep = "")), srt = 60, pos = 4, cex = 7/par("ps")/par("cex"))
	#text(i - 0.6, 9, parse(text = paste("expression(italic(", this.tf, "))", sep = "")), pos = 4)
	
}

p.adj <- matrix(p.adjust(c(female.tf.fold.enrich.p, male.tf.fold.enrich.p), method = "BH"), ncol = 6, byrow = FALSE)
female.tf.fold.enrich.p <- p.adj[, 1:3]
male.tf.fold.enrich.p <- p.adj[, 4:6]


text(which(male.tf.fold.enrich.p[, 3] < 0.05), 1, "**", cex = 7/par("ps")/par("cex"))
text(which(male.tf.fold.enrich.p[, 3] >= 0.05 & male.tf.fold.enrich.p[, 3] < 0.1), 1, "*", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 3] < 0.05), 2, "**", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 3] >= 0.05 & female.tf.fold.enrich.p[, 3] < 0.1), 2, "*", cex = 7/par("ps")/par("cex"))

text(which(male.tf.fold.enrich.p[, 2] < 0.05), 4, "**", cex = 7/par("ps")/par("cex"))
text(which(male.tf.fold.enrich.p[, 2] >= 0.05 & male.tf.fold.enrich.p[, 3] < 0.1), 4, "*", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 2] < 0.05), 5, "**", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 2] >= 0.05 & female.tf.fold.enrich.p[, 3] < 0.1), 5, "*", cex = 7/par("ps")/par("cex"))

text(which(male.tf.fold.enrich.p[, 1] < 0.05), 7, "**", cex = 7/par("ps")/par("cex"))
text(which(male.tf.fold.enrich.p[, 1] >= 0.05 & male.tf.fold.enrich.p[, 3] < 0.1), 7, "*", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 1] < 0.05), 8, "**", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 1] >= 0.05 & female.tf.fold.enrich.p[, 3] < 0.1), 8, "*", cex = 7/par("ps")/par("cex"))

text(-12.5, 7.5, expression(paste("25 ", {}*degree, "C only vs random")), cex = 7/par("ps")/par("cex"), pos = 4)
text(-12.5, 4.5, expression(paste("Common vs 25 ", {}*degree, "C only")), cex = 7/par("ps")/par("cex"), pos = 4)
text(-12.5, 1.5, expression(paste("18 ", {}*degree, "C only vs", " 25 ", {}*degree, "C only")), cex = 7/par("ps")/par("cex"), pos = 4)


text(0, 1, "\u2642", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0, 2, "\u2640", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(0, 4, "\u2642", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0, 5, "\u2640", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(0, 7, "\u2642", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0, 8, "\u2640", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("g")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

image(-10:-1, c(12, 13), matrix(c(seq(0.1, 0.9, 0.2), seq(1.5, 5.5, 1), rep(-0.1, 10)), ncol = 2), add = T, breaks = c(seq(-0.2, 1, 0.2), 2:6), col = c("white", brewer.pal(9, "Blues")[c(8, 6, 4, 2, 1)], brewer.pal(9, "Reds")[c(1, 2, 4, 6, 8)]))
segments(-10.5, 11.5, -0.5, 11.5)
segments(-10.5, 12.5, -0.5, 12.5)
segments(seq(-10.5, -0.5, 1), 11.5, seq(-10.5, -0.5, 1), 12.5)
text(seq(-10.5, -0.5, 1) + 1, 11, c(seq(0, 1, 0.2), 2:6), srt = 60, pos = 2)
text(-5.5, 13.5, "Fold enrichment", cex = 7/par("ps")/par("cex"))

text(7, 13.5, "* FDR = 0.10", cex = 7/par("ps")/par("cex"), pos = 2)
text(7, 12.25, "** FDR = 0.05", cex = 7/par("ps")/par("cex"), pos = 2)


dev.off()