# ==============================================================================
# = make figure to compare expression in females and males between 18 and 25 C =
# ==============================================================================

args <- commandArgs(TRUE)

data.file <- args[1]
symbol.file <- args[2]
library("RColorBrewer")

# load data
# ============================================================

gene.fpkm <- read.table(data.file, header = TRUE, as.is = TRUE)
gene.symbol <- read.table(symbol.file, header = FALSE, as.is = TRUE)

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width*1.1/25.4, family = args[4])
par(las = 1, tcl = -0.2, mar = c(2, 1.8, 1, 0.9), ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE), height = c(1, 1.2))

# plot female scatter plot
# ============================================================

plot(c(0, 15), c(0, 15), axes = FALSE, type = "n")

female.18c <- log2(gene.fpkm$f18c_0 + 1)
female.25c <- log2(gene.fpkm$f25c_0 + 1)

points(female.18c[abs(female.18c - female.25c) < 1 & !grepl("XLOC", gene.fpkm[, 1])], female.25c[abs(female.18c - female.25c) < 1 & !grepl("XLOC", gene.fpkm[, 1])], cex = 0.1, pch = 16, col = "grey50")
points(female.18c[abs(female.18c - female.25c) < 1 & grepl("XLOC", gene.fpkm[, 1])], female.25c[abs(female.18c - female.25c) < 1 & grepl("XLOC", gene.fpkm[, 1])], cex = 0.1, pch = 16, col = "grey20")

points(female.18c[female.18c - female.25c >= 1 & !grepl("XLOC", gene.fpkm[, 1])], female.25c[female.18c - female.25c >= 1 & !grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlGnBu")[3])
points(female.18c[female.18c - female.25c >= 1 & grepl("XLOC", gene.fpkm[, 1])], female.25c[female.18c - female.25c >= 1 & grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlGnBu")[9])

points(female.18c[female.18c - female.25c <= -1 & !grepl("XLOC", gene.fpkm[, 1])], female.25c[female.18c - female.25c <= -1 & !grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlOrRd")[3])
points(female.18c[female.18c - female.25c <= -1 & grepl("XLOC", gene.fpkm[, 1])], female.25c[female.18c - female.25c <= -1 & grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlOrRd")[9])

axis(side = 1, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0, 0))
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
box(bty = "l", lwd = 0.5)

title(xlab = expression(paste("18 ", {}*degree, "C log", {}[2], "(FPKM + 1)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("25 ", {}*degree, "C log", {}[2], "(FPKM + 1)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(12, 2, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")

# plot male scatter plot
# ============================================================

plot(c(0, 15), c(0, 15), axes = FALSE, type = "n")

male.18c <- log2(gene.fpkm$m18c_0 + 1)
male.25c <- log2(gene.fpkm$m25c_0 + 1)

points(male.18c[abs(male.18c - male.25c) < 1 & !grepl("XLOC", gene.fpkm[, 1])], male.25c[abs(male.18c - male.25c) < 1 & !grepl("XLOC", gene.fpkm[, 1])], cex = 0.1, pch = 16, col = "grey50")
points(male.18c[abs(male.18c - male.25c) < 1 & grepl("XLOC", gene.fpkm[, 1])], male.25c[abs(male.18c - male.25c) < 1 & grepl("XLOC", gene.fpkm[, 1])], cex = 0.1, pch = 16, col = "grey20")

points(male.18c[male.18c - male.25c >= 1 & !grepl("XLOC", gene.fpkm[, 1])], male.25c[male.18c - male.25c >= 1 & !grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlGnBu")[3])
points(male.18c[male.18c - male.25c >= 1 & grepl("XLOC", gene.fpkm[, 1])], male.25c[male.18c - male.25c >= 1 & grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlGnBu")[9])

points(male.18c[male.18c - male.25c <= -1 & !grepl("XLOC", gene.fpkm[, 1])], male.25c[male.18c - male.25c <= -1 & !grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlOrRd")[3])
points(male.18c[male.18c - male.25c <= -1 & grepl("XLOC", gene.fpkm[, 1])], male.25c[male.18c - male.25c <= -1 & grepl("XLOC", gene.fpkm[, 1])], cex = 0.3, pch = 16, col = brewer.pal(9, "YlOrRd")[9])

axis(side = 1, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0, 0))
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
box(bty = "l", lwd = 0.5)

title(xlab = expression(paste("18 ", {}*degree, "C log", {}[2], "(FPKM + 1)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("25 ", {}*degree, "C log", {}[2], "(FPKM + 1)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(12, 2, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")

# plot Table
# ============================================================

par(mar = c(0, 1, 0, 1), ps = 8, lwd = 0.5, xpd = TRUE)
plot(c(0, 1), c(0, 0.7), type = "n", axes = FALSE)

segments(0, 0.7, 1, 0.7, lwd = 1)
text(0.8, 0.66, "\u2640", cex = 12/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 2, family = "Arial")
text(1, 0.66, "\u2642", cex = 12/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 2, family = "Arial")
segments(0, 0.62, 1, 0.62, lwd = 0.5)


text(0, 0.575, "Total (FPKM > 0 in either temperature)", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.575, formatC(sum(female.18c > 0 | female.25c > 0), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.575, formatC(sum(male.18c > 0 | male.25c > 0), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(0, 0.525, "     Annotated", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.525, formatC(sum((female.18c > 0 | female.25c > 0) & !grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.525, formatC(sum((male.18c > 0 | male.25c > 0) & !grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(0, 0.475, "     NTR", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.475, formatC(sum((female.18c > 0 | female.25c) > 0 & grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.475, formatC(sum((male.18c > 0 | male.25c > 0) & grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))

text(0, 0.425, expression(paste({} < {}, "2 fold change")), pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.425, formatC(sum((abs(female.18c - female.25c) < 1 & ( female.18c > 0 | female.25c > 0 ))), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.425, formatC(sum((abs(male.18c - male.25c) < 1 & ( male.18c > 0 | male.25c > 0 ))), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
points(0.07, 0.375, pch = 16, cex = 1.1, col = "grey50"); text(0, 0.375, "        Annotated", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.375, formatC(sum((abs(female.18c - female.25c) < 1 & ( female.18c > 0 | female.25c > 0 ) & !grepl("XLOC", gene.fpkm[, 1]))), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.375, formatC(sum((abs(male.18c - male.25c) < 1 & ( male.18c > 0 | male.25c > 0 ) & !grepl("XLOC", gene.fpkm[, 1]))), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
points(0.07, 0.325, pch = 16, cex = 1.1, col = "grey20"); text(0, 0.325, "        NTR", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.325, formatC(sum((abs(female.18c - female.25c) < 1 & ( female.18c > 0 | female.25c > 0 ) & grepl("XLOC", gene.fpkm[, 1]))), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.325, formatC(sum((abs(male.18c - male.25c) < 1 & ( male.18c > 0 | male.25c > 0 ) & grepl("XLOC", gene.fpkm[, 1]))), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))

text(0, 0.275, expression(paste("18 ", {}*degree, "C", "/", "25 ", {}*degree, "C", {}>={}, "2 fold")), pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.275, formatC(sum(female.18c - female.25c >= 1), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.275, formatC(sum(male.18c - male.25c >= 1), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
points(0.07, 0.225, pch = 16, cex = 1.1, col = brewer.pal(9, "YlGnBu")[3]); text(0, 0.225, "        Annotated", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.225, formatC(sum(female.18c - female.25c >= 1 & !grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.225, formatC(sum(male.18c - male.25c >= 1 & !grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
points(0.07, 0.175, pch = 16, cex = 1.1, col = brewer.pal(9, "YlGnBu")[9]); text(0, 0.175, "        NTR", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.175, formatC(sum(female.18c - female.25c >= 1 & grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.175, formatC(sum(male.18c - male.25c >= 1 & grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))

text(0, 0.125, expression(paste("25 ", {}*degree, "C", "/", "18 ", {}*degree, "C", {}>={}, "2 fold")), pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.125, formatC(sum(female.18c - female.25c <= -1), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.125, formatC(sum(male.18c - male.25c <= -1), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
points(0.07, 0.075, pch = 16, cex = 1.1, col = brewer.pal(9, "YlOrRd")[3]); text(0, 0.075, "        Annotated", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.075, formatC(sum(female.18c - female.25c <= -1 & !grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.075, formatC(sum(male.18c - male.25c <= -1 & !grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
points(0.07, 0.025, pch = 16, cex = 1.1, col = brewer.pal(9, "YlOrRd")[9]); text(0, 0.025, "        NTR", pos = 4, cex = 7/par("ps")/par("cex"))
text(0.8, 0.025, formatC(sum(female.18c - female.25c <= -1 & grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))
text(1, 0.025, formatC(sum(male.18c - male.25c <= -1 & grepl("XLOC", gene.fpkm[, 1])), big.mark = ","), pos = 2, cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()

