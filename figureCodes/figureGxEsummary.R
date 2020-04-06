# ========================================
# = make figure to summarize GxE results =
# ========================================

args <- commandArgs(TRUE) # args <- c("reportData/female.adj.qgGxE.RData", "reportData/male.adj.qgGxE.RData", "reportData/female.qgGxEexample.RData", 22, "CG6293", 25, "CG14309")

library("RColorBrewer")
library("plotrix")

# load data
# ============================================================

load(args[1])
female.gxe <- gxe
load(args[2])
male.gxe <- gxe
load(args[3])
col.pal <- brewer.pal(8, "Dark2")

exp1.index <- as.numeric(args[4])
exp1.gene <- args[5]
exp2.index <- as.numeric(args[6])
exp2.gene <- args[7]

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[8], width = file.width/25.4, height = file.width/25.4, family = args[9])
par(las = 1, tcl = -0.2, mar = c(1.6, 2.2, 1, 0.5), mfrow = c(2, 2), xpd = TRUE, ps = 7, lwd = 0.5)

# plot female 1st expression
# ============================================================

# blup
plot(c(1, 2), range(c(blup.25c[exp1.index, ], blup.18c[exp1.index, ])), xlim = c(0.8, 2.2), type = "n", xlab = "", ylab = "", axes = FALSE)

for (i in 1:ncol(blup.25c)) {
  
  points(c(1, 2), c(blup.25c[exp1.index, i], blup.18c[exp1.index, i]), type = "b", lwd = 0.5, col = col.pal[i %% 8 + 1], cex = 0.5, pch = 16)
  
}

axis(side = 1, at = c(1, 2), labels = c(expression(paste("25 ", {}*degree, "C")), expression(paste("18 ", {}*degree, "C"))), lwd = 0.5, mgp = c(0.8, 0, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = parse(text = paste("paste(italic(", exp1.gene, "),", "\" expression (BLUP log\"[2], \" scale)\")")), mgp = c(1, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

text(1.5, grconvertY(0.9, from = "nfc", to = "user"), paste("GEI = ", formatC(gxe[exp1.index, 1]/(gxe[exp1.index, 1] + gxe[exp1.index, 2]), format = "f", digits = 2), sep = ""), cex = 7/par("ps")/par("cex"))
# text(1, 5.4, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 1)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot female CG14309 expression
# ============================================================

plot(c(1, 2), range(c(blup.25c[exp2.index, ], blup.18c[exp2.index, ])), xlim = c(0.8, 2.2), type = "n", xlab = "", ylab = "", axes = FALSE)

for (i in 1:ncol(blup.25c)) {
  
  points(c(1, 2), c(blup.25c[exp2.index, i], blup.18c[exp2.index, i]), type = "b", lwd = 0.5, col = col.pal[i %% 8 + 1], cex = 0.5, pch = 16)
  
}

axis(side = 1, at = c(1, 2), labels = c(expression(paste("25 ", {}*degree, "C")), expression(paste("18 ", {}*degree, "C"))), lwd = 0.5, mgp = c(0.8, 0, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = parse(text = paste("paste(italic(", exp2.gene, "),", "\" expression (BLUP log\"[2], \" scale)\")")), mgp = c(1, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

#axis.break(axis = 2, breakpos = 11.09, style = "zigzag")

text(1.5, grconvertY(0.9, from = "nfc", to = "user"), paste("GEI = ", formatC(gxe[exp2.index, 1]/(gxe[exp2.index, 1] + gxe[exp2.index, 2]), format = "f", digits = 2), sep = ""), cex = 7/par("ps")/par("cex"))
# text(2, 10.7, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 1)
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# prepare data to plot
# ============================================================

female.gei <- female.gxe[, 6]/(female.gxe[, 6] + female.gxe[, 7])
female.int.fdr <- p.adjust(as.numeric(female.gxe[, 10]), method = "BH")
female.line.fdr <- p.adjust(as.numeric(female.gxe[, 9]), method = "BH")

cat("there are ", sum(female.int.fdr < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 in females.\n")
cat("there are ", sum(female.line.fdr < 0.05, na.rm = TRUE), " genes with line FDR < 0.05 in females.\n")
cat("there are ", sum(female.int.fdr < 0.05 | female.line.fdr < 0.05, na.rm = TRUE), " genes with int | line FDR < 0.05 in females.\n")


male.gei <- male.gxe[, 6]/(male.gxe[, 6] + male.gxe[, 7])
male.int.fdr <- p.adjust(as.numeric(male.gxe[, 10]), method = "BH")
male.line.fdr <- p.adjust(as.numeric(male.gxe[, 9]), method = "BH")

cat("there are ", sum(male.int.fdr < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 in males.\n")
cat("there are ", sum(male.line.fdr < 0.05, na.rm = TRUE), " genes with line FDR < 0.05 in males.\n")
cat("there are ", sum(male.int.fdr < 0.05 | male.line.fdr < 0.05, na.rm = TRUE), " genes with int | line FDR < 0.05 in males.\n")

# female histogram
# ============================================================

female.hist.plot <- hist(female.gei[female.int.fdr < 0.05 | female.line.fdr < 0.05], breaks = seq(0, 1, 0.05), plot = FALSE)
female.hist.plot$counts[1] <- ifelse(female.hist.plot$counts[1] > 1000, female.hist.plot$counts[1] - 800, female.hist.plot$counts[1])
plot(female.hist.plot, col = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 1000))
hist(female.gei[female.int.fdr < 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Reds")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "GEI", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.3, 0, 0), cex.lab = 7/par("ps")/par("cex"))
legend("topright", pch = 22, pt.bg = c(brewer.pal(9, "Reds")[9]),
       legend = c("FDR (GEI > 0) = 0.05"), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), seq(0, 1000, 200), labels = c(seq(0, 600, 200), 1600, 1800), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
axis.break(axis = 2, breakpos = 700, style = "zigzag")
text(0.7, 400, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# male histogram
# ============================================================

male.hist.plot <- hist(male.gei[male.int.fdr < 0.05 | male.line.fdr < 0.05], breaks = seq(0, 1, 0.05), plot = FALSE)
male.hist.plot$counts[1] <- ifelse(male.hist.plot$counts[1] > 1000, male.hist.plot$counts[1] - 1200, male.hist.plot$counts[1])
plot(male.hist.plot, col = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 1000))
hist(male.gei[male.int.fdr < 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Blues")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "GEI", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.3, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("topright", pch = 22, pt.bg = c(brewer.pal(9, "Blues")[9]),
       legend = c("FDR (GEI > 0) = 0.05"), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), seq(0, 1000, 200), labels = c(seq(0, 600, 200), 1800, 2000), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
axis.break(axis = 2, breakpos = 700, style = "zigzag")

text(0.7, 400, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()
