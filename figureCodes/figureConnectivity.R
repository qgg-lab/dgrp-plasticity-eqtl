# =====================================================
# = connectivity versus GEI and stabilizing selection =
# =====================================================

args <- commandArgs(TRUE) # args <- c("reportData/female.wgcna.RData", "reportData/male.wgcna.RData",  "reportData/sigma2.RData", "reportData/female.adj.qgGxE.RData", "reportData/male.adj.qgGxE.RData")
library("RColorBrewer")

# load data
# ============================================================

load(args[1])
female.blup.25c.cor <- blup.25c.cor
female.blup.18c.cor <- blup.18c.cor
load(args[2])
male.blup.25c.cor <- blup.25c.cor
male.blup.18c.cor <- blup.18c.cor
load(args[3])
# ratio is Vm/Vg
male.sigma.ratio <- log10(sigma2.all.male[sigma2.all.male[, 9] > male.cut, 6]/sigma2.all.male[sigma2.all.male[, 9] > male.cut, 1])
female.sigma.ratio <- log10(sigma2.all.female[sigma2.all.female[, 9] > female.cut, 6]/sigma2.all.female[sigma2.all.female[, 9] > female.cut, 1])
male.sigma.ratio[male.sigma.ratio > 4] <- 4
male.sigma.ratio[male.sigma.ratio < -4] <- -4
female.sigma.ratio[female.sigma.ratio > 4] <- 4
female.sigma.ratio[female.sigma.ratio < -4] <- -4


load(args[4])
female.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)
female.gene.name <- gene.name
load(args[5])
male.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)
male.gene.name <- gene.name

female.gei <- female.gxe[, 6]/(female.gxe[, 6] + female.gxe[, 7])
female.int.fdr <- p.adjust(as.numeric(female.gxe[, 10]), method = "BH")
female.gei.gene <- female.gene.name[female.int.fdr < 0.05]

male.gei <- male.gxe[, 6]/(male.gxe[, 6] + male.gxe[, 7])
male.int.fdr <- p.adjust(as.numeric(male.gxe[, 10]), method = "BH")
male.gei.gene <- male.gene.name[which(male.int.fdr < 0.05)]

# set up plot
# ============================================================

file.width = 89
cairo_pdf(file = args[6], width = file.width/25.4, height = file.width/25.4*0.5, family = args[7])
par(mfrow = c(1, 2), las = 1, tcl = -0.2, mar = c(2, 2.5, 1, 0.1), ps = 7, lwd = 0.5)

# compute connectivity at 18C and 25C
# ============================================================

female.gene.connect.25c <- apply(abs(female.blup.25c.cor), 2, function(x) { return( (sum(x) - 1)/(length(x) - 1) ) })

male.gene.connect.25c <- apply(abs(male.blup.25c.cor), 2, function(x) { return( (sum(x) - 1)/(length(x) - 1) ) })

# plot connectivity with stabilizing selection
# ============================================================

female.common.gene <- intersect(names(female.gene.connect.25c), names(female.sigma.ratio))
male.common.gene <- intersect(names(male.gene.connect.25c), names(male.sigma.ratio))

cat("female big point genes:", sum(abs(female.sigma.ratio[female.common.gene]) >= 3.999), "/", length(female.common.gene), "\n")
cat("male big point genes:", sum(abs(male.sigma.ratio[male.common.gene]) >= 3.999), "/", length(male.common.gene), "\n")


plot(c(0.05, 0.22), c(-4, 4), axes = FALSE, type = "n", xlab = "", ylab = "")

# points(female.gene.connect.25c[female.common.gene], female.sigma.ratio[female.common.gene], pch = 1, cex = 0.05, col = "grey50", lwd = 0.2)
# lines(loess.smooth(female.gene.connect.25c[female.common.gene], female.sigma.ratio[female.common.gene], span = 0.1), col = brewer.pal(9, "Reds")[9], lwd = 1)

points(female.gene.connect.25c[setdiff(female.common.gene, female.gei.gene)], female.sigma.ratio[setdiff(female.common.gene, female.gei.gene)], pch = 1, cex = 0.05, col = "grey50", lwd = 0.2)

points(female.gene.connect.25c[intersect(female.common.gene, female.gei.gene)], female.sigma.ratio[intersect(female.common.gene, female.gei.gene)], pch = 1, cex = 0.05, col = brewer.pal(9, "Reds")[5], lwd = 0.2)

cor.test(female.gene.connect.25c[setdiff(female.common.gene, female.gei.gene)], female.sigma.ratio[setdiff(female.common.gene, female.gei.gene)], method = "spearman")
cor.test(female.gene.connect.25c[intersect(female.common.gene, female.gei.gene)], female.sigma.ratio[intersect(female.common.gene, female.gei.gene)], method = "spearman")

lines(loess.smooth(female.gene.connect.25c[setdiff(female.common.gene, female.gei.gene)], female.sigma.ratio[setdiff(female.common.gene, female.gei.gene)], span = 0.1), col = "grey50", lwd = 1)
lines(loess.smooth(female.gene.connect.25c[intersect(female.common.gene, female.gei.gene)], female.sigma.ratio[intersect(female.common.gene, female.gei.gene)], span = 0.1), col = brewer.pal(9, "Reds")[5], lwd = 1)
lines(loess.smooth(female.gene.connect.25c[female.common.gene], female.sigma.ratio[female.common.gene], span = 0.1), col = brewer.pal(9, "Reds")[9], lwd = 1)

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, at = seq(-4, 4, 2), label = c(expression(paste("<", 10^-4)), expression(10^-2), 1, expression(10^2), expression(paste(">", 10^4))), mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = "Gene connectivity", mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Stabilizing selection (", italic(V[m]), "/", italic(V[g]), ")")), mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.23, -2, "\u2640", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 2, family = "Arial")
female.cor <- formatC(cor.test(female.gene.connect.25c[female.common.gene], female.sigma.ratio[female.common.gene], method = "spearman")$estimate, format = "f", digits = 2)
female.cor.p <- formatC(cor.test(female.gene.connect.25c[female.common.gene], female.sigma.ratio[female.common.gene], method = "spearman")$p.value, format = "e", digits = 2)

# text(2, 0.2, expression(paste(italic(sigma[e]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
# text(1, 2.9, expression(paste(italic(sigma[e]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

# legend (hard code)
segments(c(0.20, 0.20, 0.20), c(0.5, 1.5, 2.5), c(0.22, 0.22, 0.22), c(0.5, 1.5, 2.5), col = c(brewer.pal(9, "Reds")[9], brewer.pal(9, "Reds")[5], "grey50"))
points(c(0.21, 0.21), c(1.5, 2.5), pch = 1, cex = 0.5, col = c(brewer.pal(9, "Reds")[5], "grey50", brewer.pal(9, "Reds")[5]))
text(c(0.21, 0.21, 0.21), c(0.5, 1.5, 2.5), c("All", "GxE", "No GxE"), pos = 2, cex = 6/par("ps")/par("cex"))

text(0.05, 3.5, parse(text = paste("paste(italic(r), \" = \", ", female.cor, ", \" (\", italic(P), \" = \", ", female.cor.p, ", \")\"", ")")), pos = 4)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


plot(c(0.05, 0.22), c(-4, 4), axes = FALSE, type = "n", xlab = "", ylab = "")

# points(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], pch = 1, cex = 0.05, col = "grey50", lwd = 0.2)
# lines(loess.smooth(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], span = 0.1), col = brewer.pal(9, "Reds")[9], lwd = 1)

points(male.gene.connect.25c[setdiff(male.common.gene, male.gei.gene)], male.sigma.ratio[setdiff(male.common.gene, male.gei.gene)], pch = 1, cex = 0.05, col = "grey50", lwd = 0.2)

points(male.gene.connect.25c[intersect(male.common.gene, male.gei.gene)], male.sigma.ratio[intersect(male.common.gene, male.gei.gene)], pch = 1, cex = 0.05, col = brewer.pal(9, "Blues")[5], lwd = 0.2)

cor.test(male.gene.connect.25c[setdiff(male.common.gene, male.gei.gene)], male.sigma.ratio[setdiff(male.common.gene, male.gei.gene)], method = "spearman")
cor.test(male.gene.connect.25c[intersect(male.common.gene, male.gei.gene)], male.sigma.ratio[intersect(male.common.gene, male.gei.gene)], method = "spearman")

lines(loess.smooth(male.gene.connect.25c[setdiff(male.common.gene, male.gei.gene)], male.sigma.ratio[setdiff(male.common.gene, male.gei.gene)], span = 0.1), col = "grey20", lwd = 1)
lines(loess.smooth(male.gene.connect.25c[intersect(male.common.gene, male.gei.gene)], male.sigma.ratio[intersect(male.common.gene, male.gei.gene)], span = 0.1), col = brewer.pal(9, "Blues")[5], lwd = 1)
lines(loess.smooth(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], span = 0.1), col = brewer.pal(9, "Blues")[9], lwd = 1)


axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, at = seq(-4, 4, 2), label = c(expression(paste("<", 10^-4)), expression(10^-2), 1, expression(10^2), expression(paste(">", 10^4))), mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = "Gene connectivity", mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Stabilizing selection (", italic(V[m]), "/", italic(V[g]), ")")), mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.23, -2, "\u2642", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 2, family = "Arial")
male.cor <- formatC(cor.test(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], method = "spearman")$estimate, format = "f", digits = 2)
male.cor.p <- formatC(cor.test(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], method = "spearman")$p.value, format = "e", digits = 2)

# text(2, 0.2, expression(paste(italic(sigma[e]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
# text(1, 2.9, expression(paste(italic(sigma[e]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
# legend (hard code)
segments(c(0.20, 0.20, 0.20), c(0.5, 1.5, 2.5), c(0.22, 0.22, 0.22), c(0.5, 1.5, 2.5), col = c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Blues")[5], "grey50"))
points(c(0.21, 0.21), c(1.5, 2.5), pch = 1, cex = 0.5, col = c(brewer.pal(9, "Blues")[5], "grey50", brewer.pal(9, "Blues")[5]))
text(c(0.21, 0.21, 0.21), c(0.5, 1.5, 2.5), c("All", "GxE", "No GxE"), pos = 2, cex = 6/par("ps")/par("cex"))

text(0.05, 3.5, parse(text = paste("paste(italic(r), \" = \", ", male.cor, ", \" (\", italic(P), \" = \", ", male.cor.p, ", \")\"", ")")), pos = 4)

text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
