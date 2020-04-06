# =====================================================
# = connectivity versus GEI and stabilizing selection =
# =====================================================

args <- commandArgs(TRUE) # args <- c("reportData/female.wgcna.RData", "reportData/male.wgcna.RData",  "reportData/sigma2.RData")
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

# set up plot
# ============================================================

file.width = 89
cairo_pdf(file = args[4], width = file.width/25.4, height = file.width/25.4*0.5, family = args[5])
par(mfrow = c(1, 2), las = 1, tcl = -0.2, mar = c(2, 2.5, 1, 0.1), ps = 7, lwd = 0.5)

# compute connectivity at 18C and 25C
# ============================================================

female.gene.connect.25c <- apply(abs(female.blup.25c.cor), 2, function(x) { return( (sum(x) - 1)/(length(x) - 1) ) })

male.gene.connect.25c <- apply(abs(male.blup.25c.cor), 2, function(x) { return( (sum(x) - 1)/(length(x) - 1) ) })

# plot connectivity with stabilizing selection
# ============================================================

female.common.gene <- intersect(names(female.gene.connect.25c), names(female.sigma.ratio))
male.common.gene <- intersect(names(male.gene.connect.25c), names(male.sigma.ratio))

plot(c(0.05, 0.22), c(-4, 4), axes = FALSE, type = "n", xlab = "", ylab = "")

points(female.gene.connect.25c[female.common.gene], female.sigma.ratio[female.common.gene], pch = 1, cex = 0.05, col = "grey50", lwd = 0.2)
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

text(0.05, 3, parse(text = paste("paste(italic(r), \" = \", ", female.cor, ", \" (\", italic(P), \" = \", ", female.cor.p, ", \")\"", ")")), pos = 4)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


plot(c(0.05, 0.22), c(-4, 4), axes = FALSE, type = "n", xlab = "", ylab = "")

points(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], pch = 1, cex = 0.05, col = "grey50", lwd = 0.2)
lines(loess.smooth(male.gene.connect.25c[male.common.gene], male.sigma.ratio[male.common.gene], span = 0.1), col = brewer.pal(9, "Reds")[9], lwd = 1)

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

text(0.05, 3, parse(text = paste("paste(italic(r), \" = \", ", male.cor, ", \" (\", italic(P), \" = \", ", male.cor.p, ", \")\"", ")")), pos = 4)

text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
