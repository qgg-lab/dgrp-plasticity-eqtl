# =========================================================
# = make figure to summarize canalization, decanalization =
# =========================================================

# 1. plot comparison of sigma2_e, similar
# 2. plot comparison of sigma2_g, different
#    emphasize the contrast between 1 and 2
# ============================================================

args <- commandArgs(TRUE) # args = c("reportData/female.adj.qgSingleTemp.RData", "reportData/male.adj.qgSingleTemp.RData", "reportData/female.adj.qgVarHet.RData", "reportData/male.adj.qgVarHet.RData", "reportData/gene.info")
library("RColorBrewer")

# load data
# ============================================================

gene.info <- read.table(args[5], header = FALSE, as.is = TRUE, row.names = 1)

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[6], width = file.width/25.4, height = 1.5*file.width/25.4, family = args[8])
layout(matrix(c(1, 1, 2, 2, 3,3, 4, 4, 5, 6, 6, 7, 8, 9, 9, 10), byrow = T, ncol = 4, nrow = 4), widths = c(0.3, 0.1, 0.1, 0.3), heights = c(0.3, 0.3, 0.2, 0.1))

### with wolbachia as covariate
# ============================================================

load(args[1])
load(args[3])
# count significant genes in female
cat("female 18c significant genes no.:", sum(p.adjust(single.temp[, 9], method = "BH") < 0.05), "\n")
cat("female 25c significant genes no.:", sum(p.adjust(single.temp[, 12], method = "BH") < 0.05), "\n")
cat("shared by both:", sum(p.adjust(single.temp[, 9], method = "BH") < 0.05 & p.adjust(single.temp[, 12], method = "BH") < 0.05), "\n")

cat("female 18c decanalized genes no.:", sum(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & (as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2 & p.adjust(single.temp[, 9], method = "BH") < 0.05), na.rm = T), "\n")
cat("female 18c canalized genes no.:", sum(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & (as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2 & p.adjust(single.temp[, 12], method = "BH") < 0.05), na.rm = T), "\n")

# plot female comparison of sigma2_e
# ============================================================

par(las = 1, tcl = -0.2, mar = c(2, 1.5, 1, 0.5), ps = 7, lwd = 0.5, xpd = F)

plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 10])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 9])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 10])), sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 9])), col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 10])), sqrt(as.numeric(var.het[, 9])), span = 0.001), col = brewer.pal(9, "Reds")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[e]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[e]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[e]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1, 2.9, expression(paste(italic(sigma[e]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot female comparison of sigma2_G
# ============================================================


plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 7])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 8])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 7])), sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 8])), col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 7])), sqrt(as.numeric(var.het[, 8])), span = 0.001), col = brewer.pal(9, "Reds")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[G]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[G]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[G]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1, 2.9, expression(paste(italic(sigma[G]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# male
# ============================================================
# ============================================================


load(args[2])
load(args[4])
# count significant genes in female
cat("male 18c significant genes no.:", sum(p.adjust(single.temp[, 9], method = "BH") < 0.05), "\n")
cat("male 25c significant genes no.:", sum(p.adjust(single.temp[, 12], method = "BH") < 0.05), "\n")
cat("shared by both:", sum(p.adjust(single.temp[, 9], method = "BH") < 0.05 & p.adjust(single.temp[, 12], method = "BH") < 0.05), "\n")

cat("male 18c decanalized genes no.:", sum(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & (as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2 & p.adjust(single.temp[, 9], method = "BH") < 0.05), na.rm = T), "\n")
cat("male 18c canalized genes no.:", sum(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & (as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2 & p.adjust(single.temp[, 12], method = "BH") < 0.05), na.rm = T), "\n")


# plot male comparison of sigma2_e
# ============================================================

plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 10])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 9])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 10])), sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 11]), method = "BH") < 0.05 & ((as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) > 2) | (as.numeric(var.het[, 10])/as.numeric(var.het[, 9]) < 0.5))), 9])), col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 10])), sqrt(as.numeric(var.het[, 9])), span = 0.001), col = brewer.pal(9, "Blues")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[e]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[e]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[e]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(1, 2.9, expression(paste(italic(sigma[e]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot male comparison of sigma2_G
# ============================================================


plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 7])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 8])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 7])), sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 12]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 7])/as.numeric(var.het[, 8]) > 2) & (p.adjust(single.temp[, 12], method = "BH") < 0.05)) | ((as.numeric(var.het[, 8])/as.numeric(var.het[, 7]) > 2) & (p.adjust(single.temp[, 9], method = "BH") < 0.05)) )), 8])), col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 7])), sqrt(as.numeric(var.het[, 8])), span = 0.001), col = brewer.pal(9, "Blues")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[G]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[G]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[G]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(1, 2.9, expression(paste(italic(sigma[G]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)






# plot small variance
# ============================================================

par(las = 1, tcl = -0.2, mar = c(1, 2, 1, 0.5), ps = 7, lwd = 0.5, xpd = T)

plot(c(-2, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")
lines(seq(-2, 2, 0.01), dnorm(seq(-2, 2, 0.01), mean = 0, sd = 0.2), lwd = 0.5)
text(0, 2.4, "Environment A", cex = 7/par("ps")/par("cex"), pos = 1)
box(bty = "l")
title(xlab = "Genetic value", mgp = c(0.05, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(0.4, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


plot(c(-2, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")

arrows(-2.5, 1.5, 1.5, 1.5, angle = 30, length = 0.06)
text(-0.5, 1.5, "Genetic", pos = 3, cex = 7/par("ps")/par("cex"))
text(-0.5, 1.5, "decanalization", pos = 1, cex = 7/par("ps")/par("cex"))


arrows(-2.5, 0, 1.5, 0, angle = 30, length = 0.06, code = 1)
text(-0.5, 0, "Genetic", pos = 3, cex = 7/par("ps")/par("cex"))
text(-0.5, 0, "canalization", pos = 1, cex = 7/par("ps")/par("cex"))


plot(c(-2, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")
lines(seq(-2, 2, 0.01), dnorm(seq(-2, 2, 0.01), mean = 0, sd = 0.6), lwd = 0.5)
text(0, 2.4, "Environment B", cex = 7/par("ps")/par("cex"), pos = 1)
box(bty = "l")
title(xlab = "Genetic value", mgp = c(0.05, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(0.4, 0, 0), cex.lab = 7/par("ps")/par("cex"))



par(las = 1, tcl = -0.2, mar = c(1, 2, 0.2, 0.5), ps = 7, lwd = 0.5, xpd = T)
plot(c(-2, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")
box()
segments(c(-0.2, 0, -0.1), c(1.6, 1, 0.4), c(-0.2, 0, -0.1) + 0.3, c(1.6, 1, 0.4))
points(c(-0.2, 0, -0.1), c(1.6, 1, 0.4), cex = 0.8)
points(c(-0.2, 0, -0.1) + 0.3, c(1.6, 1, 0.4), cex = 0.8)
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = c(0.4, 1, 1.6), label = c("c", "b", "a"), cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Phenotypic value", mgp = c(0.05, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Genotype", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

plot(c(-2, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")


plot(c(-2, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")
box()
segments(c(-1, 0.8, -0.2), c(1.6, 1, 0.4), c(-1, 0.8, -0.2) + 0.3, c(1.6, 1, 0.4))
points(c(-1, 0.8, -0.2), c(1.6, 1, 0.4), cex = 0.8)
points(c(-1, 0.8, -0.2) + 0.3, c(1.6, 1, 0.4), cex = 0.8)
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = c(0.4, 1, 1.6), label = c("c", "b", "a"), cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Phenotypic value", mgp = c(0.05, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Genotype", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))






dev.off()


# without wolbachia as covariate
# ============================================================

cairo_pdf(file = args[7], width = file.width/25.4, height = file.width/25.4, family = args[8])
par(mfrow = c(2, 2), las = 1, tcl = -0.2, mar = c(2, 1.5, 1, 0.5), ps = 7, lwd = 0.5)



# load data
# ============================================================

load(args[1])
load(args[3])
# count significant genes in female
cat("female 18c significant genes no.:", sum(p.adjust(single.temp[, 3], method = "BH") < 0.05), "\n")
cat("female 25c significant genes no.:", sum(p.adjust(single.temp[, 6], method = "BH") < 0.05), "\n")
cat("shared by both:", sum(p.adjust(single.temp[, 3], method = "BH") < 0.05 & p.adjust(single.temp[, 6], method = "BH") < 0.05), "\n")

cat("female 18c decanalized genes no.:", sum(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & (as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2 & p.adjust(single.temp[, 3], method = "BH") < 0.05), na.rm = T), "\n")
cat("female 18c canalized genes no.:", sum(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & (as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2 & p.adjust(single.temp[, 6], method = "BH") < 0.05), na.rm = T), "\n")


# plot female comparison of sigma2_e
# ============================================================

plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5))), 4])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5))), 3])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5)), 4])), sqrt(as.numeric(var.het[p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5)), 3])), col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 4])), sqrt(as.numeric(var.het[, 3])), span = 0.001), col = brewer.pal(9, "Reds")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[e]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[e]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[e]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1, 2.9, expression(paste(italic(sigma[e]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# plot female comparison of sigma2_G
# ============================================================

plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 1])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 2])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 1])), sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 2])), col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 1])), sqrt(as.numeric(var.het[, 2])), span = 0.001), col = brewer.pal(9, "Reds")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[G]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[G]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[G]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1, 2.9, expression(paste(italic(sigma[G]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# male
# ============================================================
# ============================================================


load(args[2])
load(args[4])
# count significant genes in female
cat("male 18c significant genes no.:", sum(p.adjust(single.temp[, 3], method = "BH") < 0.05), "\n")
cat("male 25c significant genes no.:", sum(p.adjust(single.temp[, 6], method = "BH") < 0.05), "\n")
cat("shared by both:", sum(p.adjust(single.temp[, 3], method = "BH") < 0.05 & p.adjust(single.temp[, 6], method = "BH") < 0.05), "\n")

cat("male 18c decanalized genes no.:", sum(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & (as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2 & p.adjust(single.temp[, 3], method = "BH") < 0.05), na.rm = T), "\n")
cat("male 18c canalized genes no.:", sum(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & (as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2 & p.adjust(single.temp[, 6], method = "BH") < 0.05), na.rm = T), "\n")


# plot male comparison of sigma2_e
# ============================================================

plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5))), 4])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5))), 3])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5)), 4])), sqrt(as.numeric(var.het[p.adjust(as.numeric(var.het[, 5]), method = "BH") < 0.05 & ((as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) > 2) | (as.numeric(var.het[, 4])/as.numeric(var.het[, 3]) < 0.5)), 3])), col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 4])), sqrt(as.numeric(var.het[, 3])), span = 0.001), col = brewer.pal(9, "Blues")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[e]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[e]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[e]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(1, 2.9, expression(paste(italic(sigma[e]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# plot male comparison of sigma2_G
# ============================================================

plot(c(0, 3), c(0, 3), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 1])), sqrt(as.numeric(var.het[!(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 2])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 1])), sqrt(as.numeric(var.het[(p.adjust(as.numeric(var.het[, 6]), method = "BH") < 0.05 & ( ((as.numeric(var.het[, 1])/as.numeric(var.het[, 2]) > 2) & (p.adjust(single.temp[, 6], method = "BH") < 0.05)) | ((as.numeric(var.het[, 2])/as.numeric(var.het[, 1]) > 2) & (p.adjust(single.temp[, 3], method = "BH") < 0.05)) )), 2])), col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(loess.smooth(sqrt(as.numeric(var.het[, 1])), sqrt(as.numeric(var.het[, 2])), span = 0.001), col = brewer.pal(9, "Blues")[9])

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[G]), " at 25 ", {}*degree, "C")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[G]), " at 18 ", {}*degree, "C")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.2, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[G]^2), " (25 ", {}*degree, "C/18 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(1, 2.9, expression(paste(italic(sigma[G]^2), " (18 ", {}*degree, "C/25 ", {}*degree, "C)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()


