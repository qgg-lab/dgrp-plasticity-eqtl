# ===========================================
# = plot example of GSEA for decanalization =
# ===========================================

args <- commandArgs(TRUE) # args <- c("reportData/female.decan.gsea.example.RData", "reportData/male.decan.gsea.example.RData")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 120
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width*0.3/25.4, family = args[4])
par(mfrow = c(1, 3), las = 1, tcl = -0.2, mar = c(2, 2.2, 1, 0.5), ps = 7, lwd = 0.5, xpd = TRUE)

# function for plot
plot.gsea <- function(p, S, score, name) {
  
  plot(c(1, length(p)), c(-1.1, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  
  # make gradient
  # ============================================================
  
  bar.height = (-1.1 - par('usr')[3])/2
  
  image(x = 1:length(p), y = c(-1.1 - bar.height*1.5, -1.1 - bar.height*0.5), z = matrix(rep(score[, 2], 2), ncol = 2), col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(21), breaks = seq(-1.05, 1.05, 0.1), add = TRUE)
  rect(0.5, -1.1 - bar.height*2, length(p) + 0.5, -1.1)
  
  # draw gene set positions
  # ============================================================
  
  segments(S, -1.1, S, -1.07)
  
  # draw ES score profile
  # ============================================================
  
  lines(1:length(p), cumsum(p))
  
  # axis
  # ============================================================
  
  axis(side = 2, at = seq(-1, 1, 0.5), lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
  axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0, 3000, 1000), cex.axis = 7/par("ps")/par("cex"))
  title(xlab = "Gene rank", mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
  title(ylab = "Enrichment score", mgp = c(1.4, 0, 0), cex.lab = 7/par("ps")/par("cex"))
  box(bty = "l")
  
  # text
  # ============================================================
  
  text(length(p)/2, 0.8, label = name, pos = 3, cex = 7/par("ps")/par("cex"))
  
}


# # load female data
# # ============================================================
#
# load(args[1])
#
# # plot GO:0003954 NADH dehydrogenase activity
# # ============================================================
#
# plot.gsea(gsea.example[[1]]$p, gsea.example[[1]]$S, gsea.example[[1]]$score, "Cell morphogenesis")
# text(0, -0.7, "\u2640", pos = 4, cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
# text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot GO:0016787	hydrolase activity
# ============================================================

# plot.gsea(gsea.example[[2]]$p, gsea.example[[2]]$S, gsea.example[[2]]$score, "Iron ion binding")
# text(0, -0.7, "\u2640", pos = 4, cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
# text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# # plot GO:0050909  sensory perception of taste
#
# plot.gsea(gsea.example[[2]]$p, gsea.example[[2]]$S, gsea.example[[2]]$score, "Sensory perception of taste")
# text(0, -0.7, "\u2640", pos = 4, cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
# text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# load male data
# ============================================================

load(args[1])

# plot GO:0005214 structural constituent of chitin-based cuticle
# ============================================================

plot.gsea(gsea.example[[1]]$p, gsea.example[[1]]$S, gsea.example[[1]]$score, "structural constituent of\nchitin-based cuticle")
text(0, -0.7, "\u2640", pos = 4, cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], family = "Arial")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# males
load(args[2])


# plot GO:0005549 odorant binding
# ============================================================

plot.gsea(gsea.example[[1]]$p, gsea.example[[1]]$S, gsea.example[[1]]$score, "odorant binding")
text(0, -0.7, "\u2642", pos = 4, cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], family = "Arial")

text(grconvertX(0.05 + file.width/3/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot GO:0003700 DNA-binding transcription factor activity
# ============================================================

plot.gsea(gsea.example[[2]]$p, gsea.example[[2]]$S, gsea.example[[2]]$score, "DNA-binding\ntranscription factor activity")
text(0, -0.7, "\u2642", pos = 4, cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], family = "Arial")

text(grconvertX(0.05 + file.width*2/3/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()

