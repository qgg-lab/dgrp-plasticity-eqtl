# ======================
# = make network plots =
# ======================

args <- commandArgs(TRUE) # args <- c("reportData/female.wgcna.RData", "reportData/male.wgcna.RData", "reportData/female.simCorr18C.RData", "reportData/male.simCorr18C.RData")
source("reportCode/alluvial.R")
library("RColorBrewer")
library("WGCNA")

# load data
# ============================================================

load(args[1])

# set up plot
file.width = 89
cairo_pdf(file = args[5], width = file.width/25.4, height = file.width*1.2/25.4, family = args[6])
par(las = 1, tcl = -0.2, mar = c(0.1, 0.1, 0.1, 0.1), ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 1, 2, 2, 3, 3, 4, 5), ncol = 2, byrow = TRUE), height = c(1, 0.18, 1, 1))

module.col <- c("grey80", brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"))

# function to re-order and get mean connectivity across module
# order based on module connectivity then connectivity within module
# ============================================================

orderCor <- function(cor.mat, module) {
	
	module <- as.character(module)
	unique.module <- unique(module)
	module.mean.corr <- numeric(length(unique.module))
	names(module.mean.corr) <- unique.module
	gene.conn <- numeric(nrow(cor.mat))
	names(gene.conn) <- rownames(cor.mat)
	
	for (i in unique.module) {
		
		this.cor <- cor.mat[which(module == i), which(module == i)]
		module.mean.corr[i] <- mean(abs(this.cor[lower.tri(this.cor, diag = FALSE)]))
		gene.conn[rownames(this.cor)] <- rowMeans(this.cor)
		
	}
	
	gene.order <- order(module.mean.corr[module], gene.conn)
	cor.mat <- cor.mat[gene.order, gene.order]
	module <- module[gene.order]
	
	# get mean cross module correlation
	for (i in 1:(length(unique.module) - 1)) {
		for (j in (i+1):length(unique.module)) {
			
			module1 <- which(module == unique.module[i])
			module2 <- which(module == unique.module[j])
			mean.cross.cor <- mean(cor.mat[module1, module2])
			cor.mat[module1, module2] <- mean.cross.cor
			cor.mat[module2, module1] <- mean.cross.cor
			
		}
	}
  	
	return(list(cor = cor.mat, mod = module, gene.name = colnames(cor.mat)))
	
}


# order based on module connectivity then connectivity within module
# ============================================================

cor.25c.order <- orderCor(blup.25c.cor, blup.25c.tree)
cor.18c.order <- orderCor(blup.18c.cor, blup.18c.tree)

plot(c(0, nrow(cor.25c.order[[1]])*3), c(0, nrow(cor.25c.order[[1]])*1.1), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
image(1:nrow(cor.25c.order[[1]]), 1:ncol(cor.25c.order[[1]]), cor.25c.order[[1]], col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(101), useRaster = T, add = TRUE)

left.unique.module <- unique(cor.25c.order[[2]])
left.module.start <- sapply(split(1:length(cor.25c.order[[2]]), factor(cor.25c.order[[2]], levels = left.unique.module)), min)
left.module.end <- sapply(split(1:length(cor.25c.order[[2]]), factor(cor.25c.order[[2]], levels = left.unique.module)), max)

for (i in 1:length(left.module.start)) {
	
	rect(left.module.start[i] - 0.5, left.module.start[i] - 0.5, left.module.end[i] + 0.5, left.module.end[i] + 0.5, lwd = 0.2)
	
}

rect(0.5, 0.5, nrow(cor.25c.order[[1]]) + 0.5, nrow(cor.25c.order[[1]]) + 0.5, lwd = 0.2)

# draw module boxes
# ============================================================

rect(nrow(cor.25c.order[[1]])*1.05, left.module.start - 0.5, nrow(cor.25c.order[[1]])*1.1, left.module.end + 0.5, col = module.col[1:length(left.unique.module)], border = NA)

# 18C
# ============================================================
#
image((2*nrow(cor.25c.order[[1]])):(3*nrow(cor.25c.order[[1]]) - 1), 1:ncol(cor.25c.order[[1]]), cor.18c.order[[1]][nrow(cor.25c.order[[1]]):1, ], col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(101), useRaster = T, add = TRUE)

right.unique.module <- unique(cor.18c.order[[2]])
right.module.start <- sapply(split(1:length(cor.18c.order[[2]]), factor(cor.18c.order[[2]], levels = right.unique.module)), min)
right.module.end <- sapply(split(1:length(cor.18c.order[[2]]), factor(cor.18c.order[[2]], levels = right.unique.module)), max)

for (i in 1:length(right.module.start)) {
	
	rect(3*nrow(cor.25c.order[[1]]) - right.module.end[i] + 0.5, right.module.start[i] - 0.5, 3*nrow(cor.25c.order[[1]]) - right.module.start[i] + 1.5, right.module.end[i] + 0.5, lwd = 0.2)
	
}

rect(2*nrow(cor.25c.order[[1]]) - 0.5, 0.5, 3*nrow(cor.25c.order[[1]]) + 0.5, nrow(cor.25c.order[[1]]) + 0.5, lwd = 0.2)

rect(nrow(cor.25c.order[[1]])*1.9, right.module.start - 0.5, nrow(cor.25c.order[[1]])*1.95, right.module.end + 0.5, col = module.col[1:length(right.unique.module)], border = NA)

# add alluvial ribbons
# ============================================================

left.module.sum <- rep(0, length(left.module.start)); names(left.module.sum) <- names(left.module.start)
right.module.sum <- rep(0, length(right.module.start)); names(right.module.sum) <- names(right.module.start)

for (i in 1:length(left.module.sum)) {
	for (j in 1:length(right.module.sum)) {

		this.count <- length(intersect(
      cor.25c.order[[3]][cor.25c.order[[2]] == left.unique.module[i]],
      cor.18c.order[[3]][cor.18c.order[[2]] == right.unique.module[j]]))
		left.start <- left.module.start[i] + left.module.sum[i] - 0.5
		left.end <- left.module.start[i] + left.module.sum[i] - 0.5 + this.count

		right.start <- right.module.start[j] + right.module.sum[j] - 0.5
		right.end <- right.module.start[j] + right.module.sum[j] - 0.5 + this.count

		xspline(c(nrow(cor.25c.order[[1]])*1.1, nrow(cor.25c.order[[1]])*1.3, nrow(cor.25c.order[[1]])*1.5, nrow(cor.25c.order[[1]])*1.7, nrow(cor.25c.order[[1]])*1.9, nrow(cor.25c.order[[1]])*1.9, nrow(cor.25c.order[[1]])*1.7, nrow(cor.25c.order[[1]])*1.5, nrow(cor.25c.order[[1]])*1.3, nrow(cor.25c.order[[1]])*1.1, nrow(cor.25c.order[[1]])*1.1), c(left.start, left.start, (left.start + right.start)/2, right.start, right.start, right.end, right.end, (left.end + right.end)/2, left.end, left.end, left.start), shape = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0), open = FALSE, col = rgb(col2rgb(module.col[i])[1, 1], col2rgb(module.col[i])[2, 1], col2rgb(module.col[i])[3, 1], 75, max = 255), border = NA)

		left.module.sum[i] = left.module.sum[i] + this.count
		right.module.sum[j] = right.module.sum[j] + this.count

	}
}

# add text
# ============================================================

text(nrow(cor.25c.order[[1]])*1.5, 300, "Not in modules", cex = 7/par("ps")/par("cex"))
text(nrow(cor.25c.order[[1]])*0.5, nrow(cor.25c.order[[1]])*0.98, expression(paste("Modules at 25 ", {}*degree, "C")), pos = 3, cex = 7/par("ps")/par("cex"))
text(nrow(cor.25c.order[[1]])*1.5, nrow(cor.25c.order[[1]])*0.95, "\u2640", cex = 12/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(nrow(cor.25c.order[[1]])*2.5, nrow(cor.25c.order[[1]])*0.98, expression(paste("Modules at 18 ", {}*degree, "C")), pos = 3, cex = 7/par("ps")/par("cex"))
text(nrow(cor.25c.order[[1]])*(-0.01), nrow(cor.25c.order[[1]])*0.9, paste("n[genes] = ", formatC(nrow(cor.25c.order[[1]]), big.mark = ","), sep = ""), cex = 7/par("ps")/par("cex"), pos = 4)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# plot correlation colors
# ============================================================
cor100 <- seq(-1, 1, by = 0.02)
n.col <- length(cor100)
par(las = 1, tcl = -0.2, mar = c(1, 10, 0, 10), ps = 7, lwd = 0.5)

# plot the scale bar
plot(c(0, length(cor100) + 1), c(0, 3), type = "n", axes = FALSE, xlab = "", ylab = "")
image(1:length(cor100), 1:2, matrix(rep(cor100, 2), ncol = 2, byrow = F), col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(101), add = TRUE, useRaster = T)
rect(0.5, 0.5, length(cor100) + 0.5, 2.5, lwd = 0.5)
text(0, 1.5, "Pearson correlation", pos = 2, cex = 7/par("ps")/par("cex"), xpd = TRUE)
# add ticks
segments((1:length(cor100))[round(cor100, 2) %in% round(seq(-1, 1, 0.2), 2)], 0.5, (1:length(cor100))[round(cor100, 2) %in% round(seq(-1, 1, 0.2), 2)], 0.8, lwd = 0.5)
text((1:length(cor100))[round(cor100, 2) %in% round(seq(-1, 1, 0.2), 2)] + 10, -0.05, seq(-1, 1, 0.2), xpd = TRUE, srt = 60, pos = 2, cex = 6/par("ps")/par("cex"))


##### for males
load(args[2])
par(las = 1, tcl = -0.2, mar = c(0.1, 0.1, 0.1, 0.1), ps = 7, lwd = 0.5)

# order based on module connectivity then connectivity within module
# ============================================================

cor.25c.order <- orderCor(blup.25c.cor, blup.25c.tree)
cor.18c.order <- orderCor(blup.18c.cor, blup.18c.tree)

plot(c(0, nrow(cor.25c.order[[1]])*3), c(0, nrow(cor.25c.order[[1]])*1.1), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
image(1:nrow(cor.25c.order[[1]]), 1:ncol(cor.25c.order[[1]]), cor.25c.order[[1]], col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(101), useRaster = T, add = TRUE)

left.unique.module <- unique(cor.25c.order[[2]])
left.module.start <- sapply(split(1:length(cor.25c.order[[2]]), factor(cor.25c.order[[2]], levels = left.unique.module)), min)
left.module.end <- sapply(split(1:length(cor.25c.order[[2]]), factor(cor.25c.order[[2]], levels = left.unique.module)), max)

for (i in 1:length(left.module.start)) {
	
	rect(left.module.start[i] - 0.5, left.module.start[i] - 0.5, left.module.end[i] + 0.5, left.module.end[i] + 0.5, lwd = 0.2)
	
}

rect(0.5, 0.5, nrow(cor.25c.order[[1]]) + 0.5, nrow(cor.25c.order[[1]]) + 0.5, lwd = 0.2)

# draw module boxes
# ============================================================

rect(nrow(cor.25c.order[[1]])*1.05, left.module.start - 0.5, nrow(cor.25c.order[[1]])*1.1, left.module.end + 0.5, col = module.col[1:length(left.unique.module)], border = NA)

# 18C
# ============================================================
#
image((2*nrow(cor.25c.order[[1]])):(3*nrow(cor.25c.order[[1]]) - 1), 1:ncol(cor.25c.order[[1]]), cor.18c.order[[1]][nrow(cor.25c.order[[1]]):1, ], col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(101), useRaster = T, add = TRUE)

right.unique.module <- unique(cor.18c.order[[2]])
right.module.start <- sapply(split(1:length(cor.18c.order[[2]]), factor(cor.18c.order[[2]], levels = right.unique.module)), min)
right.module.end <- sapply(split(1:length(cor.18c.order[[2]]), factor(cor.18c.order[[2]], levels = right.unique.module)), max)

for (i in 1:length(right.module.start)) {
	
	rect(3*nrow(cor.25c.order[[1]]) - right.module.end[i] + 0.5, right.module.start[i] - 0.5, 3*nrow(cor.25c.order[[1]]) - right.module.start[i] + 1.5, right.module.end[i] + 0.5, lwd = 0.2)
	
}

rect(2*nrow(cor.25c.order[[1]]) - 0.5, 0.5, 3*nrow(cor.25c.order[[1]]) + 0.5, nrow(cor.25c.order[[1]]) + 0.5, lwd = 0.2)

rect(nrow(cor.25c.order[[1]])*1.9, right.module.start - 0.5, nrow(cor.25c.order[[1]])*1.95, right.module.end + 0.5, col = module.col[1:length(right.unique.module)], border = NA)

# add alluvial ribbons
# ============================================================

left.module.sum <- rep(0, length(left.module.start)); names(left.module.sum) <- names(left.module.start)
right.module.sum <- rep(0, length(right.module.start)); names(right.module.sum) <- names(right.module.start)

for (i in 1:length(left.module.sum)) {
	for (j in 1:length(right.module.sum)) {

		this.count <- length(intersect(
      cor.25c.order[[3]][cor.25c.order[[2]] == left.unique.module[i]],
      cor.18c.order[[3]][cor.18c.order[[2]] == right.unique.module[j]]))
		left.start <- left.module.start[i] + left.module.sum[i] - 0.5
		left.end <- left.module.start[i] + left.module.sum[i] - 0.5 + this.count

		right.start <- right.module.start[j] + right.module.sum[j] - 0.5
		right.end <- right.module.start[j] + right.module.sum[j] - 0.5 + this.count

		xspline(c(nrow(cor.25c.order[[1]])*1.1, nrow(cor.25c.order[[1]])*1.3, nrow(cor.25c.order[[1]])*1.5, nrow(cor.25c.order[[1]])*1.7, nrow(cor.25c.order[[1]])*1.9, nrow(cor.25c.order[[1]])*1.9, nrow(cor.25c.order[[1]])*1.7, nrow(cor.25c.order[[1]])*1.5, nrow(cor.25c.order[[1]])*1.3, nrow(cor.25c.order[[1]])*1.1, nrow(cor.25c.order[[1]])*1.1), c(left.start, left.start, (left.start + right.start)/2, right.start, right.start, right.end, right.end, (left.end + right.end)/2, left.end, left.end, left.start), shape = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0), open = FALSE, col = rgb(col2rgb(module.col[i])[1, 1], col2rgb(module.col[i])[2, 1], col2rgb(module.col[i])[3, 1], 75, max = 255), border = NA)

		left.module.sum[i] = left.module.sum[i] + this.count
		right.module.sum[j] = right.module.sum[j] + this.count

	}
}

# add text
# ============================================================

text(nrow(cor.25c.order[[1]])*1.5, 300, "Not in modules", cex = 7/par("ps")/par("cex"))
text(nrow(cor.25c.order[[1]])*0.5, nrow(cor.25c.order[[1]])*0.98, expression(paste("Modules at 25 ", {}*degree, "C")), pos = 3, cex = 7/par("ps")/par("cex"))
text(nrow(cor.25c.order[[1]])*1.5, nrow(cor.25c.order[[1]])*0.95, "\u2642", cex = 12/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(nrow(cor.25c.order[[1]])*2.5, nrow(cor.25c.order[[1]])*0.98, expression(paste("Modules at 18 ", {}*degree, "C")), pos = 3, cex = 7/par("ps")/par("cex"))
text(nrow(cor.25c.order[[1]])*(-0.01), nrow(cor.25c.order[[1]])*0.9, paste("n[genes] = ", formatC(nrow(cor.25c.order[[1]]), big.mark = ","), sep = ""), cex = 7/par("ps")/par("cex"), pos = 4)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# female permutation test
# ============================================================

load(args[3])
par(mar = c(2, 2, 1, 0.5))

plot(c(-1, 1), c(0, 4), type = "n", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 4))
sim.sd <- numeric(1000)

for (i in 1:1000) {
  lines(sim.corr.diff.gxe.dens[[i]], col = "grey50", lwd = 0.5)
	sim.sd[i] <- sim.corr.diff.gxe.dens[[i]]$sd
}

cat(sum(sim.sd < corr.diff.gxe.dens$sd), "\n")

lines(corr.diff.gxe.dens, col = brewer.pal(9, "Reds")[9], lwd = 1)
text(0.4, 4, "Observed", col = brewer.pal(9, "Reds")[9], cex = 7/par("ps")/par("cex"))
text(-0.6, 2, "Simulated", col = "grey50", cex = 7/par("ps")/par("cex"))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))

title(xlab = expression(paste(Delta, italic(r), " (18 ", {}*degree, "C - 25 ", {}*degree, "C)")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

text(0.9, 1.5, "\u2640", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 1, family = "Arial")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# male
# ============================================================

load(args[4])
par(mar = c(2, 2, 1, 0.5))

plot(c(-1, 1), c(0, 4), type = "n", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 4))
sim.sd <- numeric(1000)

for (i in 1:1000) {
  lines(sim.corr.diff.gxe.dens[[i]], col = "grey50", lwd = 0.5)
	sim.sd[i] <- sim.corr.diff.gxe.dens[[i]]$sd
}

cat(sum(sim.sd < corr.diff.gxe.dens$sd), "\n")

lines(corr.diff.gxe.dens, col = brewer.pal(9, "Blues")[9], lwd = 1)
text(0.4, 4, "Observed", col = brewer.pal(9, "Blues")[9], cex = 7/par("ps")/par("cex"))
text(-0.6, 2, "Simulated", col = "grey50", cex = 7/par("ps")/par("cex"))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))

title(xlab = expression(paste(Delta, italic(r), " (18 ", {}*degree, "C - 25 ", {}*degree, "C)")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

text(0.9, 1.5, "\u2642", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 1, family = "Arial")

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()
