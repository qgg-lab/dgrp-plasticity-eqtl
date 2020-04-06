# =======================================================
# = make figure of quantiles of standardized expression =
# =======================================================

args <- commandArgs(TRUE) # args <- c("reportData/prelim.norm.summary.RData", "reportData/cel.headers.RData")

data.file <- args[1]
cel.headers.file <- args[2]
library(RColorBrewer)
library("lubridate")

# load data
# ============================================================

load(data.file)
load(cel.headers.file)

rownames(cel.headers.18c) <- cel.headers.18c[, 1]
rownames(cel.headers.25c) <- cel.headers.25c[, 1]

cel.date.18c <- parse_date_time(cel.headers.18c[, 2], "mdy HMS")
cel.date.25c <- parse_date_time(cel.headers.25c[, 2], "mdy HMS")

names(cel.date.18c) <- cel.headers.18c[, 1]
names(cel.date.25c) <- cel.headers.25c[, 1]

female.18c.quantile <- female.18c.quantile[, order(cel.date.18c[colnames(female.18c.quantile)])]
male.18c.quantile <- male.18c.quantile[, order(cel.date.18c[colnames(male.18c.quantile)])]

female.25c.quantile <- female.25c.quantile[, order(cel.date.25c[colnames(female.25c.quantile)])]
male.25c.quantile <- male.25c.quantile[, order(cel.date.25c[colnames(male.25c.quantile)])]

# prepare file to plot
# ============================================================

plot.col = c(brewer.pal(n = 9, name = "Blues")[c(9, 7, 4)], "grey80", brewer.pal(n = 9, name = "Reds")[c(4, 7, 9)])

file.width = 120
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width * 1.6/25.4, family = args[4])
par(las = 1, tcl = -0.2, mar = c(2, 2, 1, 0.9), ps = 7, lwd = 0.5, mfrow = c(4, 1))

plot(c(1, 384), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(female.18c.quantile[i, ], col = plot.col[i])
}

bad.array <- which(female.18c.quantile[1, ] < -6 | female.18c.quantile[7, ] > 6)

if (length(bad.array) > 0) { 
  
  for (i in 1:length(bad.array)) { 
    
    k <- ifelse(abs(female.18c.quantile[1, bad.array[i]]) > abs(female.18c.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], female.18c.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(gsub(".CEL", "", names(bad.array[i])), split = "-"))
    text(bad.array[i], female.18c.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))
  
  }

}
  
legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Females (18 ", {}*degree, "C)")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw date

female.18c.date <- date(cel.date.18c[colnames(female.18c.quantile)])
female.18c.date.unique <- unique(date(cel.date.18c[colnames(female.18c.quantile)]))

for (i in 1:length(female.18c.date.unique)) {

  segments(min(which(female.18c.date == female.18c.date.unique[i])), -15, max(which(female.18c.date == female.18c.date.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(female.18c.date == female.18c.date.unique[i])), ifelse(i %% 2, -14, -12.5), female.18c.date.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  
}


box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)





plot(c(1, 384), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(female.25c.quantile[i, ], col = plot.col[i])
}

bad.array <- which(female.25c.quantile[1, ] < -6 | female.25c.quantile[7, ] > 6)

if (length(bad.array) > 0) { 
  
  for (i in 1:length(bad.array)) { 
    
    k <- ifelse(abs(female.25c.quantile[1, bad.array[i]]) > abs(female.25c.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], female.25c.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(gsub("T", "", gsub(".CEL", "", names(bad.array[i]))), split = "-"))
    text(bad.array[i], female.25c.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))
  
  }

}
  
legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Females (25 ", {}*degree, "C)")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw date

female.25c.date <- date(cel.date.25c[colnames(female.25c.quantile)])
female.25c.date.unique <- unique(date(cel.date.25c[colnames(female.25c.quantile)]))

for (i in 1:length(female.25c.date.unique)) {

  segments(min(which(female.25c.date == female.25c.date.unique[i])), -15, max(which(female.25c.date == female.25c.date.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(female.25c.date == female.25c.date.unique[i])), ifelse(i %% 2, -14, -12.5), female.25c.date.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  
}


box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)





plot(c(1, 384), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(male.18c.quantile[i, ], col = plot.col[i])
}

bad.array <- which(male.18c.quantile[1, ] < -6 | male.18c.quantile[7, ] > 6)

if (length(bad.array) > 0) { 
  
  for (i in 1:length(bad.array)) { 
    
    k <- ifelse(abs(male.18c.quantile[1, bad.array[i]]) > abs(male.18c.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], male.18c.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(gsub(".CEL", "", names(bad.array[i])), split = "-"))
    text(bad.array[i], male.18c.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))
  
  }

}
  
legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Males (18 ", {}*degree, "C)")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw date

male.18c.date <- date(cel.date.18c[colnames(male.18c.quantile)])
male.18c.date.unique <- unique(date(cel.date.18c[colnames(male.18c.quantile)]))

for (i in 1:length(male.18c.date.unique)) {

  segments(min(which(male.18c.date == male.18c.date.unique[i])), -15, max(which(male.18c.date == male.18c.date.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(male.18c.date == male.18c.date.unique[i])), ifelse(i %% 2, -14, -12.5), male.18c.date.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  
}


box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)





plot(c(1, 384), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(male.25c.quantile[i, ], col = plot.col[i])
}

bad.array <- which(male.25c.quantile[1, ] < -6 | male.25c.quantile[7, ] > 6)

if (length(bad.array) > 0) { 
  
  for (i in 1:length(bad.array)) { 
    
    k <- ifelse(abs(male.25c.quantile[1, bad.array[i]]) > abs(male.25c.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], male.25c.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(gsub("T", "", gsub(".CEL", "", names(bad.array[i]))), split = "-"))
    text(bad.array[i], male.25c.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))
  
  }

}
  
legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Males (25 ", {}*degree, "C)")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw date

male.25c.date <- date(cel.date.25c[colnames(male.25c.quantile)])
male.25c.date.unique <- unique(date(cel.date.25c[colnames(male.25c.quantile)]))

for (i in 1:length(male.25c.date.unique)) {

  segments(min(which(male.25c.date == male.25c.date.unique[i])), -15, max(which(male.25c.date == male.25c.date.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(male.25c.date == male.25c.date.unique[i])), ifelse(i %% 2, -14, -12.5), male.25c.date.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  
}


box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
