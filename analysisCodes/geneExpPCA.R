# ================================
# = PCA for gene expression data =
# ================================

args <- commandArgs(TRUE) # args <- c("tilingInt/female.gene.exp.RData", "female.gene.exp", "tilingInt/cel.headers.RData")
rdata.file <- args[1]
data.set <- args[2]
library("lubridate")

# load data
# ============================================================

load(rdata.file)
eval(parse(text = paste("gene.exp = ", data.set, sep = "")))
load(args[3])

nqt <- function(x) {
  
  return(qnorm(rank(x)/(length(x) + 1), sd = mad(x)) + median(x))
  
}

# date
# ============================================================

rownames(cel.headers.18c) <- cel.headers.18c[, 1]
rownames(cel.headers.25c) <- cel.headers.25c[, 1]

cel.date.18c <- parse_date_time(cel.headers.18c[, 2], "mdy HMS")
cel.date.25c <- parse_date_time(cel.headers.25c[, 2], "mdy HMS")

names(cel.date.18c) <- cel.headers.18c[, 1]
names(cel.date.25c) <- cel.headers.25c[, 1]

# find line id and temperature
# ============================================================

temp <- rep("25C", ncol(gene.exp))
temp[grep("18C", colnames(gene.exp))] <- "18C"
temp <- factor(temp, levels = c("25C", "18C"))
line.id <- factor(paste("line_", gsub("T", "", gsub("-.*", "", colnames(gene.exp))), sep = ""))
scan.date <- date(c(cel.date.18c, cel.date.25c)[colnames(gene.exp)])

# overall PCA across temperature
# ============================================================

gene.exp.nqt <- matrix(NA, ncol = ncol(gene.exp), nrow = nrow(gene.exp))
for (i in 1:nrow(gene.exp)) { gene.exp.nqt[i, ] <- nqt(gene.exp[i, ]) }
colnames(gene.exp.nqt) <- colnames(gene.exp)
all.pca <- prcomp(t(gene.exp.nqt), scale = T)

# PCA within each temperature
# ============================================================

gene.exp.18c <- gene.exp[, temp == "18C"]
gene.exp.18c.nqt <- gene.exp
for (i in 1:nrow(gene.exp.18c)) { gene.exp.18c.nqt[i, ] <- nqt(gene.exp.18c[i, ]) }

gene.exp.25c <- gene.exp[, temp == "25C"]
gene.exp.25c.nqt <- gene.exp
for (i in 1:nrow(gene.exp.25c)) { gene.exp.25c.nqt[i, ] <- nqt(gene.exp.25c[i, ]) }

pca.18c <- prcomp(t(gene.exp.18c.nqt), scale = T)
pca.25c <- prcomp(t(gene.exp.25c.nqt), scale = T)

all.pca.top2 <- predict(all.pca)[, 1:2]
pca.18c.top2 <- predict(pca.18c)[, 1:2]
pca.25c.top2 <- predict(pca.25c)[, 1:2]

all.pca.eigen <- all.pca$sdev^2
pca.18c.eigen <- pca.18c$sdev^2
pca.25c.eigen <- pca.25c$sdev^2

save(all.pca.top2, pca.18c.top2, pca.25c.top2, all.pca.eigen, pca.18c.eigen, pca.25c.eigen, temp, line.id, scan.date, cel.date.18c, cel.date.25c, file = args[4])

# session info
# ============================================================

sessionInfo()
