# =================================================
# = summarize TF bindign site permutation results =
# =================================================

args <- commandArgs(TRUE) # args <- c("eqtl/tf.site.list", "eqtl/male.tf.perm.count", 10000)

tf.list <- scan(args[1], what = "")
tf.perm.count <- read.table(args[2], header = FALSE, as.is = TRUE)
n.perm <- as.numeric(args[3])

tf.fold.enrich <- matrix(nrow = length(tf.list), ncol = 3)
tf.fold.enrich.p <- matrix(nrow = length(tf.list), ncol = 3)

# loop through all tfs
# ============================================================

for (i in 1:length(tf.list)) {
  
  tf <- tf.list[i]
  
  this.tf.data <- tf.perm.count[tf.perm.count[, 1] == tf, ]
  this.tf.obs <- this.tf.data[this.tf.data[, 3] == "obs", ]
  this.tf.obs.count <- c(sum(this.tf.obs[, 2] == "25c"), (sum(this.tf.obs[, 2] == "both") + 0.5)/(sum(this.tf.obs[, 2] == "25c") + 0.5), (sum(this.tf.obs[, 2] == "18c") + 0.5)/(sum(this.tf.obs[, 2] == "25c") + 0.5))
  
  this.tf.perm.count <- matrix(nrow = n.perm, ncol = 3)
  
  this.tf.perm.count[, 1] <- sapply(split(this.tf.data[this.tf.data[, 3] != "obs", 2], factor(this.tf.data[this.tf.data[, 3] != "obs", 3], levels = paste("perm", 1:n.perm, sep = ""))), function(x) { return (sum(x == "25c"))})
  this.tf.perm.count[, 2] <- sapply(split(this.tf.data[this.tf.data[, 3] != "obs", 2], factor(this.tf.data[this.tf.data[, 3] != "obs", 3], levels = paste("perm", 1:n.perm, sep = ""))), function(x) { return ((sum(x == "both") + 0.5)/(sum(x == "25c") + 0.5)) })
  this.tf.perm.count[, 3] <- sapply(split(this.tf.data[this.tf.data[, 3] != "obs", 2], factor(this.tf.data[this.tf.data[, 3] != "obs", 3], levels = paste("perm", 1:n.perm, sep = ""))), function(x) { return ((sum(x == "18c") + 0.5)/(sum(x == "25c") + 0.5)) })
  
  tf.fold.enrich[i, 1] <- this.tf.obs.count[1]/mean(this.tf.perm.count[, 1])
  tf.fold.enrich[i, 2] <- this.tf.obs.count[2]/mean(this.tf.perm.count[, 2])
  tf.fold.enrich[i, 3] <- this.tf.obs.count[3]/mean(this.tf.perm.count[, 3])
  
	if (tf.fold.enrich[i, 1] >= 1) {
		tf.fold.enrich.p[i, 1] <- sum(this.tf.perm.count[, 1] >= this.tf.obs.count[1])/n.perm
	} else {
		tf.fold.enrich.p[i, 1] <- sum(this.tf.perm.count[, 1] < this.tf.obs.count[1])/n.perm
	}
	
	if (tf.fold.enrich[i, 2] >= 1) {
	  tf.fold.enrich.p[i, 2] <- sum(this.tf.perm.count[, 2] >= this.tf.obs.count[2])/n.perm
	} else {
	  tf.fold.enrich.p[i, 2] <- sum(this.tf.perm.count[, 2] < this.tf.obs.count[2])/n.perm
	}
	
	if (tf.fold.enrich[i, 3] >=1) {
	  tf.fold.enrich.p[i, 3] <- sum(this.tf.perm.count[, 3] >= this.tf.obs.count[3])/n.perm
	} else {
	  tf.fold.enrich.p[i, 3] <- sum(this.tf.perm.count[, 3] < this.tf.obs.count[3])/n.perm
	}
  
}

save(tf.list, tf.fold.enrich, tf.fold.enrich.p, file = args[4])

