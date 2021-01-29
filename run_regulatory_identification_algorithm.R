DATA_PATH <- "./data"
RESULT_PATH <- "./result"

cohorts <- c("KIRC", "KIRP")
set.seed(623)

###load first cohort
load(sprintf("%s/TCGA-%s.RData", DATA_PATH, toupper(cohorts[1])))  
common_patients <- intersect(rownames(TCGA$mirna), rownames(TCGA$mrna))

X1 <- TCGA$mirna[common_patients,]
X1 <- X1[, which(colMeans(X1 > 0) >= 0.5)]
X1 <- scale(log2(X1 + 1))

Y1 <- TCGA$mrna[common_patients,]
Y1 <- Y1[, which(colMeans(Y1 > 0) >= 0.5)]
Y1 <- scale(log2(Y1 + 1))

###load second cohort
load(sprintf("%s/TCGA-%s.RData", DATA_PATH, toupper(cohorts[2])))
common_patients <- intersect(rownames(TCGA$mirna), rownames(TCGA$mrna))

X2 <- TCGA$mirna[common_patients,]
X2 <- X2[, which(colMeans(X2 > 0) >= 0.5)]
X2 <- scale(log2(X2 + 1))

Y2 <- TCGA$mrna[common_patients,]
Y2 <- Y2[, which(colMeans(Y2 > 0) >= 0.5)]
Y2 <- scale(log2(Y2 + 1))

###filter shared miRNAs and mRNAs
common_miRNAs <- intersect(colnames(X1), colnames(X2))
common_mRNAs <- intersect(colnames(Y1), colnames(Y2))
X1 <- X1[, common_miRNAs]
Y1 <- Y1[, common_mRNAs]
X2 <- X2[, common_miRNAs]
Y2 <- Y2[, common_mRNAs]

##########run the algorithm
ranks <- list()
ranks$R1_u <- 10 # upper bound for first cohort factors 
ranks$R2_u <- 10 # upper bound for second cohort factors
ranks$Rs_u <- 10 # upper bound for shared factors

lam1 <- 300
lam2 <- 25
lam3 <- 300
lam4 <- 30

source("MSRFR.R")
out <- MSRFR(X1, X2, Y1, Y2, lam1, lam2, lam3, lam4, ranks)
Wx1 <- as.matrix(data.frame(out[1]))
Wx2 <- as.matrix(data.frame(out[2]))
Wsx <- as.matrix(data.frame(out[3]))
Wy1 <- as.matrix(data.frame(out[4]))
Wy2 <- as.matrix(data.frame(out[5]))
Wsy <- as.matrix(data.frame(out[6]))
iteration <- data.frame(out[8])
obj_vals <- data.frame(out[9])
ranks <- data.frame(out[7])

######normalization
colnames(Wy1) <- colnames(Y1, do.NULL = TRUE)
colnames(Wy2) <- colnames(Y2, do.NULL = TRUE)
colnames(Wsy) <- colnames(Y1, do.NULL = TRUE)
rownames(Wx1)<-colnames(X1,do.NULL = TRUE)
rownames(Wx2)<-colnames(X2,do.NULL = TRUE) 
rownames(Wsx)<-colnames(X1,do.NULL = TRUE) 
program_norms_normalized_1 <- sqrt(colSums(Wx1^2) * rowSums(Wy1^2))
program_norms_normalized_2 <- sqrt(colSums(Wx2^2) * rowSums(Wy2^2))
program_norms_normalized_s <- sqrt(colSums(Wsx^2) * rowSums(Wsy^2))
sorted_indices_normalized_1 <- sort(program_norms_normalized_1, decreasing = TRUE, index.return = TRUE)$ix
sorted_indices_normalized_2 <- sort(program_norms_normalized_2, decreasing = TRUE, index.return = TRUE)$ix
sorted_indices_normalized_s <- sort(program_norms_normalized_s, decreasing = TRUE, index.return = TRUE)$ix

summary <- list()
summary$weights_normalized_1 <- program_norms_normalized_1[sorted_indices_normalized_1]
summary$weights_normalized_2 <- program_norms_normalized_2[sorted_indices_normalized_2]
summary$weights_normalized_3 <- program_norms_normalized_s[sorted_indices_normalized_s]
summary$Wx1_normalized <- Wx1[,sorted_indices_normalized_1]
summary$Wx2_normalized <- Wx2[,sorted_indices_normalized_2]
summary$Wsx_normalized <- Wsx[,sorted_indices_normalized_s]
for (r in 1:ncol(Wx1)) {
  summary$Wx1_normalized[,r] <- summary$Wx1_normalized[,r] / sqrt(sum(summary$Wx1_normalized[,r]^2))
}
for (r in 1:ncol(Wx2)) {
  summary$Wx2_normalized[,r] <- summary$Wx2_normalized[,r] / sqrt(sum(summary$Wx2_normalized[,r]^2))
}
for (r in 1:ncol(Wsx)) {
  summary$Wsx_normalized[,r] <- summary$Wsx_normalized[,r] / sqrt(sum(summary$Wsx_normalized[,r]^2))
}
rownames(summary$Wx1_normalized) <- rownames(Wx1)
rownames(summary$Wx2_normalized) <- rownames(Wx2)
rownames(summary$Wsx_normalized) <- rownames(Wsx)
summary$Wy1_normalized <- Wy1[sorted_indices_normalized_1,]
summary$Wy2_normalized <- Wy2[sorted_indices_normalized_2,]
summary$Wsy_normalized <- Wsy[sorted_indices_normalized_s,]
for (r in 1:nrow(Wy1)) {
  summary$Wy1_normalized[r,] <- summary$Wy1_normalized[r,] / sqrt(sum(summary$Wy1_normalized[r,]^2))
}
for (r in 1:nrow(Wy2)) {
  summary$Wy2_normalized[r,] <- summary$Wy2_normalized[r,] / sqrt(sum(summary$Wy2_normalized[r,]^2))
}
for (r in 1:nrow(Wsy)) {
  summary$Wsy_normalized[r,] <- summary$Wsy_normalized[r,] / sqrt(sum(summary$Wsy_normalized[r,]^2))
}
colnames(summary$Wy1_normalized) <- colnames(Wy1)
colnames(summary$Wy2_normalized) <- colnames(Wy2)
colnames(summary$Wsy_normalized) <- colnames(Wsy)

#######selection schema
summary$Wx1_thresholded <- Wx1[,sorted_indices_normalized_1] * (abs(summary$Wx1_normalized) > 2 / sqrt(nrow(summary$Wx1_normalized)))
summary$Wy1_thresholded <- Wy1[sorted_indices_normalized_1,] * (abs(summary$Wy1_normalized) > 2 / sqrt(ncol(summary$Wy1_normalized)))
summary$Wx2_thresholded <- Wx2[,sorted_indices_normalized_2] * (abs(summary$Wx2_normalized) > 2 / sqrt(nrow(summary$Wx2_normalized)))
summary$Wy2_thresholded <- Wy2[sorted_indices_normalized_2,] * (abs(summary$Wy2_normalized) > 2 / sqrt(ncol(summary$Wy2_normalized)))
summary$Wsx_thresholded <- Wsx[,sorted_indices_normalized_s] * (abs(summary$Wsx_normalized) > 2 / sqrt(nrow(summary$Wsx_normalized)))
summary$Wsy_thresholded <- Wsy[sorted_indices_normalized_s,] * (abs(summary$Wsy_normalized) > 2 / sqrt(ncol(summary$Wsy_normalized)))

########### Summary 
summary$Wy1 <- Wy1
summary$Wy2 <- Wy2
summary$Wsy <- Wsy
summary$Wx1 <- Wx1
summary$Wx2 <- Wx2
summary$Wsx <- Wsx
summary$iteration <- iteration
summary$obj_vals <- obj_vals

save(summary, file = sprintf("%s/MSRFR_summary.RData", RESULT_PATH))

