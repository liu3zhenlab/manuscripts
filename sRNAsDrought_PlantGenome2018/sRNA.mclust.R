################################################################################
### Clustering analysis
### Sanzhen Liu and Haiyan Wang
### Kansas State University
### 10/5/2017
################################################################################
### experimental design
expdesign <- read.delim("exp.design.txt")
expdesign$CondDay <- paste0(expdesign$Condition, expdesign$Day)

################################################################################
### 1. pre-clustering: organize input data
################################################################################
signorm <- read.delim("sig.sRNA.normalized.counts.txt", stringsAsFactors = F)
signorm <- signorm[, as.character(expdesign$Sample)]
signorm2 <- (signorm[, 1:16*2] + signorm[, 1:16*2-1])/2
head(signorm2)
colnames(signorm2) <- gsub("\\.[0-9]+", "", colnames(signorm2))
rownames(signorm2) <- signorm$sRNA  ### label rows with sRNAs
### log2 ratios, 1 was added to avoid 0 counts
ratios <- log2((signorm2[, 1:8] + 1) / (signorm2[, 9:16] + 1))

################################################################################
### clustering
################################################################################
library(mclust)
mc <- Mclust(ratios, G=1:20, modelNames ='VVV') ### why do we select "VVV"
mc.summary <- summary(mc, parameter = TRUE, classification = TRUE)
plot(mc,'BIC')
title("Relative expression")
G <- unlist(mc.summary[4])
meancurve <- matrix(unlist(mc.summary[11]), ncol=G)
covariance.between.time.points <- mc.summary[12]
class.member <- data.frame(mc.summary[15])
rownames(meancurve) <- colnames(ratios)
colnames(meancurve) <- paste('Group', 1:G, sep="")

