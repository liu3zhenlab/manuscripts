################################################################################
### Clustering analysis
### Haiyan Wang and Sanzhen
### Modified by Cheng
### 01/04/2018
################################################################################
### experimental design
expdesign <- read.delim("exp.design.leafwater.batch_rm11d.txt")

### sample design information in the DESeq format
mm <- model.matrix(~as.factor(expdesign$Day) + as.factor(expdesign$Condition)
                   + as.factor(expdesign$Day)*as.factor(expdesign$Condition) + expdesign$Batch)
dim(mm)
### coefficients from the modeling
beta <- read.delim("../2a_mrna_interaction/results/sig.interaction.mRNAs.coefficients.output.txt")
head(beta)
xbeta.raw <- as.matrix(beta[,-ncol(beta)]) %*% as.matrix(t(mm)) ### remove the last col of beta (6022*17) and time with t(mm) (17*32) to form a matrix xbeta.raw (6022*32)
colnames(xbeta.raw) <- expdesign$Sample
head(xbeta.raw)
xbeta <- (xbeta.raw[, 2*(1:16) - 1]  + xbeta.raw[, 2*(1:16)])/2  ### remove duplication by calculate the mean of 2 replicates
colnames(xbeta) <- gsub("\\.[0-9]+", "", colnames(xbeta))
colnames(xbeta)
head(xbeta)

### log2(q) = xbeta
### genometric mean of q
ratios <- xbeta[, 1:8] - xbeta[, 9:16]
nrow(ratios)
head(ratios)
rownames(ratios) <- beta$mRNA
################################################################################
### clustering
################################################################################
library(mclust)
mc_select <- mclustBIC(ratios) ### perform cluster with all models
pdf("./figures/mclust_model_select.pdf", width=9, height=9)
plot(mc_select) ### model with the highest BIC should be the best cluster model
dev.off()

mc <- Mclust(ratios, G=1:15, modelNames ='VVV') ### G is the total number of components
mc.summary <- summary(mc, parameter = TRUE, classification = TRUE)

pdf("./figures/mclust_BIC.pdf", width=9, height=9)
plot(mc,"BIC")
dev.off()

G <- unlist(mc.summary[4]) ### number of cluster components
head(G)
meancurve <- matrix(unlist(mc.summary[11]), ncol=G)
covariance.between.time.points <- mc.summary[12]
class.member <- data.frame(mc.summary[15])

rownames(meancurve) <- 3:10 ### day of drought
colnames(meancurve)=paste('Group', 1:G, sep="") ### cluster groups (components)

write.table(meancurve, file="meanXbeta.Ratio.Curve.data.txt")

write.table(data.frame(mRNA=rownames(class.member), group=class.member[,1]),
            file="clustering.relative.ratio.txt", row.names=F, sep="\t", quote=F)


### extract data
d <- read.delim("clustering.relative.ratio.txt")
head(d)
pdf("./figures/mRNA.interaction.cluster_cheng.pdf", width=9, height=9)
par(mar=c(2, 2, 1.5, 0.5), mfrow=c(4,3))
for (group.num in 1:9) {  ### as much as number of components
  mRNAs <- d[d$group==group.num, "mRNA"]
  counts <- length(mRNAs)
  ### yaxis range
  ymin <- min(ratios[mRNAs, ], na.rm=T)
  ymax <- max(ratios[mRNAs, ], na.rm=T)
  
  ### plotting
  plot(NULL, NULL, xlim=c(3,10), ylim=c(ymin, ymax), cex.main=0.85,
       main=paste0("group", group.num, "; N=", counts),
       cex.main=1)

  linecol <- "light grey"
  ### dominant lines
  for (eachmRNA in mRNAs) {
    lines(3:10, ratios[eachmRNA, 1:8],
          col=linecol, lwd=0.1)
  }
  
  Alld.mean <- apply(ratios[mRNAs, ], 2, mean)
  lines(3:10, Alld.mean, col="red", lwd=1.5)
}

dev.off()


### manually grouping
pdf("./figures/mRNA.interaction.cluster.2major.groups_cheng.pdf", width=6, height=5)

up.group <- c(3, 4, 9)
down.group <- c(5, 6, 7)
other.group <- c(1, 2, 8)
par(mar=c(2, 2, 1.5, 0.5), mfrow=c(3,2))
for (group.num in c(up.group, down.group)) {

  linecol <- "lavenderblush3"
  if (group.num %in% up.group) {
    linecol <- "light blue"
  } else {
    if (group.num %in% down.group) {
      linecol <- c("navajowhite")
    }
  }
  mRNAs <- d[d$group==group.num, "mRNA"]
  counts <- length(mRNAs)
  ### yaxis range
  ymin <- min(ratios[mRNAs, ], na.rm=T)
  ymax <- max(ratios[mRNAs, ], na.rm=T)
  
  ### plotting
  plot(NULL, NULL, xlim=c(3,10), ylim=c(ymin, ymax), cex.main=0.85,
       main=paste0("group", group.num, "; N=", counts),
       cex.main=1)
  
  ### dominant lines
  for (eachmRNA in mRNAs) {
    lines(3:10, ratios[eachmRNA, 1:8],
          col=linecol, lwd=0.1)
  }
  
  Alld.mean <- apply(ratios[mRNAs, ], 2, mean)
  lines(3:10, Alld.mean, col="red", lwd=1.5)
}

dev.off()
