################################################################################
### mRNA Drought*Day interaction analysis
### Modified by Cheng He
### Sanzhen Liu
### 1/2/2018
################################################################################
### required packages
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

### RC data (rm 11d data)
d <- read.delim("rc_mrna_b73_rm11d.txt")
### experimental design (rm 11d data)
expdesign <- read.delim("exp.design.leafwater.batch_rm11d.txt")

### subset (8 time points in each condition)
# minimum total reads
minCounts <- 64 ## average 2 read counts each sample
d2 <- d[, c("Gene", as.character(expdesign$Sample))]
### require at least "minCounts" total counts from all the samples
d2 <- d2[rowSums(d2[, -1])>=minCounts, ]
nrow(d2)

### sample design information in the DESeq format
sample.info <- data.frame(row.names=expdesign$Sample,
                          day=as.factor(expdesign$Day),
                  			  drought=as.factor(expdesign$Condition),
      			      			  batch=expdesign$Batch)

### DESeq analysis
mRNA <- d2[, 1] ### mRNA sequence
dds <- DESeqDataSetFromMatrix(countData=as.matrix(d2[, -1]),
                              colData=sample.info,
                              formula(~ day + drought + batch + day*drought))

### leave effect, including interaction with "day", DEGs in different models
dds <- DESeq(dds, "LRT",
             full=formula(~ day + drought + batch + day*drought),
             reduced=formula(~ day + drought + batch))

### test result
res <- results(dds, cooksCutoff=F,independentFiltering=T, alpha=0.1)
summary(res)
res$mRNA <- as.character(mRNA)
out.coef <- coef(dds, SE=F) ### model coefficients
out.coef <- cbind(as.data.frame(out.coef),mRNA)
out.coef[1,]
out.se <- coef(dds, SE=T) ### standard errors
out.se <- cbind(as.data.frame(out.se),mRNA)
out.se[1,]
out.disp <- data.frame(mRNA=mRNA, Dispersion=dispersions(dds))
out.disp[1,]

#### filter
sum(!is.na(res$padj))
res2 <- res[!is.na(res$padj), ]
nrow(res2)
summary(res2)

fdr.cutoff <- 0.01

sum(!is.na(res$padj) & res$padj<fdr.cutoff)
sig.mRNAs <- res[!is.na(res$padj) & res$padj<fdr.cutoff, "mRNA"]
summary(res[!is.na(res$padj) & res$padj<fdr.cutoff,])
sig.out.coef <- out.coef[out.coef$mRNA %in% sig.mRNAs, ]
sig.out.se <- out.se[out.se$mRNA %in% sig.mRNAs, ]
sig.out.dispersion <- out.disp[out.disp$mRNA %in% sig.mRNAs, ]

save.image("mRNA.interaction.RData")

### the "alpha" value has been adjusted to obtain more a uniform p-value
### histogram shape at the high p-value side
hist(res2$pvalue, ylim=c(0,12000),xlab="p-values", ylab="Number of mRNAs", main="interaction effect",col="gray")

### plotMA: red points represent genes with padj < 0.01
plotMA(results(dds, cooksCutoff=F,independentFiltering=T, alpha=0.01), ylim=c(-15,15))

### output
write.table(res2, "mRNA.interaction.raw.output.txt", quote=F, row.names=F, sep="\t")  ## raw mRNA interaction output without NA
write.table(sig.out.coef, "sig.interaction.mRNAs.coefficients.output.txt",
            quote=F, row.names=F, sep="\t")                                           ## coefficient of mRNAs with padj < 0.01
write.table(sig.out.se, "sig.interaction.mRNAs.se.output.txt",
            quote=F, row.names=F, sep="\t")                                           ## standard error of mRNAs with padj < 0.01
write.table(sig.out.dispersion, "sig.interaction.mRNAs.dispersion.output.txt",
            quote=F, row.names=F, sep="\t")                                           ## dispersion of mRNAs with padj < 0.01

