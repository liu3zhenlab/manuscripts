################################################################################
### mRNA Drought*Day interaction analysis
### Modified by Cheng He
### Sanzhen Liu
### 1/3/2018
################################################################################
### required packages
library(DESeq2)

### RC data (rm 11d data)
d <- read.delim("rc_mrna_b73_rm11d.txt")

### experimental design (rm 11d data)
expdesign <- read.delim("exp.design.leafwater.batch_rm11d.txt")

middle <- expdesign[expdesign$Day >= 6 & expdesign$Day <= 8, ]

### subset (8 time points in each condition)
# minimum total reads
minCounts <- 2*nrow(middle)
d2 <- d[, c("Gene", as.character(middle$Sample))]
### require at least "minCounts" total counts from all the samples
d2 <- d2[rowSums(d2[, -1])>=minCounts, ]
nrow(d2)
head(d2)
### sample design information in the DESeq format
sample.info <- data.frame(row.names=middle$Sample,
                          day=as.factor(middle$Day),
            			  drought=as.factor(middle$Condition),
						  batch=as.factor(middle$Batch))

### DESeq analysis
mRNA <- d2[, 1] ### mRNA sequence
dds <- DESeqDataSetFromMatrix(countData=as.matrix(d2[, -1]),
                              colData=sample.info,
                              formula(~day+drought+day*drought+batch))

### leave effect, including interaction with "day"
dds <- DESeq(dds, "LRT",
             full=formula(~day+drought+day*drought+batch),
             reduced=formula(~day+drought+batch))

### test result
res <- results(dds, cooksCutoff=F,independentFiltering=T, alpha=0.1)
summary(res)
res$mRNA <- as.character(mRNA)
out.coef <- coef(dds, SE=F) ### model coefficients
out.coef <- cbind(as.data.frame(out.coef),mRNA)
out.se <- coef(dds, SE=T) ### standard errors
out.se <- cbind(as.data.frame(out.se),mRNA)
out.disp <- data.frame(mRNA=mRNA, Dispersion=dispersions(dds))

#### filter
sum(!is.na(res$padj))
res2 <- res[!is.na(res$padj), ]
nrow(res2)

fdr.cutoff <- 0.05

sum(!is.na(res$padj) & res$padj<fdr.cutoff)
sig.mRNAs <- res[!is.na(res$padj) & res$padj<fdr.cutoff, "mRNA"]
sig.out.coef <- out.coef[out.coef$mRNA %in% sig.mRNAs, ]
sig.out.se <- out.se[out.se$mRNA %in% sig.mRNAs, ]
sig.out.dispersion <- out.disp[out.disp$mRNA %in% sig.mRNAs, ]

save.image("middle.mRNA.interaction.RData")

### the "alpha" value has been adjusted to obtain more a uniform p-value
### histogram shape at the high p-value side
pdf(file="middle.interaction.hist.pdf",width = 6,height = 5)
hist(res2$pvalue, xlab="p-values", ylab="Number of mRNAs", main="interaction effect")
dev.off()

### plotMA: red points represent genes with padj < 0.05
pdf(file="middle.interaction.plotma.pdf",width = 6,height = 5)
plotMA(results(dds, cooksCutoff=F,independentFiltering=T, alpha=0.05), ylim=c(-15,15))
dev.off()

### output
write.table(res2, "middle.mRNA.interaction.raw.output.txt", quote=F, row.names=F, sep="\t")
write.table(sig.out.coef, "middle.sig.interaction.mRNAs.coefficients.output.txt",
            quote=F, row.names=F, sep="\t")
write.table(sig.out.se, "middle.sig.interaction.mRNAs.se.output.txt",
            quote=F, row.names=F, sep="\t")
write.table(sig.out.dispersion, "middle.sig.interaction.mRNAs.dispersion.output.txt",
            quote=F, row.names=F, sep="\t")

