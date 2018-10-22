################################################################################
### sRNA Drought*Day interaction analysis
### Sanzhen Liu
### 10/5/2017
################################################################################
### required packages
library(DESeq2)

### threshold
fdr_cutoff <- 0.05

### data
d <- read.delim("sRNA.counts.txt", stringsAsFactors = F)  ### sRNA counting data

### experimental design
expdesign <- read.delim("exp.design.txt")

### sample design information in the DESeq format
sample.info <- data.frame(row.names=expdesign$Sample,
                          condition=expdesign$Condition,
                          day=as.factor(expdesign$Day))

### DESeq analysis
sRNA <- d[, 1] ### sRNA sequence (reverse complementary)
dds <- DESeqDataSetFromMatrix(countData=as.matrix(d[, -1]),
                              colData=sample.info,
                              formula(~ condition + day + condition*day))

dds <- DESeq(dds, "LRT", full=formula(~ condition + day + condition*day),
             reduced=formula(~ condition + day))

### result
### the "alpha" value has been adjusted to obtain more a uniform p-value
res <- results(dds, cooksCutoff=F, independentFiltering=T, alpha=0.1)  ### test result
sum(!is.na(res$padj) & res$padj < fdr_cutoff)

