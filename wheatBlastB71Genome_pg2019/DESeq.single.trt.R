# analysis for a single treatment RNA-seq data:
DESeq.single.trt <- function (input.matrix, comparison,
                             fdr=0.05, min.mean.reads=5,
                             min.positive.sample=2) {
  # Use the package of DESeq
  library(DESeq)
  # input.matrix: see the example below, row names are gene IDs 
  #               the input matrix should have the header each of 
  #               which matches one of element of the comparison
  # comparison: a vector of comparison information, e.g., c("w","m"), 
  # if group1.col and group2.col are not specified, these columns will be automatically
  # determined by "comparison". In this case, the column names should contain the features
  # specified in the comparison
  # fdr: false discovery rate
  # min.mean.reads: the minimum average reads required across all samples
  # min.positive.sample: minimum number of samples with positive read counts
  
  # treatment: a vector of treatment information
  if (!is.null(group1.col)) {
  	group1.col <- grep(comparison[1],colnames(input.matrix))
  }

  if (!is.null(group2.col)) {
  	group2.col <- grep(comparison[2],colnames(input.matrix))
  }

  repnum1 <- length(group1.col)
  repnum2 <- length(group2.col)
  treatment <- c(rep(1,repnum1),rep(2,repnum2))
  indata <- input.matrix[, c(group1.col, group2.col)]
  # filter low expressed genes:
  count.positive <- function (x) { sum(x>0, na.rm=T) }
  indata <- indata[apply(indata, 1, count.positive) >= min.positive.sample &
    apply(indata, 1, sum) >= min.mean.reads, ] # filtering
  # experimental design:
  design <- data.frame(row.names=colnames(indata), condition=treatment)
  condition <- factor(design$condition)
  cds <- newCountDataSet(indata, condition)
  # normalization:
  cds <- estimateSizeFactors(cds)
  sizeFactors(cds)
  #normalizedcount <- counts(cds, normalized=TRUE)
  cds <- estimateDispersions(cds)

  # DE:
  res <- nbinomTest(cds, "1", "2")
  res <- res[,c("id", "log2FoldChange", "pval", "padj")]
  cp <- paste(comparison[2],"-",comparison[1],"_",sep="")
  colnames(res) <- c("GeneID", 
                     paste(cp,"log2FC",sep=""),
                     paste(cp,"pvalue",sep=""),
                     paste(cp,"qvalue",sep=""))
  res$GeneID <- rownames(indata)
  res
}

