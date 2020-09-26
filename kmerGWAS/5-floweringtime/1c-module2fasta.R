setwd("/bulk/liu3zhen/research/projects/kGWAS/ft50k/2-CON/")
options(stringsAsFactors=F)
con <- read.delim("../0-data_fromCheng/CON_DTS_p_2.326453e-08_P14M150_36_modules.txt")
con <- con[, c(1,2,3)]
colnames(con) <- c("id", "CON", "pvalue")

conmodules <- unique(con$CON)
conmodules <- conmodules[conmodules != "grey"]

for (em in conmodules) {
  kmers <- con[con$CON==em, 1]
  pvalues <- con[con$CON==em, 3]
  outfile <- paste0("1o-clust_", em, ".fasta")
  outdata <- paste0(">k", 1:length(kmers), "_", pvalues, "\n", kmers, "\n")
  cat(outdata, sep="", file=outfile)
}
