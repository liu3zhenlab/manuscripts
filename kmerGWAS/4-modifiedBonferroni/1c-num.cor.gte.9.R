setwd("/bulk/liu3zhen/research/projects/kGWAS/modifiedBonferroni")
args <- commandArgs(trailingOnly=T)
first <- args[1]
last <- args[2]
cat("first:", first, "\n", sep="")
cat("last:", last, "\n", sep="")

d <- read.delim("DTS_random_pick_1M_kmer_RC.txt", stringsAsFactors=F)
d <- d[, -1]
dp <- d[first:last, ]

################################################
corfun <- function (x, data, cor.cutoff=0.9) {
# one to the whole set and return # of cor greater than the cutoff
# the correlation to itself was excluded.
  x <- as.numeric(x)
  cor.vals <- apply(data, 1, cor, y=x)
  sum(cor.vals >= cor.cutoff) - 1
}
################################################

out <- apply(dp, 1, corfun, data=d)

outdf <- data.frame(NumHighCor=out)

# outfile
outfile <- paste0("1o-num.cor.gte.9.", first, "-", last)
write.table(outdf, outfile, row.names=F, quote=F, sep="\t")

