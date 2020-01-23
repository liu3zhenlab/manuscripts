x <- read.table("B73_drought_GWAS_368lines_16338snps_popstructure_110607_+_maize368_drought_phe_for_tassel_+_368_55w_v3_imp_stats.txt",sep="\t",header=TRUE)
pos <- read.table("gene_pos_v3.txt",sep="\t",header=TRUE)
head(x)
head(pos)
c <- c()
i <- 1
for (i in 1:length(pos[,1]))
{
  snp <- x[which(x$Chr==pos$chr[i] & x$Pos >= (pos$start[i]-1000) & x$Pos < (pos$end[i]+1000)),]
  if (length(snp$p) > 0)
  {
    p <- min(snp$p)
    c <- c(c,p)
  }else
  {
    c <- c(c,"NA")
  }
}
out <- cbind(pos,c)
names(out)[5] <- "min_pvalue"
write.table(out,file="B73_drought_GWAS_gene_1k.txt",sep="\t",quote=FALSE,row.names=FALSE)