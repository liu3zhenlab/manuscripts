#################################################################################################
### sRNA_target_mRNAs
### Sanzhen Liu
### Modified by Cheng
### 1/06/2018
#################################################################################################

### module
source("enrich.test.R")
### packages
library(GO.db)
### function to extract GO terms:
term.fun <- function (term) {
  ### library GO.db was required
  term.find.test <- try(Term(GOTERM[[term]]), silent=T)
  if (!is(term.find.test, "try-error")) {
    term.text <- Term(GOTERM[[term]])
  } else {
    term.text <- NA
  }
  return(term.text)
}


###
### data
###
pairs <- read.delim("./inputs/psRNATarge_sRNA_merge_miRNA_pairs.txt")
head(pairs)

sRNAmt <- read.delim("./inputs/all.gt72.sRNA.sub.spRNATarget_cheng.txt")
colnames(sRNAmt)
table(sRNAmt$TimeSeriesGroup)

genemt <- read.delim("./inputs/mRNA_drought_timeseries_output1_cheng.txt")
dim(genemt)
head(genemt)
genemt <- as.data.frame(genemt)
is.list(genemt)
#############################################################################################
sRNA_target_genes <- function (group) {
  #group <- "DOWNonWS"
  gsRNAs <- sRNAmt[!is.na(sRNAmt$TimeSeriesGroup) & sRNAmt$TimeSeriesGroup == group, "sRNA"]
  cat("time-series sRNA group is ", group, "\n")
  cat("Number of sRNAs of this group is ", length(gsRNAs), "\n")

  gsRNAs.genes <- pairs[pairs$sRNA %in% gsRNAs, "GeneID"]  ### sRNA-targeted genes
  gsRNAs.target.gt30genes <- sum(gsRNAs.genes == ">30")
  cat("Number of sRNAs targeting >30 genes: ", gsRNAs.target.gt30genes, "\n")
  gsRNAs.uniq.genes <- unique(gsRNAs.genes[gsRNAs.genes != ">30"])
  cat("Number of unique genes targeted by this group of sRNAs: ", length(gsRNAs.uniq.genes), "\n")

  uniq.genes.stat <- data.frame(table(genemt[genemt$GeneID %in% gsRNAs.uniq.genes, "RegulationOnDS"]))
  colnames(uniq.genes.stat) <- c("GenePattern", "Number of genes")
  print(uniq.genes.stat)
  
  as.character(gsRNAs.uniq.genes)
}
#############################################################################################

down.sRNA2genes <- sRNA_target_genes("DOWNonWS")
up.sRNA2genes <- sRNA_target_genes("UPonWS")
shock.sRNA2genes <- sRNA_target_genes("SHOCKonWS")


#####################
### GO analysis: down
#####################
go <- read.delim("./inputs/ZmB73_v3.gene2go.txt")
gf <- data.frame(GeneID=unique(go$genes))
gf$down.sRNA.target <- "no"
gf$down.sRNA.target[gf$GeneID %in% down.sRNA2genes] <- "down.sRNA.target"
table(gf$down.sRNA.target) 
go.result <- enrich.test(gene2feature=gf, gene2cat=go,
                         gene2bias=NULL, test.feature="down.sRNA.target",
                         pval.cutoff=0.01, enrich.cutoff=2,
                         method="Sampling", sampling=5000, outsave=T,
                         outfile="./enrichments/down.sRNA.targeted.genes.GO.enrichment.txt")

sigresult <- go.result$down.sRNA.target[go.result$down.sRNA.target$pvals<0.01 & go.result$down.sRNA.target$enrich>2, ]
nrow(sigresult)
sigresult$Term <- sapply(as.character(sigresult$Cat), term.fun)
sigresult
write.table(sigresult, "down.psRNATarget.enriched.GO.genes_cheng.txt", sep="\t", quote=F, row.names=F)

######################
### barplot GO enrichment ter of down.sRNA.targeted.genes
######################
pdf("down.psRNATarget.genes.sig.GO.enrich_cheng.pdf", width=6, height=6)

plotdata <- sigresult$Cat_in_down.sRNA.target_gene_num
names(plotdata) <- sigresult$Term
plotdata <- plotdata[!is.na(names(plotdata))]
par(mar=c(4,18,1,1))
barplot(plotdata, horiz=T, las=2, col="cornflowerblue", axes=F,cex.names=0.6, 
        xlab="Number of significant genes", cex.lab=1)
abline(v=0)
#axis(side=1, at = seq(0, max(plotdata), by=round(max(plotdata)/6, 0)))
axis(side=1, at =c((0:10)*10))

dev.off()




###################
### GO analysis: up
###################
go <- read.delim("./inputs/ZmB73_v3.gene2go.txt")
gf <- data.frame(GeneID=unique(go$genes))
gf$up.sRNA.target <- "no"
gf$up.sRNA.target[gf$GeneID %in% up.sRNA2genes] <- "up.sRNA.target"
table(gf$up.sRNA.target) 
go.result <- enrich.test(gene2feature=gf, gene2cat=go,
                         gene2bias=NULL, test.feature="up.sRNA.target",
                         pval.cutoff=0.01, enrich.cutoff=2,
                         method="Sampling", sampling=5000, outsave=T,
                         outfile="./enrichments/up.sRNA.targeted.genes.GO.enrichment.txt")

sigresult <- go.result$up.sRNA.target[go.result$up.sRNA.target$pvals<0.01 & go.result$up.sRNA.target$enrich>2, ]
nrow(sigresult)
sigresult$Term <- sapply(as.character(sigresult$Cat), term.fun)
sigresult
write.table(sigresult, "up.psRNATarget.enriched.GO.genes_cheng.txt", sep="\t", quote=F, row.names=F)

######################
### barplot GO enrichment ter of up.sRNA.targeted.genes
######################
pdf("up.psRNATarget.sig.GO.enrich_cheng.pdf", width=6, height=6)

plotdata <- sigresult$Cat_in_up.sRNA.target_gene_num
names(plotdata) <- sigresult$Term
plotdata <- plotdata[!is.na(names(plotdata))]
par(mar=c(4,18,1,1))
barplot(plotdata, horiz=T, las=2, col="cornflowerblue", axes=F,cex.names=0.6, 
        xlab="Number of significant genes", cex.lab=1)
abline(v=0)
#axis(side=1, at = seq(0, max(plotdata), by=round(max(plotdata)/6, 0)))
axis(side=1, at =c((0:10)*10))

dev.off()


#####################
### GO analysis shock
#####################
go <- read.delim("./inputs/ZmB73_v3.gene2go.txt")
gf <- data.frame(GeneID=unique(go$genes))
gf$shock.sRNA.target <- "no"
gf$shock.sRNA.target[gf$GeneID %in% shock.sRNA2genes] <- "shock.sRNA.target"
table(gf$shock.sRNA.target) 
go.result <- enrich.test(gene2feature=gf, gene2cat=go,
                         gene2bias=NULL, test.feature="shock.sRNA.target",
                         pval.cutoff=0.01, enrich.cutoff=2,
                         method="Sampling", sampling=5000, outsave=T,
                         outfile="./enrichments/shock.sRNA.targeted.genes.GO.enrichment.txt")

sigresult <- go.result$shock.sRNA.target[go.result$shock.sRNA.target$pvals<0.01 & go.result$shock.sRNA.target$enrich>2, ]
nrow(sigresult)
sigresult$Term <- sapply(as.character(sigresult$Cat), term.fun)
sigresult
write.table(sigresult, "shock.psRNATarget.enriched.GO.genes_cheng.txt", sep="\t", quote=F, row.names=F)

######################
### barplot GO enrichment ter of shock.sRNA.targeted.genes
######################
pdf("shock.psRNATarget.sig.GO.enrich_cheng.pdf", width=6, height=6)

plotdata <- sigresult$Cat_in_shock.sRNA.target_gene_num
names(plotdata) <- sigresult$Term
plotdata <- plotdata[!is.na(names(plotdata))]
par(mar=c(4,18,1,1))
barplot(plotdata, horiz=T, las=2, col="cornflowerblue", axes=F,cex.names=0.6, 
        xlab="Number of significant genes", cex.lab=1)
abline(v=0)
#axis(side=1, at = seq(0, max(plotdata), by=round(max(plotdata)/6, 0)))
axis(side=1, at =c((0:10)*10))

dev.off()


