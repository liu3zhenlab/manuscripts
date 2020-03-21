################################################################################
### Identification of the association between sRNA and mRNA
### Sanzhen Liu
### Modified by Cheng
### 01/07/2018
################################################################################
source("num.fuzzy.group.R")
######################################################
### correlation function
######################################################
cor.cal <- function(data, scols, mcols) {
  sRNA.rc <- as.numeric(as.character(data[scols]))
  gene.rc <- as.numeric(as.character(data[mcols]))
  
  sRNA.logcounts <- log(sRNA.rc+1)
  gene.logcounts <- log(gene.rc+1)
  
  mscor <- cor(sRNA.logcounts, gene.logcounts)
  return(mscor)
}
######################################################


######################################################
### correlation function
######################################################
mscor.test <- function(data, scols, mcols) {
  sRNA.rc <- as.numeric(as.character(data[scols]))
  gene.rc <- as.numeric(as.character(data[mcols]))
  
  sRNA.logcounts <- log(sRNA.rc+1)
  gene.logcounts <- log(gene.rc+1)
  
  mscor <- cor.test(sRNA.logcounts, gene.logcounts, alternative="less")[[3]]
  return(mscor)
}
######################################################


######################################################
### to test (beta < 0)
######################################################
linear.test <- function(data, scols, mcols) {
  sRNA.rc <- as.numeric(as.character(data[scols]))
  gene.rc <- as.numeric(as.character(data[mcols]))
  sRNA.logcounts <- log(sRNA.rc + 1)
  gene.logcounts <- log(gene.rc + 1)
  lm(gene.logcounts ~ sRNA.logcounts)
  mscor <- cor.test(sRNA.logcounts, gene.logcounts, alternative="less")[[3]]
  return(mscor)
}
######################################################


#### mRNA data:
mrd <- read.delim("mRNA.normalized.counts_cheng.txt")
colnames(mrd)[1] <- "GeneID"
colnames(mrd)
head(mrd)
dim(mrd)

### sRNA data:
srd <- read.delim("sRNA.sub.normalized.counts.txt")
nrow(srd)

### matched gene and sRNA
#mslink <- read.delim("../../sRNA_analysis/analysis2_sRNA_mRNA_pairs_otherspecies/conserved.sRNA_mRNA_minus_match.pairs.txt")
mslink <- read.delim("psRNATarge_sRNA_merge_miRNA_pairs.txt", stringsAsFactors=F)

dim(mslink)
head(mslink)

sRNAs <- read.delim("all.gt72.sRNA.sub.spRNATarget_cheng.txt")
dim(sRNAs)

mslink <- mslink[mslink$sRNA %in% sRNAs$sRNA, ]

sum(mslink$GeneID == ">30")  ### sRNA targeting many (>30) genes

tail(mslink)
# remove sRNAs with high number of hits
mslink2 <- mslink[mslink$GeneID != ">30", ]  ### high number of hits
# remove duplicated pairs of sRNA-mRNAs
sg <- paste(mslink2[, 1], mslink2[, 2])
mslink3 <- mslink2[!duplicated(sg), ]

####################################################################
### characterization of sRNA-mRNA pairs
####################################################################
nrow(mslink3)


####################################################################

### determine the neighboring positions:
pos.groups <- num.fuzzy.group(category=mslink3$GeneID, num=mslink3$TranscriptPos, distance = 10)
pos.groups <- pos.groups[!duplicated(paste(pos.groups$category, pos.groups$data)), ]
nrow(pos.groups)
head(pos.groups,10)

mslink3 <- merge(mslink3, pos.groups, by.x=c("GeneID", "TranscriptPos"), by.y=c("category", "data"), all.x=T)
head(mslink3)
dim(mslink3)

###################################

alldata <- merge(mslink3, srd, by="sRNA")
nrow(alldata)
sRNA.cols <- grep("^D[SRW]", colnames(alldata))
colnames(alldata)[sRNA.cols]
alldata$rcsum <- rowSums(alldata[, sRNA.cols])
maxcounts <- aggregate(rcsum ~ GeneID + group, data=alldata, max)
head(maxcounts)
alldata <- merge(alldata, maxcounts, by=colnames(maxcounts))
alldata <- merge(alldata, mrd, by="GeneID") 
nrow(alldata)
head(alldata)
length(unique(alldata$sRNA))
length(unique(alldata$GeneID))
sRNA.cols <- grep(".x", colnames(alldata))
colnames(alldata)[sRNA.cols]
mRNA.cols <- grep(".y", colnames(alldata))
colnames(alldata)[mRNA.cols]

###
### correlation and correlation test
###
allcor <- apply(alldata, 1, cor.cal, scols=sRNA.cols, mcols=mRNA.cols)
hist(allcor, nclass=20, xlim=c(-1, 1))

allcor.pvals <- apply(alldata, 1, mscor.test, scols=sRNA.cols, mcols=mRNA.cols)
hist(allcor.pvals)
allcor.qvals <- p.adjust(allcor.pvals, method="BH")
sum(allcor.qvals < 0.05)
cor.set <- alldata[,1:9]
cor.set$Cor <- allcor
head(alldata)
write.table(cor.set,file="all.srna.mrna.cor.txt",sep="\t",quote=FALSE,row.names = FALSE)
neg.cor.set <- alldata[allcor.qvals < 0.05, 1:9]
neg.cor.set$Cor <- allcor[allcor.qvals < 0.05]
neg.cor.set
nrow(neg.cor.set)

sig.mRNAs <- read.delim("../../../2_mRNA_analysis/2a_mrna_interaction/results/sig.interaction.mRNAs.coefficients.output.txt")
sum(neg.cor.set$GeneID %in% sig.mRNAs$mRNA)
neg.cor.set[neg.cor.set$GeneID %in% sig.mRNAs$mRNA, ]

###
### combined with other information
###
### gene functional annotation:
annot <- read.delim("gene_v3_annotation_pang_nodup_brief.txt")
head(annot)
annot[annot$gene_name %in% intersect(sig.mRNAs$mRNA, neg.cor.set$GeneID), ]
annot[annot$gene_name %in% neg.cor.set$GeneID, ]

TF <- read.delim("maize.TF.family.txt")
TF <- TF[, 1:2]
tfset <- merge(TF, neg.cor.set, by="GeneID", all.y=T)
head(tfset)
tfset <- tfset[, c("GeneID", "TF", "sRNA", "Transcript", "TranscriptPos", "Match", "Editdistance", "Cor")]
nrow(tfset)
colnames(tfset) <- c("GeneID", "TF", "sRNA", "Transcript", "TranscriptPos", "Match", "Editdistance", "gene_sRNA_correlation")
tfset$gene_sRNA_correlation <- round(tfset$gene_sRNA_correlation, 3)
nrow(tfset)

### sRNA master table
smt <- read.delim("all.gt72.sRNA.sub.spRNATarget_cheng.txt")
smt2 <- smt[, c("sRNA", "miRNA_perfect_hits", "Group", "TimeSeriesGroup", "psRNATarget_GeneID", "psRNATarget_mirna")]
head(smt2)
smt2[smt2$sRNA %in% tfset$sRNA, ]

tfset <- merge(tfset, smt2, by="sRNA", all.x=T)
head(tfset)
colnames(tfset)
nrow(tfset)

table(tfset$miRNA_perfect_hits)
table(tfset$Group)

################################################
### sig sRNA characterization
################################################
tfset.nr <- tfset[!duplicated(tfset$sRNA), ]
nrow(tfset.nr)
sum(grepl("MIR", tfset.nr$miRNA_perfect_hits))
table(tfset.nr$Group)
table(smt2$Group)

### chisq.test:
miRNAs <- c(11, 16, 2861-11, sum(table(smt2$Group[smt2$Group != "others"])) - 2861 - 16)
chisq.test(as.matrix(miRNAs, nrow=2, byrow=T))

snoRNAs <- c(7, 20, 5630 - 7, sum(table(smt2$Group[smt2$Group != "others"])) - 5630 - 20)
chisq.test(as.matrix(snoRNAs, nrow=2, byrow=T))

tfset.nr$TimeSeriesGroup[grep("snoRNA", tfset.nr$Group)]
tfset.nr$TimeSeriesGroup[grep("miRNA", tfset.nr$Group)]
tfset.nr$TimeSeriesGroup[grep("tRNA", tfset.nr$Group)]

################################################


### master table
mt <- read.delim("fd0.05_mRNA_master_table.txt", stringsAsFactors=F)
tfset2 <- merge(tfset, mt[, c("GeneID", "chr", "start", "end", "strand", "Symbol", "Gene.Name", "gene_description", "a.th_ortholog", "rice_ortholog", "interaction.padj", "earlySig", "middleSig", "lateSig", "RegulationOnDS","CNV_SW","CNV_Lai")], all.x=T)
tfset2
nrow(tfset2)

################################################
### mRNA characterization
################################################
table(tfset2$gene_description)
colnames(tfset2)

tfset2.defined <- tfset2[!grepl("^Uncharacterized protein", ignore.case = FALSE, tfset2$gene_description), ]
table(tfset2.defined$Group, tfset2.defined$gene_description)
table(tfset2.defined$gene_description, tfset2.defined$Group)
table(tfset2.defined$gene_description)

tfnum <- table(as.character(tfset2$TF))
sum(tfnum)

table(as.character(tfset2$gene_description))

dim(mt)

###
### GO term analysis
###
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
go <- read.delim("ZmB73_v3.gene2go.txt")
geneset <- mt[1]
geneset$Target <- "no"
geneset$Target[geneset$GeneID %in% unique(tfset2$GeneID)] <- "yes"

### go enrichment
go.result <- enrich.test(gene2feature=geneset, gene2cat=go,
                         gene2bias=NULL, test.feature="yes",
                         pval.cutoff=0.01, enrich.cutoff=2,
                         method="Sampling", sampling=5000, outsave=T,
                         outfile="sRNA.negative.correlation.genes.GO.enrichment.txt")

sigresult <- go.result$yes[go.result$yes$pvals<0.01 & go.result$yes$enrich>1.5, ]
nrow(sigresult)
sigresult$Term <- sapply(as.character(sigresult$Cat), term.fun)
sigresult
nrow(sigresult)

write.table(sigresult, "allGO.sRNA.negative.correlation.genes_psRNATarget.txt", sep="\t", quote=F, row.names=F)

###
### TF
###
tfset2.nr <- tfset2[!duplicated(tfset2$GeneID), ]
table(as.character(tfset2.nr$TF))

table(as.character(tfset2$TF), tfset2$Group)
################################################


################################################
### conserved sRNAs
################################################
conserved <- read.delim("conserved.sRNA_mRNA_minus_match.pairs.txt")
head(conserved)

sig.gs <- paste(tfset2$GeneID, tfset2$sRNA)
conserved.gs <- paste(conserved$GeneID, conserved$sRNA)

conserved.sig <- tfset2[sig.gs %in% conserved.gs, ]

nrow(conserved.sig)
length(unique(conserved.sig$sRNA))
length(unique(conserved.sig$GeneID))

tfset2$Conserved <- "no"
tfset2$Conserved[sig.gs %in% conserved.gs] <- "yes"
################################################

###
### output
###
write.table(tfset2, "fdr0.05.all.gene_sRNA.sig.negative.cor.psRNATarget.txt", quote=F, row.names=F, sep="\t")



### others
known <- read.delim("PMTED.gene_sRNA.sig.negative.cor.txt")

sum(paste(known$GeneID, known$miRNA_seq) %in% paste(tfset2$GeneID, tfset2$sRNA))

