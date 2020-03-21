#############################################################################
### mRNA-CNV
### Sanzhen liu
### Modified by Cheng
### 01/08/2018
#############################################################################
proj <- "B73Drought"
proc <- "CNV"
fig.d <- paste0(wd, "figures")
table.d <- paste0(wd, "tables")

pid <- paste0(proj, "_PROC", proc)  # procedure ID
logfile <- paste0(pid, ".log")
cat(pid, "\n", file = logfile)
outlog <- function (text) cat(text, "\n", file = logfile, append = T)  # define outlog function
setwd(wd)

#############################################################################
### modules
#############################################################################
source("output.standard.R")
#png.sl("test", outpath = fig.d, prefix = pid)
#write.table.sl(test, "test", outpath = table.d, prefix = pid)

#############################################################################
### data
#############################################################################
# drought response gene:
dr1 <- read.delim("DRG_master_table_v2.txt")
head(dr1)
colnames(dr1)[1] <- "GeneID"

# cnv genes from Ruth 2010 Genome Research
cnv <- read.delim("Swanson-Wagner2010_zea_CNV.tableS2.txt")
head(cnv)
nrow(cnv)
cnv2 <- cnv[, c(1, 4)]
colnames(cnv2) <- c("GeneID", "CNV")
table(cnv2$CNV)

# merge by drought genes and CNV table
dr2 <- merge(dr1, cnv2, by = "GeneID", all.x  = T) 
dr2$CNV <- as.character(dr2$CNV)
dr2$CNV[is.na(dr2$CNV)] <- "noCNV"
head(dr2)

# define DR response genes and CNV genes
dr2$CNV2 <- "CNV"
dr2$CNV2[dr2$CNV == "noCNV"] <- "noCNV"

dr2$DRgenes <- "no"
dr2$DRgenes[!is.na(dr2$Cluster)] <- "yes"

head(dr2)

# Up/Down Regulation vs. CNV
table(dr2$CNV, dr2$RegulationOnDS)
table(dr2$RegulationOnDS)
table(dr2$CNV, dr2$RegulationOnDS)[2, ] / table(dr2$RegulationOnDS)

# Cluster vs. CNV
table(dr2$CNV, dr2$Cluster)
table(dr2$Cluster)
table(dr2$CNV, dr2$Cluster)[2, ]/table(dr2$Cluster)


### early, middle, late vs. CNV
table(dr2$CNV, dr2$DE)
chisq.test(matrix(c(518, 970, 8108, 16554), nrow = 2, byrow = T))
table(dr2$CNV, dr2$earlySig)
chisq.test(matrix(c(19, 1469, 174, 24483), nrow = 2, byrow = T))
table(dr2$CNV, dr2$middleSig)
chisq.test(matrix(c(175, 1279, 2135, 23292), nrow = 2, byrow = T))
table(dr2$CNV, dr2$lateSig)
chisq.test(matrix(c(156, 1277, 2345, 21435), nrow = 2, byrow = T))

table(dr2$CNV, dr2$earlySig)[2, ]/table(dr2$earlySig)
table(dr2$CNV, dr2$middleSig)[2, ]/table(dr2$middleSig)
table(dr2$CNV, dr2$lateSig)[2, ]/table(dr2$lateSig)


table(dr2$CNV, paste(dr2$earlySig, dr2$middleSig, dr2$lateSig))
table(paste(dr2$earlySig, dr2$middleSig, dr2$lateSig))
sig3 <- paste(dr2$earlySig, dr2$middleSig, dr2$lateSig)
table(dr2$CNV, sig3)[2, ]/table(sig3)

table(as.character(is.na(dr2$CNV)), sig3)[1, ]/table(sig3)
dr2[sig3 == "yes no no", ]

table(dr2$CNV, paste(dr2$earlySig, dr2$middleSig, dr2$lateSig))
table(paste(dr2$earlySig, dr2$middleSig, dr2$lateSig))
sig4 <- paste(dr2$earlySig, dr2$middleSig, dr2$lateSig, dr2$RegulationOnDS)
table(dr2$CNV, sig4)[2, ]/table(sig4)


enrich.test <- function(data2x2, contrast.col = T, enrich.factor = 1) {
  chi_pval <- chisq.test(data2x2)[[3]]
  if (contrast.col) {
    factor.counts <- data2x2[enrich.factor, ]
    all.counts <- colSums(data2x2)
  } else {
    factor.counts <- data2x2[, enrich.factor]
    all.counts <- rowSums(data2x2)
  }
  freq <- round(factor.counts / all.counts, 4)
  out <- c(chi_pval, paste0(freq, " (", factor.counts, "/", all.counts, ")" ))
  names(out) <- c("Chisq.pval", names(all.counts))
  out
}

#chisq.test(table(dr2$CNV2, dr2$DRgenes))
enrich.test(table(dr2$CNV2, dr2$DRgenes))
enrich.test(table(dr2$CNV2, dr2$earlySig))
enrich.test(table(dr2$CNV2, dr2$middleSig))
enrich.test(table(dr2$CNV2, dr2$lateSig))

dr.cnv <- dr2[dr2$DRgenes == "yes" & dr2$CNV2 == "CNV", c("GeneID", "interaction.padj", "Symbol", "Gene.Name", "gene_description", "a.th_ortholog", "rice_ortholog")]
nrow(dr.cnv)
write.table(dr.cnv,file="dr.cnv.anno.list_cheng.txt",sep="\t",quote=FALSE,row.names = FALSE)

head(dr2)
earlysig.cnv <- dr2[dr2$earlySig == "yes" & dr2$CNV2 == "CNV", c("GeneID", "interaction.padj", "Symbol", "Gene.Name", "gene_description", "a.th_ortholog", "rice_ortholog")]
nrow(earlysig.cnv)
write.table(earlysig.cnv,file="earlysig.cnv.anno.list_cheng.txt",sep="\t",quote=FALSE,row.names = FALSE)

middlesig.cnv <- dr2[dr2$middleSig == "yes" & dr2$CNV2 == "CNV", c("GeneID", "interaction.padj", "Symbol", "Gene.Name", "gene_description", "a.th_ortholog", "rice_ortholog")]
nrow(middlesig.cnv)
write.table(middlesig.cnv,file="middlesig.cnv.anno.list_cheng.txt",sep="\t",quote=FALSE,row.names = FALSE)

latesig.cnv <- dr2[dr2$lateSig == "yes" & dr2$CNV2 == "CNV", c("GeneID", "interaction.padj", "Symbol", "Gene.Name", "gene_description", "a.th_ortholog", "rice_ortholog")]
nrow(latesig.cnv)
write.table(latesig.cnv,file="latesig.cnv.anno.list_cheng.txt",sep="\t",quote=FALSE,row.names = FALSE)
