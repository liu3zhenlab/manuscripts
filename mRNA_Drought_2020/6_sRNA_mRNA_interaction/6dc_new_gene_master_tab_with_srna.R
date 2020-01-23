setwd("~/Documents/Cheng_project/121917_B73_drought/Analysis/5_new_sRNA_target_mRNA_predict/5c_psRNATarget_sRNA_mRNA_interaction/input/")
oldpair <- read.delim("sRNA_mRNA_minus_match.pairs.txt")
pspair <- read.delim("psRNATarge_sRNA_merge_miRNA_pairs.txt")
mirna <- read.delim("miRNA_predicted_target.txt")
head(oldpair)
olevel <- levels(oldpair$GeneID)[-1]
plevel <- levels(pspair$GeneID)[-1]
#i <- 1801
oldout <- c()
for (i in 1:length(olevel))
{
  a <- as.character(oldpair$sRNA[oldpair$GeneID %in% olevel[i]])
  c <- length(unique(a))
  if (c > 30)
  {
    b <- ">30"
  } else {
    b <- paste(unique(a),collapse = ";",sep="")
  }
  d <- c(olevel[i],b,c)
  oldout <- rbind(oldout,d)
}
colnames(oldout) <- c("GeneID","Targed_sRNA","Num of sRNAs")
#write.table(oldout,file="oldpair_gene_master_tab.txt",sep="\t",quote=FALSE,row.names = FALSE)
head(oldout)


psout <- c()
for (i in 1:length(plevel))
{
  a <- as.character(pspair$sRNA[pspair$GeneID %in% plevel[i]])
  c <- length(unique(a))
  if (c > 30)
  {
    b <- ">30"
  } else {
    b <- paste(unique(a),collapse = ";",sep="")
  }
  d <- c(plevel[i],b,c)
  psout <- rbind(psout,d)
}
colnames(psout) <- c("GeneID","psRNATarged_sRNA","Num of psRNATarget")
#write.table(psout,file="pspair_gene_master_tab.txt",sep="\t",quote=FALSE,row.names = FALSE)
head(psout)

mout <- merge(oldout,psout,by="GeneID",all=TRUE)

head(mirna)
gene <- gsub("_T[0-9][0-9]","",mirna$Transcript)
gene <- gsub("_FGT","_FG",gene)
mirna <- cbind(mirna,gene)
milevel <- levels(mirna$gene)
miout <- c()
for (i in 1:length(milevel))
{
  a <- as.character(mirna$miRNA[mirna$gene %in% milevel[i]])
  b <- paste(unique(a),collapse = ";",sep="")
  c <- length(unique(a))
  d <- c(milevel[i],b,c)
  miout <- rbind(miout,d)
}
colnames(miout) <- c("GeneID","known_miRNA","Num of miRNA")
fmout <- merge(mout,miout,by="GeneID",all=TRUE)
#write.table(fmout,file="new_gene_master_tab_with_sRNA.txt",sep="\t",quote=FALSE,row.names = FALSE)

mmast <- read.delim("mRNA_drought_timeseries_output1_cheng.txt")
head(mmast)
out <- merge(mmast,fmout,by="GeneID",all.x=TRUE)
write.table(out,file="new_gene_master_tab_with_sRNA.txt",sep="\t",quote=FALSE,row.names = FALSE)