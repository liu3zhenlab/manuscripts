smaster <- read.table("all.gt72.sRNA.sub.mt3_cheng.txt",sep="\t",header=TRUE)
ps <- read.table("psRNATarge_sRNA_merge_miRNA_pairs.txt",sep="\t",header=TRUE)
srna_ps <- unique(ps$sRNA)
#i <- 10
unique <- c()
for (i in 1:length(srna_ps))
{
  sample <- ps[ps$sRNA %in% srna_ps[i],]
  psRNATarget_GeneID <- paste(unique(sample$GeneID),collapse = ";")
  row <- c(as.character(srna_ps[i]),as.character(psRNATarget_GeneID))
  unique <- rbind(unique,row)
}
colnames(unique) <- c("sRNA","psRNATarget_GeneID")
write.table(unique,file="unique_sRNA_origin.txt",sep="\t",quote=FALSE,row.names = FALSE)
unique_v2 <- read.delim("unique_sRNA_origin.txt")
mirna <- ps[which(ps$Editdistance=="miRNA_known"),]
mirna_unique <- unique(mirna$sRNA)
psRNATarget_mirna <- rep("NA",length(unique$sRNA))
psRNATarget_mirna[unique_v2$sRNA %in% mirna_unique] <- "miRNA_known"
unique_v2 <- cbind(unique_v2,psRNATarget_mirna)
head(unique_v2)
write.table(unique_v2,file="unique_sRNA.txt",sep="\t",quote = FALSE,row.names = FALSE)
smaster_new <- merge(smaster,unique_v2,by="sRNA",all.x=TRUE)
head(smaster_new)
write.table(smaster_new,file="all.gt72.sRNA.sub.spRNATarget_cheng.txt",sep="\t",quote=FALSE,row.names = FALSE)
