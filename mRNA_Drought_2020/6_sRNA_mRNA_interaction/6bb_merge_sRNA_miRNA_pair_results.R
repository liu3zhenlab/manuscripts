srna <- read.delim("psRNATarge_cuf_1_sRNA_pairs.txt")
mirna <- read.delim("psRNATarge_cuf_2_miRNA_pairs.txt")
a <- cbind(paste(srna$sRNA,srna$Transcript,sep="_"),c(1:length(srna$sRNA)))
colnames(a) <- c("srna","flag_srna")
b <- cbind(paste(mirna$sRNA,mirna$Transcript,sep="_"),c(1:length(mirna$sRNA)))
colnames(b) <- c("mirna","flag_mirna")
d <- a[a[,1] %in% b[,1],]
merge <- rbind(srna[-as.numeric(d[,2]),],mirna)
head(merge)
merge_sort <- merge[order(merge$sRNA),]
write.table(merge_sort,file="psRNATarge_sRNA_merge_miRNA_pairs.txt",sep="\t",quote=FALSE,row.names = FALSE)