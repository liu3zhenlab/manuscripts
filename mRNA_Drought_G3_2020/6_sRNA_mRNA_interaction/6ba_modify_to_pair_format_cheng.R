x <- read.delim("cutoff_2_top200.txt")
head(x)
srna <- x$miRNA_Acc.
gene <- gsub("_T[0-9][0-9]","",x$Target_Acc.)
gene <- gsub("_FGT","_FG",gene)
trans <- x$Target_Acc.
pos <- x$Target_start
match <- x$Expectation
stand <- "-"
ed <- "NA"
y <- cbind(as.character(srna),as.character(gene),as.character(trans),pos,match,stand,ed)
head(y)
colnames(y) <-  c("sRNA", "GeneID", "Transcript", "TranscriptPos", "Match", "Strand", "Editdistance")
t <- table(x$miRNA_Acc.)
flag <- c(1:length(x$miRNA_Acc.))
l <- names(t)[t > 30]
#i <- l[1]
all <- c()
for (i in l)
{
  na <- flag[x$miRNA_Acc. %in% i]
  add <- y[na[1],]
  add[2:7] <- c(">30",">30","NA","NA","NA","NA")
  all <- rbind(all,add)
}
rm <- flag[x$miRNA_Acc. %in% all[,1]]
z <- y[-rm,]
z <- rbind(z,all)
z <- z[order(z[,1]),]
write.table(z,file="psRNATarge_cuf_2_sRNA_pairs.txt",sep="\t",quote=FALSE,row.names = FALSE)
