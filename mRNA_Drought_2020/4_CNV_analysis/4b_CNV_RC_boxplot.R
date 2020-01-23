x <- read.table("cor_fdr0.05_mRNA_eml_new.txt",sep="\t",header=TRUE)
head(x)
drg <- x[which(x$interaction.padj < 0.05),]
nodrg <- x[which(x$interaction.padj > 0.05),]
e <- x[which(x$earlySig=="yes"),]
m <- x[which(x$middleSig=="yes"),]
l <- x[which(x$lateSig=="yes"),]
drgp <- drg$cor_RC_kmernf
nodrgp <- nodrg$cor_RC_kmernf
ep <- e$cor_RC_kmernf
mp <- m$cor_RC_kmernf
lp <- l$cor_RC_kmernf
boxplot(drgp,nodrgp,ep,mp,lp,ylim=c(0.87,1.02),col=c("brown","gray","red","blue","green"),outline = FALSE,range=0.5,names=c("DRG","noDRG","Early","Middle","Late"),ylab="R2",boxwex=0.7)

t.test(drgp,nodrgp)
t.test(ep,lp)
t.test(ep,mp)
t.test(mp,lp)

?boxplot()