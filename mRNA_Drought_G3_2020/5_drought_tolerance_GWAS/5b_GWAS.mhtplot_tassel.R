######################################################
###############The variables will be used##################
######################################################
chr.length <- c(301354135,237068873,232140174,241473504,217872852,
                169174353,176764762,175793759,156750706,150189435)
#use a 'for' cycle to produce several useful vectors.
for(i.chr in 1:10){
    if(i.chr==1){
        pos.elment1 <- sum(chr.length[1:i.chr])
        pos.elment2 <- pos.elment1/2
        pos <- pos.elment1
        pos.chr <- pos.elment2
    }else{
        pos.elment2 <- pos.elment1+chr.length[i.chr]/2
        pos.elment1 <- sum(chr.length[1:i.chr])
        pos <- c(pos,pos.elment1)
        pos.chr <- c(pos.chr,pos.elment2)}
}
file.no <- list.files()
title.names <- paste(file.no,"_mhtplot",sep="")
data.file <- file.no
i <- 2
######################################################
#############Using 'for' cycles to draw mahtplot#############
######################################################
for(i in 1:length(file.no)){
   gwas.data <- read.table(data.file[i],sep="\t",header=TRUE)
   p <- gwas.data$p
   ps <- sort(p)
   order <- c(1:length(ps))/length(ps)*0.05
   bhcutoff <- ps[which((ps-order) > 0)[1]]
   mht.data.pro <- gwas.data[,c(3,4,7)]
   dimnames(mht.data.pro)[[2]] <- c("chr","pos","p")
#use a second 'for' cycle to modify data structure.
   for(j in 2:10){
       if(j==2){
           chr.data <- mht.data.pro[mht.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           mht.data <- rbind(mht.data.pro[mht.data.pro[,1]==1,],chr.data)
       }else{
           chr.data <- mht.data.pro[mht.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           mht.data <- rbind(mht.data,chr.data)
            }
   }
#set some variables for mhtplot.
   pngfile.names <- paste(file.no[i],"_GWAS.png",sep="")
   title.name <- title.names[i]
   threshold <- -log10(bhcutoff)
   threshold_l <- 4
   label.name<-c(paste("chr.",1:10,sep=""))
   color.array<-rep(c("red","blue"),5)
  
#open a png file for plotting.
   png(pngfile.names,res=200,height=1500,width=2500)
   y <- -log10(mht.data[,3])
   x <- mht.data[,2]
   plot(x,y,type="p",cex=0,xlab="Maize chromosomes region.",ylab="-log10(p-value)",xlim=c(0,pos[10]),
        ylim=c(0,10),xaxs="i",yaxs="i",xaxt="n",family="serif",font=2)       
#use a 'for' cycle to draw points with differnt colors.
   for(k in 1:10){
       x.zc <- mht.data[mht.data[,1]==k,2]
       y.zc <- -log10(mht.data[mht.data[,1]==k,3])
       points(x.zc,y.zc,col=color.array[k],pch=19,cex=0.5)
   }      
#modify and embellish the mhtplot.
   axis(1,at=pos.chr,labels=paste("chr",1:10,sep=""),family="serif",font=2)
   lines(x=c(0,pos[10]),y=c(threshold,threshold),lty=2,type="l",col="black",lwd=2)
   lines(x=c(0,pos[10]),y=c(threshold_l,threshold_l),type="l",col="black",lwd=2)
   lines(x=c(0,pos[10]),y=c(0,0),type="l",col="black",lwd=1)
   title(title.name,cex.main=2)
   par(family="serif",font.lab=2,font.axis=2)
   dev.off()
}
