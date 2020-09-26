setwd("/bulk/liu3zhen/research/projects/kGWAS/ft50k/2-CON")
options(stringsAsFactors=F)

##############################################################################
# mapping information of a cluster
##############################################################################
clustinfo <- function(alnfile, min_major_chr=0.8, max_region_range=10000000) {
  clust <- gsub(".*\\/", "", alnfile)
  clust <- gsub("3o\\-clust_", "", clust)
  clust <- gsub("\\..*", "", clust)
  
  finfo = file.info(alnfile)
  if (finfo$size == 0) {
    out <- c(clust, rep(NA, 5))
  } else {
    
    sam <- read.delim(alnfile, header=F)
    ctgkmers.file <- paste0("1o-clust_", clust, ".fasta.cap.ace.contig.kmers")
    
    sam <- sam[, c(1,3,4)]
    colnames(sam) <- c("ctg", "chr", "pos")
    
    ### mapping information of the cluster:
    clust.chr <- NA
    clust.start <- NA
    clust.end <- NA
    if (nrow(sam) > 0) {
      # major chr
      chr.counts <- table(sam$chr)
      major_chr <- names(which.max(chr.counts))
      prop_major_chr <- max(chr.counts) / nrow(sam)
      pos <- sam[sam$chr == major_chr, "pos"]
      min_pos <- min(pos)
      max_pos <- max(pos)
      if ((max_pos - min_pos + 1) > max_region_range & nrow(sam) >= 5) {
        min_pos <- quantile(pos, 0.1)
        max_pos <- quantile(pos, 0.9)
      }
      
      if (prop_major_chr >= min_major_chr & (max_pos - min_pos + 1) <= max_region_range) {
        # high prop of major chr and the range is small enought
        clust.chr <- major_chr
        clust.start <- min_pos
        clust.end <- max_pos
      }
    }
    
    ################################################
    ctgkmers <- read.delim(ctgkmers.file, header=F)
    colnames(ctgkmers) <- c("ctg", "bp", "kmer")
    ctgkmers$pvalue <- gsub(".*_", "", ctgkmers$kmer)
    ctgkmers$pvalue <- as.numeric(ctgkmers$pvalue)
    ctgpval <- tapply(ctgkmers$pvalue, ctgkmers$ctg, mean)
    clust.pval <- mean(ctgpval)
    clust.pval <- format(clust.pval, scientific=T, digits=2)
    clust.ctglens <- sum(ctgkmers[!duplicated(ctgkmers$ctg), "bp"])
    
    out <- c(clust, clust.ctglens, clust.pval, clust.chr, clust.start, clust.end)
  }
  out
}

##############################################################################
### plotting of mapping result on chromosomes
##############################################################################
pval.oneplot <- function (infile, yrange = NULL, xlab.text="chromosome",
                          mar.set =c(4.5, 4.5, 0.5, 0.5), mgp.set = c(2.5, 0.8, 0),
                          cols = c("grey90", "grey80"), pvalue.cutoff = NULL,
                          ylab.text = "-log10(p-values)", na.rm = T, pval_nonzero = T,
                          main.text="", axisline.width=1.5, chr.set, plot.cex.lab = 1.3,
                          chr.size, order.by.chrsize=F, label.rm=NULL, cexaxis=1.3,
                          cex.line.width = 0.5, saveplot=FALSE, plot.width=3, plot.heigh=2.5, 
                          plot.path=".", plot.filename="default.png", ...) {
  
  plot.path <- gsub("/$", "", plot.path)
  input <- read.delim(infile, stringsAsFactors=F)
  colnames(input) <- c("Cluster", "length", "Pval", "Chr", "Start", "End")
  input$Pos <- round((input$Start + input$End) / 2)
  input$Chr <- as.character(input$Chr)
  input <- input[input$Chr %in% chr.set, ]
  
  if (na.rm) {
    input <- input[!is.na(input$Pval), ]
  }
  
  if (pval_nonzero) {
    input[input$Pval == 0, "Pval"] = min(input[input$Pval != 0, "Pval"])
  }
  
  input$neglogP <- -log10(input$Pval)
  
  print(input)
  
  # chromosome length:
  colnames(chr.size) <- c("Chr", "Size")
  chr.size$Chr <- as.character(chr.size$Chr)
  chr.size <- chr.size[chr.size$Chr %in% chr.set, ]
  
  if (order.by.chrsize) {
    chr.size <- chr.size[order(chr.size$Size, decreasing=T), ]
  }
  
  # plotting
  accum <- 0
  all.col <- NULL
  all.chr <- NULL
  centers <- NULL
  gap <- sum(chr.size$Size)/80
  if (is.null(yrange)) {
    ymax <- max(input$neglogP)
    yrange <- c(0, ymax)
  } else {
    ymax <- max(yrange)
  }
  
  if (saveplot) {
    pdf(paste(plot.path, plot.filename, sep="/"), width=plot.width, heigh=plot.heigh, useDingbats=F)
  }
  
  ### plot
  par(mar = mar.set, mgp = mgp.set)
  plot(NULL, NULL, ylim = yrange,
       xlim = c(0, gap * nrow(chr.size) + sum(chr.size$Size)),
       xaxt = "n", yaxt = "n",
       xlab = xlab.text, ylab = ylab.text,
       main = main.text, cex.lab = plot.cex.lab,
       bty = "n", cex.axis = cexaxis, lwd = cex.line.width)
  
  axis(side = 2, cex.axis = cexaxis, lwd = cex.line.width)
  box(lwd = cex.line.width)
  
  
  all.accum <- NULL
  if (length(cols) < nrow(chr.size)) {
    cols <- rep(cols, ceiling(nrow(chr.size) / length(cols)))
  }
  
  for (i in 1:(nrow(chr.size))) {
    all.accum <- c(all.accum, accum)
    pre.accum <- accum
    chr <- chr.size[i, "Chr"]
    len <- chr.size[i, "Size"]
    plot.col <- input[input$Chr==chr, "Cluster"]
    plot.size <- input[input$Chr==chr, "length"]
    pos <- input[input$Chr==chr, "Pos"]
    prob <- input[input$Chr==chr, "neglogP"]
    polygon(x=c(accum, accum, accum+len, accum+len), border=NA,
            y=c(yrange[1], yrange[2], yrange[2], yrange[1]), col=cols[i %% 2 + 1])
    points(accum+pos, prob, lwd=cex.line.width, pch=21, col="black", bg=plot.col, ...)
    accum <- accum + len + gap
    center.point <- (pre.accum + accum - gap)/2
    all.col <- c(all.col, plot.col)
    all.chr <- c(all.chr, chr)
    centers <- c(centers, center.point)
  }
  
  if (!is.null(label.rm)) {
    for (each.label.rm in label.rm) {
      all.chr <- gsub(each.label.rm, "", all.chr)
    }
  }
  
  # axis
  axis(side = 1, at = centers, labels=all.chr, tick=F, cex.axis=cexaxis, lwd = cex.line.width)
  
  if (!is.null(pvalue.cutoff)) {
    abline(h = -log10(pvalue.cutoff), col = "grey", lwd = cex.line.width, lty = 2)
  }
  
  if (saveplot) dev.off()
}

##############################################################################
### A188
##############################################################################
samfiles <- dir(pattern="A188uniq.sam")
allcinfo <- NULL
for (samfile in samfiles) {
  cinfo <- clustinfo(alnfile=samfile, min_major_chr=0.6, max_region_range=10000000)
  if (is.null(allcinfo)) {
    allcinfo <- cinfo
  } else {
    allcinfo <- rbind(allcinfo, cinfo)
  }
}

allcinfo <- data.frame(allcinfo)
colnames(allcinfo) <- c("Cluster", "ctgBp", "meanPval", "Chr", "Start", "End")
allcinfo
outfile <- "5o-A188.ctg.mapping.summary"
write.table(allcinfo, outfile, row.names=F, quote=F, sep="\t")
chr.size <- read.delim("~/references/A188Ref1/genome/A188Ref1.length", header=F)

# plot
ymax <- max(-log10(as.numeric(as.character(allcinfo$meanPval[!is.na(allcinfo$Chr)]))), na.rm=T)
pval.oneplot(infile=outfile, yrange = c(6, ymax+2), xlab.text="chromosome",
             mar.set =c(4.5, 4.5, 2.5, 0.5), mgp.set = c(2.5, 0.8, 0),
             pvalue.cutoff = NULL, ylab.text = "-log10(p)",
             cex=1.2, main.text="Associated loci on A188",
             axisline.width=1.2, chr.set=1:10, plot.cex.lab = 1.3,
             chr.size=chr.size, order.by.chrsize=F, label.rm=NULL, cexaxis=1.3,
             cex.line.width = 0.5, saveplot=T, plot.width=7, plot.heigh=3, 
             plot.path=".", plot.filename="5o-A188.ft.mapping.result.pdf")


##############################################################################
### B73
##############################################################################
samfiles <- dir(pattern="B73uniq.sam")
allcinfo <- NULL
for (samfile in samfiles) {
  cinfo <- clustinfo(alnfile=samfile, min_major_chr=0.8, max_region_range=10000000)
  if (is.null(allcinfo)) {
    allcinfo <- cinfo
  } else {
    allcinfo <- rbind(allcinfo, cinfo)
  }
}

allcinfo <- data.frame(allcinfo)
colnames(allcinfo) <- c("Cluster", "ctgBp", "meanPval", "Chr", "Start", "End")
allcinfo
outfile <- "5o-B73.ctg.mapping.summary"
write.table(allcinfo, outfile, row.names=F, quote=F, sep="\t")

# plot
ymax <- max(-log10(as.numeric(as.character(allcinfo$meanPval[!is.na(allcinfo$Chr)]))), na.rm=T)
pval.oneplot(infile=outfile, yrange = c(6, ymax+2), xlab.text="chromosome",
             mar.set =c(4.5, 4.5, 2.5, 0.5), mgp.set = c(2.5, 0.8, 0),
             pvalue.cutoff = NULL, ylab.text = "-log10(p)",
             cex=1.2, main.text="Associated loci on B73",
             axisline.width=1.5, chr.set=1:10, plot.cex.lab = 1.3,
             chr.size=chr.size, order.by.chrsize=F, label.rm=NULL, cexaxis=1.3,
             cex.line.width = 0.5, saveplot=T, plot.width=7, plot.heigh=3, 
             plot.path=".", plot.filename="5o-B73.ft.mapping.result.pdf")

