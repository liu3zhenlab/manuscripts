nucmer.twosets.plot <- function(datapath=".", datafile, lend.turnoff=F,
                        aln.min.size=100, aln.min.identity=90,
                        xlabel.rm = "chr", ylabel.rm = "chr",
                        tableout=F, pdfout=F, outpath=".", line.width.factor=5,
                        xrange=NULL, yrange=NULL,
			pdf.width = 10, pdf.height = 10,
                        tableoutfile="seq.nucmer.txt", imageoutfile="seq.nucmer.pdf", ...) {
  ### datafile is the nucmer output .delta
  tmp.out.file <- paste(".nucmer.aln.", gsub(" ", "", date()), ".table.txt", sep="")
  indata <- paste(datapath, "/", datafile, sep="")
  sc.cmd <- paste("show-coords -clT -q -L", aln.min.size, "-I",
                  aln.min.identity, indata, ">", tmp.out.file)
  system(sc.cmd) ### run show-coords to generate alignment output
  co <- read.delim(tmp.out.file, comment.char="#", header=F,stringsAsFactor=F, skip=4)
  # [S1]  [E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[TAGS]
  colnames(co) <- c("firsts", "firste", "seconds", "seconde",
                    "first.align.l", "second.align.l", "ident",
                    "firstl", "secondl", "first.cov",
                    "second.cov", "first", "second")
  co$firsts <- as.numeric(as.character(co$firsts))
  co <- co[order(co$first, co$firsts, co$second, co$seconds), ]

  ### determine the total size of all the contigs in each set
  co$first <- as.character(co$first)
  co$second <- as.character(co$second)
  first.contigs.size <- tapply(co$firstl, co$first, max)
  second.contigs.size <- tapply(co$secondl, co$second, max)
  second.contigs.size <- second.contigs.size[unique(co$second)]
  
  print(second.contigs.size)
  
  
  if (tableout) {
    write.table(co, tableoutfile, quote=F, row.names=F, sep="\t")
  } else {
    print(co)
  }

  system(paste("rm", tmp.out.file))
  ### maximum values for each contig:
  
  nfirst <- length(first.contigs.size)
  first.accum <- rep(0, times=nfirst)
  if (nfirst>1) {
    for (i in 2:nfirst) {
      first.accum[i] <- sum(first.contigs.size[1:(i-1)])   
    }
  } else {
    first.accum <- 0
  }
  names(first.accum) <- names(first.contigs.size)
  print(head(first.accum))

  ### second set:
  nsecond <- length(second.contigs.size)
  second.accum <- rep(0, times=nsecond)
  if (nsecond>1) {
    for (i in 2:nsecond) {
      second.accum[i] <- sum(second.contigs.size[1:(i-1)])   
    }
  } else {
    second.accum <- 0
  }
  names(second.accum) <- names(second.contigs.size)
  print(head(second.accum))
  
  max.size <- max(sum(first.contigs.size), sum(second.contigs.size))
  max.xsize <- max(sum(first.contigs.size))
  max.ysize <- max(sum(second.contigs.size))

  co$first.accum.s <- co$firsts + as.numeric(first.accum[co$first])
  co$first.accum.e <- co$firste + as.numeric(first.accum[co$first])
  co$second.accum.s <- co$seconds + as.numeric(second.accum[co$second])
  co$second.accum.e <- co$seconde + as.numeric(second.accum[co$second])

  ###
  ### plot
  ###
  if (pdfout) {
    outpdffile <- paste(outpath, "/", imageoutfile, sep="")
    pdf(outpdffile, width=pdf.width, height=pdf.height)
  }
  
  if (is.null(xrange)) {
  	xrange <- c(-sum(first.contigs.size) / 50, sum(first.contigs.size))
  }

  if (is.null(yrange)) {
  	yrange <- c(0, sum(second.contigs.size))
  }

  ### plot
  plot(NULL, NULL, type="n", xlim=xrange, ylim=yrange, ...)
  col1 <- rgb(0, 0.4, 0, 0.5)
  col2 <- rgb(1, 0, 0, 0.5)
  ctg.num <- length(first.contigs.size)
  col.db <- data.frame(Col=rep(c(col1, col2), ctg.num)[1:ctg.num], Contig=names(first.contigs.size))
  
  abline(v=first.accum[1], lwd=0.6)
  abline(v=first.accum[-1], lwd=1, col="light grey")
  abline(v=sum(first.contigs.size), lwd=0.6)
  
  abline(h=second.accum[1], lwd=0.6)
  abline(h=second.accum, lwd=1, col="light grey")
  abline(h=sum(second.contigs.size), lwd=0.6)
  
  for (i in 1:nrow(co)) {
    plot.col <- as.character(col.db$Col[col.db$Contig==co[i, "first"]])
    distance <- co$second.accum.e[i] - co$second.accum.s[i]
    lend.val <- 2
    if (distance > max.size/50) {
      lend.val <- 1
    }
    if (lend.turnoff) {
      lend.val <- 1
    }
    
    lines(c(co$first.accum.s[i], co$first.accum.e[i]),
          c(co$second.accum.s[i], co$second.accum.e[i]),
          lwd=lend.val*line.width.factor, col= plot.col, lend=lend.val)
  }
  
  first.coord <- (c(first.accum[-1], sum(first.contigs.size)) + first.accum) / 2
  cat(first.coord)
  second.coord <- (c(second.accum[-1], sum(second.contigs.size)) + second.accum) / 2
  cat(second.coord)
  xlabels <- names(first.accum)
  xlabels <- gsub(xlabel.rm, "", xlabels)
  ylabels <- names(second.accum)
  cat("--", ylabels, "\n")
  ylabels <- gsub(ylabel.rm, "", ylabels)
  cat("==", ylabels, "\n")
  text(x=first.coord, y= - max.ysize / 50, labels=xlabels, cex = 0.8, xpd = T)
  text(x= - max.xsize / 40, y=second.coord, labels=ylabels, cex = 0.8, xpd = T)
  if (pdfout) {
    dev.off()
  }
}

