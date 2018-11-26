################################################
### BSR-Seq module
### 8/9/2016
## This is the code for the BSR-Seq analysis
## Author: Sanzhen Liu and Dan Nettleton
## Note: core code was developed by Dan Nettleton
################################################

#======== Start: data subject to change =========#
dan.bsa <- function (ac, out.path=".", out.file,
                     chr.colname="CHR", pos.colname="POS",
                     mut.ref.colname="mutREF", mut.alt.colname="mutALT",
                     wt.ref.colname="wtREF", wt.alt.colname="wtALT",
                     wt.ref.min=3, wt.alt.min=3, mut.total.min=5,
		     mut.ind.num, total.genetic.length=2338,
                     genetic.interval=20) {
## filter SNPs:
  ac <- ac[!is.na(ac[, wt.ref.colname]) &
           !is.na(ac[, wt.alt.colname]) &
	   !is.na(ac[, mut.ref.colname]) &
	   !is.na(ac[, mut.alt.colname]), ]
  ac <- ac[ac[, wt.ref.colname] >= wt.ref.min
          & ac[, wt.alt.colname] >= wt.alt.min
          & (ac[, mut.ref.colname] + ac[, mut.alt.colname]) >= mut.total.min, ]
  cat("After SNP filtering, the number of remaining SNPs is: ")
  cat(nrow(ac))

  ##Find estimated probability that a randomly selected
  ##SNP is in complete linkage disequilibrium with mutant gene.
  dtor <- function(x) {
  # convert genetic distance to recombination rate:
    r = 0.5 * (1 - exp(-2 * x)) # Haldane's mapping function, x is unit in Morgan
    return(r)
  }

  cM <- c(seq(0, genetic.interval, by=0.01)) # with 20cM interval
  rf <- dtor(cM/100)  # recombination rate
  pnr <- (1-rf)^(2*mut.ind.num)  # the prior probability of no recombination 
                                # between the SNP and the causal gene

  ##This is numerical integration of P(no recomb|dist)P(dist)
  ptheta0=sum(2*(0.01/total.genetic.length)*(pnr[-length(pnr)]+pnr[-1])/2)

  ##Find prior distribution of theta
  wt.total <- apply(ac[, c(wt.ref.colname, wt.alt.colname)], 1, sum)
  thetahat <- ac[, wt.ref.colname]/wt.total # prior allele frequency
  thetahat <- c(thetahat, 1-thetahat)

  ### function:
  getppp <- function(x) {
    nm=sum(x[3:4]) # Mutant alleles
    nw=sum(x[1:2]) # Wt alleles
    nonma=which.min(x[3:4]) # allele with the smaller number
    xm=min(x[3:4])
    xw=(x[1:2])[nonma]
    ##R="some recombinaton between SNP and causal gene in mutant pool"
    if (xm>0) {
      ppmt=0
      ppwt=0
    } else { # no recombination
    ##Find P(x1,x2|R)
    ##Integrate P(x1=0 or n|theta)prior(theta)
      p1=mean((1-thetahat)^nm+thetahat^nm)
    ##Find P(R|x1,x2)
      ppmt=(1/(1+p1*((1-ptheta0)/ptheta0)))
    ##Find P(wm,w) by integrating P(wm,w|theta)prior(theta) 
      p2=mean(dbinom(xw,nw,thetahat))
    ##Find P(wm,w|R)
      p3=mean(dbinom(xw,nw,thetahat[thetahat>=0.5]))
    ##Find P(R|wm,w)
      ppwt=.5*p3/p2
    ##Find the product of the posterior probabilities.
    }
    ppp=ppmt*ppwt
    c(ppmt,ppwt,ppp)
  } # end of the core function code
 
  o = t(apply(ac[, c(wt.ref.colname, wt.alt.colname,
                     mut.ref.colname, mut.alt.colname)], 1, getppp))
##Columns are 
##posterior probability from mutant data,
##posterior probability from wild type data,
##product of the posterior probabilities.
  ac2 <- ac[, c(chr.colname, pos.colname, mut.ref.colname,
         mut.alt.colname, wt.ref.colname, wt.alt.colname)]
  out = data.frame(ac2, o)
  names(out)[(ncol(out)-2):ncol(out)]=c("ppmt","ppwt","ppp")
  
  out =  out[order(out$ppp, decreasing = T), ]
  
  cat("Five top SNPs with the highest ppp values:\n")
  print(out[1:5, ], )
 
  out.path = gsub("/$", "", out.path)
  write.table(out, paste(out.path, "/", out.file, sep=""), 
              row.names=F, quote=F, sep="\t")
  invisible(out)
}

