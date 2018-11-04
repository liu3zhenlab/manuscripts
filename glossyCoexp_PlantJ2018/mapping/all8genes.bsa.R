####################################################################
### manuscript glossy
####################################################################
snps <- read.delim("snp.variants.select.txt")
cn <- colnames(snps)
cat(cn, "\n")

### loading modules:
source("bsa.R")

### main
genes <- paste("gl", c(1, 2, 3, 4, 6, 8, 26, 28), sep="")
for (gene in genes) {
  for (rep_num in c(1, 2, 3)) {
    out.file <- paste(gene, ".rep", rep_num, ".bsa.txt", sep="")
	mut <- paste0(gene, "MUT.rep", rep_num)
	wt <- paste0(gene, "WT.rep", rep_num)
	cat("\n")
	if (sum(grepl(mut, cn)) == 2 & sum(grepl(wt, cn)) == 2) {
		bsa <- dan.bsa(ac=snps, out.path=".", out.file=out.file,
        	           total.genetic.length=2000, genetic.interval=20,
            	       chr.colname="CHR", pos.colname="POS",
     	       	       mut.ref.colname=paste0(mut,  "_REF"),
				       mut.alt.colname=paste0(mut,  "_ALT"),
          	  	       wt.ref.colname=paste0(wt,  "_REF"),
				       wt.alt.colname=paste0(wt,  "_ALT"),
				       wt.ref.min=3, wt.alt.min=3,
				       mut.total.min=5, mut.ind.num=30)
	   
	}
  }
}

