starRCmerger <- function(datapath, outfile, suffix = "ReadsPerGene.out.tab") {
### script to merge all counts data from STAR output
### Sanzhen Liu
### 5/4/2017
### datapath: directory for read count files
### suffix: suffix of file names
### function: starRCmerger("datapath", "outfile)
  count.files <- dir(path = datapath, pattern = suffix)
  
  ### merge all counts
  allcounts <- NULL
  for (cf in count.files) {
    counts <- read.delim(paste0(datapath, "/", cf), header = F, stringsAsFactors = F, skip = 4)
    base <- gsub(suffix, "", cf)
    counts <- counts[, 1:2]
    colnames(counts) <- c("Gene", base)
    
    ### merge data
    if (is.null(allcounts)) {
      allcounts <- counts
    } else {
      allcounts <- merge(allcounts, counts, by = "Gene")
    }
  }
  
  ### output
  write.table(allcounts,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
}
