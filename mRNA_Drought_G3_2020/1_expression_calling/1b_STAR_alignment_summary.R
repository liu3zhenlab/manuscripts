kstarALNsummary <- function(datapath, outfile, suffix = "Log.final.out", simple.out = F) {
  ### script to summarize STAR alignments
  ### Sanzhen Liu
  ### 5/4/2017
  ### datapath: directory of STAR outputs
  ### suffix: suffix of file names
  ### version 0.01
  ### function: kstarALNsummary("datapath","outfile")
  
  summary.files <- dir(path = datapath, pattern = suffix)
  
  allsummary <- NULL
  for (sf in summary.files) {
    data <- read.delim(paste0(datapath, "/", sf), header = F)
    datainfo <- gsub(suffix, "", sf)
    colnames(data) <- c("Item", datainfo)
    
    data <- data[grep("\\|$", data$Item), ]
    data$Item <- gsub(" \\|", "", data$Item)
    
    if (is.null(allsummary)) {
      allsummary <- data
    } else {
      allsummary[, datainfo] <- data[, 2]
    }
  }
  
  rownames(allsummary) <- 1:nrow(allsummary)
  
  if (simple.out) {
  ### simple output, only containing four columns
    allsummary <- allsummary[c(5, 7, 8), ]
    rownames(allsummary) <- c("TotalInputReads", "UniquelyMappedReads", "Unique_Percentage")
    allsummary <- t(allsummary[, -1])
    allsummary <- data.frame(allsummary)
    allsummary$Sample <- rownames(allsummary)
    allsummary <- allsummary[, c("Sample", "TotalInputReads", "UniquelyMappedReads", "Unique_Percentage")]
    for (i in 2:3) {
      allsummary[, i] <- as.numeric(as.character(allsummary[, i]))
    }
  }
  
  ### output
  write.table(allsummary,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
}
