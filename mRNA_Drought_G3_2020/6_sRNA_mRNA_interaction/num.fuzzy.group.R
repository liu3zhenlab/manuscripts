num.fuzzy.group <- function (category, num, distance=10) {
  
  sortdiff <- function(x) {
    group <- 1
    curgroup <- 1
    if (length(x) > 1) {
      ### sort
      x <- sort(x)
      ### diff
      xdiff <- diff(x)
      gap.points <- (xdiff > distance)
      ### group
      for (i in 2:length(x)) {
        if (gap.points[i-1]) {
          curgroup <- curgroup + 1
        }
        group <- c(group, curgroup)
      }
    }
    data.frame(data=x, group=group)
  }
  
  ### by category
  num.groups <- tapply(num, category, sortdiff)
  nnum.groups <- do.call(rbind.data.frame, num.groups)
  nnum.groups$category <- gsub("\\.[0-9]+$", "", rownames(nnum.groups))
  rownames(nnum.groups) <- 1:nrow(nnum.groups)
  nnum.groups[, c("category", "data", "group")]
}
