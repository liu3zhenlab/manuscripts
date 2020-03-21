png.sl <- function (figname, outpath = ".", prefix = "", dwidth = 1.5, dheight = 1.5, suffix=".png", ...) {
  out <- paste0(outpath, "/", prefix, "_", figname, suffix)
  png(out, width = dwidth, height = dheight, res = 1200, units = "in", pointsize = 4, ...)
}

write.table.sl <- function(variable, tablename, outpath = ".", prefix = "", suffix=".txt", ...) {
  out <- paste0(outpath, "/", prefix, "_", tablename, suffix)
  write.table(variable, out, quote = F, row.names = F, sep = "\t", ...)
}
