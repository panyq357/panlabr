#' Export Granges as narrowPeak file
#' @param gr GRanges imported from a narrowPeak file
#' @param file path to output file
#' @export
write_narrowPeak <- function(gr, file) {
  df <- as.data.frame(gr)[c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")]
  levels(df$strand)[levels(df$strand) == "*"] <- "."
  df$start <- df$start - 1  # In BED file, ranges are 0-indexed, left-closed, right-open.
  write.table(df, file, quote = FALSE, sep = "\t", na = ".", row.names = FALSE, col.names = FALSE)
}
