write_bed <- function(gr, file, columns) {
  df <- as.data.frame(gr)[columns]
  levels(df$strand)[levels(df$strand) == "*"] <- "."
  df$start <- df$start - 1  # In BED file, ranges are 0-indexed, left-closed, right-open.
  write.table(df, file, quote = FALSE, sep = "\t", na = ".", row.names = FALSE, col.names = FALSE)
}


#' Export Granges as narrowPeak file
#' @param gr GRanges imported from a narrowPeak file using [rtracklayer::import()]
#' @param file path to output file
#' @export
write_narrowPeak <- function(gr, file) {
  write_bed(gr, file, c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
}


#' Export Granges as narrowPeak file
#' @param gr GRanges imported from a narrowPeak file using [rtracklayer::import()]
#' @param file path to output file
#' @export
write_broadPeak <- function(gr, file) {
  write_bed(gr, file, c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue"))
}
