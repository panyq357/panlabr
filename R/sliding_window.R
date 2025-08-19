#' @export
sliding_window <- function(df, window, step, fun, chr_col = "#CHROM", pos_col = "POS") {
  df_by_chr <- split(df, df[[chr_col]])

  chr_range_list <- lapply(df_by_chr, function(x) range(x$POS))

  res_list <- list()

  for (chr in names(df_by_chr)) {
    window_start <- seq(
      from = (chr_range_list[[chr]][[1]] %/% step) * step,
      to = (((chr_range_list[[chr]][[2]] - window) %/% step) + 1) * step,
      by = step
    )
    window_end <- window_start + window

    res_list[[chr]] <- data.frame(
      CHROM = chr,
      BIN_START = window_start,
      BIN_END = window_end,
      VALUE = slider::hop_index_vec(df_by_chr[[chr]], df_by_chr[[chr]][[pos_col]], window_start, window_end, fun)
    )
  }

  res_df <- do.call(rbind, res_list)
  row.names(res_df) <- NULL

  return(res_df)
}


#' @export
plot_sliding_window <- function(
    df,
    chr_col = "CHROM",
    start_col = "BIN_START",
    end_col = "BIN_END",
    stat_col = "VALUE") {
  ggplot2::ggplot(df) +
    ggplot2::geom_line(ggplot2::aes(x = (.data[[start_col]] + .data[[end_col]]) / 2, y = .data[[stat_col]])) +
    ggplot2::facet_wrap(~ .data[[chr_col]], nrow = 1, scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1E8, by = 1E7), labels = sprintf("%dM", seq(0, 100, by = 10)))
}
