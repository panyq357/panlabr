range_transform <- function(x, from_range, to_range) {
  rel_x <- (x - min(from_range)) / (max(from_range) - min(from_range))
  out_x <- rel_x * (max(to_range) - min(to_range)) + min(to_range)
  return(out_x)
}


#' @export
plot_hap_table <- function(hap_stat, genome_gr, cell_padding=1, line_length = 5) {

  nchar_list <- apply(hap_stat, 2, function(col) max(nchar(col)))

  cell_width_list <- nchar_list + 2 * cell_padding

  vertical_border_list <- c(0, cumsum(cell_width_list))

  cell_x_center_list <- (vertical_border_list[-1] + vertical_border_list[-length(vertical_border_list)]) / 2

  cell_x_center_list <- range_transform(cell_x_center_list, range(vertical_border_list), c(start(genome_gr), end(genome_gr)))

  vertical_border_list <- range_transform(vertical_border_list, range(vertical_border_list), c(start(genome_gr), end(genome_gr)))

  cell_height_list <- -(rep(1, nrow(hap_stat)) + 2 * cell_padding)

  horizontal_border_list <- c(0, cumsum(cell_height_list))

  cell_y_center_list <- (horizontal_border_list[-1] + horizontal_border_list[-length(horizontal_border_list)]) / 2

  label_coord_list <- list()

  for (i in seq_len(nrow(hap_stat))) {
    for (j in seq_len(ncol(hap_stat))) {
      label_coord_list[[length(label_coord_list) + 1]] <- list(
        label = as.character(hap_stat[[i, j]]),
        y = cell_y_center_list[i],
        x = cell_x_center_list[j]
      )
    }
  }

  label_df <- data.table::rbindlist(label_coord_list)

  pos_idx <- grep("Pos", names(hap_stat))

  header <- names(hap_stat) |>
    sub("Pos\\.", "", x = _) |>
    sub("Group\\.", "", x = _)

  header_height <- max(nchar(names(hap_stat))) + 2 * cell_padding

  horizontal_border_list[[length(horizontal_border_list) + 1]] <- header_height

  header_df <- data.frame(x = cell_x_center_list, y = header_height / 2, label = header)

  pos_x <- cell_x_center_list[pos_idx]
  pos_xend <- names(hap_stat)[pos_idx] |>
    sub("Pos\\.", "", x = _) |>
    as.integer()

  plt_obj <- ggplot() +
    geom_segment(aes(x = vertical_border_list, xend = vertical_border_list, y = min(horizontal_border_list), yend = max(horizontal_border_list))) +
    geom_segment(aes(x = min(vertical_border_list), xend = max(vertical_border_list), y = horizontal_border_list, yend = horizontal_border_list)) +
    geom_text(data = label_df, aes(x = x, y = y, label = label)) +
    geom_text(data = header_df, aes(x = x, y = y, label = label), angle = 90) +
    geom_segment(aes(x = pos_x, xend = pos_xend, y = max(horizontal_border_list), yend = max(horizontal_border_list) + line_length)) +
    theme_classic() +
    scale_y_continuous(breaks = NULL, labels = NULL) +
    labs(x = "Pos", y = "")

  return(plt_obj)
}
