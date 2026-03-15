#' Plot scatter pie map.
#'
#' @param plot_data a long data.frame, contains these columns: SampleID, HapGroup, lon, lat.
#' @param sf_obj a sf object read from a shapefile by [sf::read_sf()].
#' @return a ggplot plot object.
#' @importFrom scatterpie geom_scatterpie geom_scatterpie_legend
#' @export
plot_scatter_pie_map <- function(plot_data, sf_obj) {

  plot_data$pos <- paste(plot_data$lon, plot_data$lat)

  hap_names <- unique(plot_data$HapGroup)

  pie_data <- split(plot_data$HapGroup, plot_data$pos) |>
    lapply(table) |>
    lapply(as.list) |>
    data.table::rbindlist(fill=TRUE, idcol="pos") |>
    as.data.frame() |>
    (function(x) {
      x[is.na(x)] <- 0
      return(x)
    })() |>
    dplyr::left_join(dplyr::distinct(plot_data[c("pos", "lon", "lat")]), by="pos")

  pie_data$radius <- log10(rowSums(pie_data[hap_names])) + 1

  plt_obj <- ggplot() +
    geom_sf(data = sf_obj) +
    geom_scatterpie(data=pie_data, aes(x=lon, y=lat, r=radius/1.5), cols=hap_names) +
    theme_bw() +
    geom_scatterpie_legend(log10(hap_pie_data$N), breaks=(log10(c(1, 10, 100)) + 1) / 1.5, x=73, y=-30, labeller=function(x) 10 ^ (x*1.5-1), n=3) +
    labs(x = "Longitude", y = "Latitude", fill = "Haplotype")

  return(plt_obj)
}
