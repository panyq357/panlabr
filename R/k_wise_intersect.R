all_intersect <- function(gr_list) {
  if (length(gr_list) == 1) {
    return(gr_list[[1]])
  } else {
    return(Reduce(GenomicRanges::intersect, gr_list))
  }
}

all_union <- function(gr_list) {
  if (length(gr_list) == 1) {
    return(gr_list[[1]])
  } else {
    return(Reduce(GenomicRanges::union, gr_list))
  }
}


#' K-wise intersection of GRanges
#' @param gr_list a list of GRanges
#' @param k e.g. a 2-wise intersection of 3 GRanges means 3 choose 2.
#' @export
k_wise_intersect <- function(gr_list, k) {
  combn(gr_list, k, all_intersect, simplify=FALSE) |> all_union()
}
