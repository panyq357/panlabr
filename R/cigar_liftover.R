#' @title
#' Convert CIGAR string to named character vector
#'
#' @param cigar A CIGAR string from SAM file. e.g. `c("2M1I6M")`
#'
#' @return a named character vector. e.g. `c(M=2, I=1, M=6)`
#'
#' @export
#'
parse_cigar <- function(cigar) {
  number <- as.integer(strsplit(cigar, "[A-Z]")[[1]])
  letter <- strsplit(cigar, "\\d+")[[1]][-1]
  return(setNames(number, letter))
}


#' @title
#' Use CIGAR string to shift position
#'
#' @param pos_vect Positions that need to be shifted.
#' @param query_start Start position of query sequence in query genome (1-based index).
#' @param aln A GAlignments object (from GenomicAlignments) contains CIGAR field.
#'
#' @return A integer vector of new position match `pos_vect`.
#'
#' @details
#' First, find matched windows between query sequence and target sequence.
#' Then, check each window pair, shift pos from query to target.
#'
#' @export
#'
cigar_liftover <- function(pos_vect, query_start, aln) {

  cigar_vect <- parse_cigar(aln@cigar)

  matched_interval_list <- list()

  target_pointer <- aln@start
  query_pointer <- query_start

  for (i in seq_along(cigar_vect)) {
    if (names(cigar_vect)[i] == "M") {
      matched_interval_list[[i]] <- c(
        target_left = target_pointer,
        target_right = target_pointer + unname(cigar_vect[i]) - 1,
        query_left = query_pointer,
        query_right = query_pointer + unname(cigar_vect[i]) - 1
      )
      target_pointer <- target_pointer + cigar_vect[i]
      query_pointer <- query_pointer + cigar_vect[i]
    } else if (names(cigar_vect)[i] == "I") {
      query_pointer <- query_pointer + cigar_vect[i]
    } else if (names(cigar_vect)[i] %in% c("D", "N")) {
      target_pointer <- target_pointer + cigar_vect[i]
    }
  }

  matched_df <- do.call(rbind, matched_interval_list) |>
    as.data.frame()

  new_pos_vect <- rep(-1, length(pos_vect))

  for (i in seq_len(nrow(matched_df))) {
    in_interval <- matched_df$target_left[i] <= pos_vect & pos_vect <= matched_df$target_right[i]
    new_pos_vect[in_interval] <- pos_vect[in_interval] - matched_df$target_left[i] + matched_df$query_left[i]
  }

  return(new_pos_vect)
}
