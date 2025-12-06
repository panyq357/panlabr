#' Wrapper for running beagle to impute vcf
#'
#' @param gt path to .vcf.gz file.
#' @param out output prefix, output file will be `prefix.vcf.gz`.
#' @param nthreads number of threads beagle can use.
#'
#' @export
beagle_wrapper <- function(gt, out, nthreads = parallel::detectCores(), seed = 12345) {
  system2("beagle", args = c(
    sprintf("gt=%s", gt),
    sprintf("out=%s", out),
    sprintf("nthreads=%d", nthreads),
    sprintf("seed=%d", seed)
  ))
}


#' Convert vcf df to haplotype matrix
#'
#' Row is haplotype, column is variant.
#' @param vcf vcf data.frame or tibble.
#'
#' @export
vcf_to_hap_mat <- function(vcf) {
  if (length(vcf) < 10) {
    stop("No genotype in vcf, stop.")
  } else if (length(grep("/", vcf[10:ncol(vcf)])) > 0) {
    stop("vcf not phased, stop.")
  }

  # Split haplotype by `|`, into two matrix.
  # And transpose, row is sample, column is variant.
  hap1 <- sub("([^|]+)\\|([^|]+)", "\\1", vcf_to_gt_mat(vcf)) |> t()
  hap2 <- sub("([^|]+)\\|([^|]+)", "\\2", vcf_to_gt_mat(vcf)) |> t()

  # Add suffix to sample names.
  rownames(hap1) <- paste(rownames(hap1), "1", sep = ".")
  rownames(hap2) <- paste(rownames(hap2), "2", sep = ".")

  # Combine two matrix interleaving.
  hap_mat <- matrix(nrow = nrow(hap1) * 2, ncol = ncol(hap1))
  hap_mat[seq(1, nrow(hap_mat), by = 2), ] <- hap1
  hap_mat[seq(2, nrow(hap_mat), by = 2), ] <- hap2

  # Add rownames interleaving.
  rownames(hap_mat) <- character(nrow(hap_mat))
  rownames(hap_mat)[seq(1, nrow(hap_mat), by = 2)] <- rownames(hap1)
  rownames(hap_mat)[seq(2, nrow(hap_mat), by = 2)] <- rownames(hap2)
  colnames(hap_mat) <- colnames(hap1)

  return(hap_mat)
}

#' Count haplotypes per group
#'
#' @param hap_mat haplotype matrix produced by `vcf_to_hap_mat`.
#' @param group_info a named vector or factor that indicates each haplotype's group.
#'
#' @export
hap_mat_to_hap_long <- function(hap_mat, group_info) {
  colnames(hap_mat) <- paste("Pos", colnames(hap_mat), sep = ".")
  hap_long <- cbind(
    data.frame(
      SampleID = sub("(.*)\\.[12]$", "\\1", rownames(hap_mat)),
      HapID = rownames(hap_mat)
    ),
    hap_mat
  )

  if (class(group_info) == "factor") {
    for (group in levels(group_info)) {
      hap_long[[sprintf("Group.%s", group)]] <- group_info[hap_long$SampleID] == group
    }
  } else {
    for (group in unique(group_info)) {
      hap_long[[sprintf("Group.%s", group)]] <- group_info[hap_long$SampleID] == group
    }
  }

  hap_long[startsWith(names(hap_long), "Group.")][is.na(hap_long[startsWith(names(hap_long), "Group.")])] <- FALSE

  gt_str_vect <- apply(hap_long[substr(names(hap_long), 1, 4) == "Pos."], 1, function(row) paste(row, collapse="|"))

  gt_str_freq <- table(gt_str_vect) |> sort(decreasing = TRUE)

  hap_long[["HapGroup"]] <- sprintf("Hap%d", match(gt_str_vect, names(gt_str_freq)))

  return(hap_long)
}

#' @export
hap_long_to_hap_stat <- function(hap_long) {
  hap_long_stat <- split(hap_long, hap_long$HapGroup) |> lapply(function(df) {
    stat_df <- df[startsWith(names(df), "Group.")] |>
      apply(2, sum) |>
      as.data.frame() |>
      t() |>
      as.data.frame()
    stat_df <- cbind(HapGroup = df$HapGroup[[1]], df[1, startsWith(names(df), "Pos.")], stat_df, Count = rowSums(stat_df))
    return(stat_df)
  }) |>
    do.call(rbind, args=_)
  return(hap_long_stat[order(sub("[^0-9]+([0-9]+).*", "\\1", hap_long_stat$HapGroup) |> as.integer()), ])
}


#' @export
hap_stat_to_nexus_df <- function(hap_stat) {

  if (nrow(hap_stat) >= 20) {
    warning("haplotype number is ", nrow(hap_stat), ", might need filter first.")
  }

  pos_mat <- hap_stat[startsWith(names(hap_stat), "Pos.")]

  for (j in seq_len(ncol(pos_mat))) {
    if (max(nchar(pos_mat[, j])) <= 1) {
      next
    }
    if (length(table(pos_mat[, j])) > 4) {
      warning(colnames(pos_mat)[j], " have more then 4 allele, set this variants to N")
    }

    message(colnames(pos_mat)[j], " is indel, convert it to lower case a/t/c/g")

    pos_mat[, j] <- c("a", "t", "c", "g")[match(pos_mat[, j], table(pos_mat[, j]) |> sort(decreasing = TRUE) |> names())]
  }

  seq <- apply(pos_mat, 1, function(row) paste(row, collapse=""))

  nexus_df <- cbind(ID=hap_stat$HapGroup, Seq=seq, hap_stat[startsWith(names(hap_stat), "Group.")])
  row.names(nexus_df) <- NULL

  return(nexus_df)
}


#' @export
write_nexus <- function(nexus_df, out) {
  trait_cols <- grep("Group.", names(nexus_df))

  ntax <- nrow(nexus_df)
  taxlabels <- paste(nexus_df$ID, collapse = "\n")
  nchar_ <- nchar(nexus_df$Seq[[1]])
  seq_matrix <- sprintf("%s %s", nexus_df$ID, nexus_df$Seq) |> paste(collapse = "\n")
  ntraits <- length(trait_cols)
  trait_labels <- paste(names(nexus_df)[trait_cols], collapse = " ")
  trait_matrix <- paste(
    sprintf("%s %s", nexus_df$ID, apply(nexus_df[trait_cols], 1, paste, collapse = ",")),
    collapse = "\n"
  )
  nexus_string <- glue::glue(
    "#NEXUS
BEGIN TAXA;
DIMENSIONS NTAX={ntax};
TAXLABELS
{taxlabels}
;
END;

BEGIN CHARACTERS;
DIMENSIONS NCHAR={nchar_};
FORMAT DATATYPE=DNA MISSING=? GAP=- ;
MATRIX
{seq_matrix}
;
END;

BEGIN TRAITS;
Dimensions NTRAITS={ntraits};
Format labels=yes missing=? separator=Comma;
TraitLabels {trait_labels};
Matrix
{trait_matrix}
;
END;
"
  )
  write(nexus_string, out)
}

#' @export
plot_hap_table <- function(hap_stat, range) {

  hap_stat_to_label_df <- function(hap_stat, range) {

    # TODO: Split cells by nchar(cell)

    col_rel_width <- apply(hap_stat, 2, function(col) max(nchar(col)))

    inc <- (max(range) - min(range)) / sum(col_rel_width)

    # e.g. if 3 cell in a row, split row into 6 segment,
    # let cell center letter in the center of two segment.

    label_list <- list()
    x_list <- list()
    y_list <- list()

    for (i in seq_len(nrow(hap_stat))) {
      for (j in seq_len(ncol(hap_stat))) {
        if (j == 1) {
          x_list[[length(x_list)+1]] <- min(range) +
            col_rel_width[j] / 2 * inc
        } else {
          x_list[[length(x_list)+1]] <- min(range) +
            sum(col_rel_width[1:j-1]) * inc +
            col_rel_width[j] / 2 * inc
        }
        y_list[[length(y_list)+1]] <- -i
        label_list[[length(label_list)+1]] <- hap_stat[[i, j]]
      }
    }

    label_df <- data.frame(
      X = unlist(x_list),
      Y = unlist(y_list),
      Label = unlist(label_list)
    )

    return(label_df)
  }

  hap_stat_to_header_df <- function(hap_stat, range) {
    col_rel_width <- apply(hap_stat, 2, function(col) max(nchar(col)))
    inc <- (max(range) - min(range)) / sum(col_rel_width)
    label_list <- list()
    x_list <- list()
    y_list <- list()

    for (j in seq_len(ncol(hap_stat))) {
      if (j == 1) {
        x_list[[length(x_list)+1]] <- min(range) +
          col_rel_width[j] / 2 * inc
      } else {
        x_list[[length(x_list)+1]] <- min(range) +
          sum(col_rel_width[1:j-1]) * inc +
          col_rel_width[j] / 2 * inc
      }
      y_list[[length(y_list)+1]] <- 0
      label <- sub("^(Pos.|Group.)", "", names(hap_stat)[j])
      label_list[[length(label_list)+1]] <- label
    }

    header_df <- data.frame(
      X = unlist(x_list),
      Y = unlist(y_list),
      Label = unlist(label_list)
    )

    return(header_df)
  }

  hap_stat_to_line_df <- function(hap_stat, range, y=3, yend=5) {
    col_rel_width <- apply(hap_stat, 2, function(col) max(nchar(col)))
    inc <- (max(range) - min(range)) / sum(col_rel_width)

    x_list <- list()
    y_list <- list()
    xend_list <- list()
    yend_list <- list()

    for (j in seq_len(ncol(hap_stat))) {

      name <- names(hap_stat)[j]

      if (!startsWith(name, "Pos.")) {
        next
      }

      if (j == 1) {
        x_list[[length(x_list)+1]] <- min(range) +
          col_rel_width[j] / 2 * inc
      } else {
        x_list[[length(x_list)+1]] <- min(range) +
          sum(col_rel_width[1:j-1]) * inc +
          col_rel_width[j] / 2 * inc
      }

      xend_list[[length(xend_list)+1]] <- sub("Pos.", "", name) |> as.integer()
      y_list[[length(y_list)+1]] <- y
      yend_list[[length(yend_list)+1]] <- yend
    }

    line_df <- data.frame(
      X = unlist(x_list),
      Y = unlist(y_list),
      XEnd = unlist(xend_list),
      YEnd = unlist(yend_list)
    )

    return(line_df)
  }

  hap_stat_to_grid_df <- function(hap_stat, range, top=3) {

    col_rel_width <- apply(hap_stat, 2, function(col) max(nchar(col)))
    inc <- (max(range) - min(range)) / sum(col_rel_width)

    x_list <- list(min(range))
    for (j in seq_along(col_rel_width)) {
      x_list[[length(x_list)+1]] <- min(range) +
        sum(col_rel_width[1:j]) * inc
    }

    vline_df <- data.frame(
      X = unlist(x_list),
      XEnd = unlist(x_list),
      Y = top,
      YEnd = - nrow(hap_stat) - 0.5
    )

    y_list <- c(top, - seq_len(nrow(hap_stat)) + 0.5, - nrow(hap_stat) - 0.5)

    hline_df <- data.frame(
      Y = unlist(y_list),
      YEnd = unlist(y_list),
      X = min(range),
      XEnd = max(range)
    )

    return(rbind(vline_df, hline_df))
  }

  label_df <- hap_stat_to_label_df(hap_stat, range)
  header_df <- hap_stat_to_header_df(hap_stat, range)
  line_df <- hap_stat_to_line_df(hap_stat, range)
  grid_df <- hap_stat_to_grid_df(hap_stat, range)

  hap_table <- ggplot() +
    geom_text(data=label_df, mapping=aes(x=X, y=Y, label=Label)) +
    geom_text(data=header_df, mapping=aes(x=X, y=Y, label=Label), angle=90, hjust=0, vjust=0.5) +
    geom_segment(data=line_df, mapping=aes(x=X, y=Y, xend=XEnd, yend=YEnd)) +
    geom_segment(data=grid_df, mapping=aes(x=X, y=Y, xend=XEnd, yend=YEnd)) +
    theme_void() +
    scale_x_continuous(limits=range, expand=c(0, 0))

  return(hap_table)
}
