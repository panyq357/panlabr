#' Read VCF file as tibble.
#' @param vcf_path Path to VCF file.
#' @export
setGeneric("read_vcf",
  function(vcf_path, chromosome, start, end) {
    standardGeneric("read_vcf")
  }
)


setMethod("read_vcf",
  signature(
    vcf_path = "character",
    chromosome = "missing",
    start = "missing",
    end = "missing"
  ),
  function(vcf_path) {
    header <- readr::read_lines(vcf_path)
    header <- header[1:(grep("#CHROM", header) - 1)]

    vcf <- readr::read_tsv(
      file = vcf_path,
      comment = "##",
      col_types = readr::cols(
        POS = readr::col_integer(),
        .default = readr::col_character()
      )
    )

    attr(vcf, "header") <- header

    return(vcf)
  }
)


setMethod("read_vcf",
  signature(
    vcf_path = "character",
    chromosome = "character",
    start = "numeric",
    end = "numeric"
  ),
  function(vcf_path, chromosome, start, end) {

    if (!file.exists(paste0(vcf_path, ".tbi"))) {
      system2("tabix", args = vcf_path)
    }

    tmp_vcf <- tempfile(fileext = ".vcf")
    system2("tabix",
      args = c("-h", vcf_path, sprintf("%s:%d-%d", chromosome, start, end)),
      stdout = tmp_vcf
    )
    vcf <- read_vcf(tmp_vcf)
    unlink(tmp_vcf)
    return(vcf)
  }
)


#' Write VCF tibble as bgzipped file.
#' @param vcf A VCF tibble.
#' @param file Path to output bgzipped file.
#' @importFrom readr write_lines write_tsv
#' @export
write_vcf_bgzip <- function(vcf, file) {

  tmp_vcf <- tempfile(fileext = ".vcf")

  readr::write_lines(attr(vcf, "header"), tmp_vcf)
  readr::write_tsv(vcf, tmp_vcf, append = TRUE, col_names = TRUE)

  system2("bgzip", stdin = tmp_vcf, stdout = file)

  unlink(tmp_vcf)
}


#' Get genotype character matrix from vcf tibble.
#' @param vcf_df vcf tibble or data.frame.
#' @param simplify e.g. "T/T" > "T", "AG/AG" > "AG", won't simplify heterozygote.
#' @export
vcf_to_gt_mat <- function(vcf_df, simplify = FALSE) {

  # Empty matrix.
  gt_mat <- matrix(nrow = nrow(vcf_df), ncol = ncol(vcf_df) - 9)
  colnames(gt_mat) <- names(vcf_df)[10:length(vcf_df)]
  rownames(gt_mat) <- vcf_df$POS

  # Add to matrix row by row.
  for (i in seq_len(nrow(vcf_df))) {
    ref_alt <- c(vcf_df$REF[i], strsplit(vcf_df$ALT[i], ",")[[1]])
    names(ref_alt) <- 0:(length(ref_alt) - 1)

    gt <- vcf_df[i, 10:length(vcf_df)] |>
      unlist() |>
      sub("([^:]+).*", "\\1", x = _)  # To remove fields after genotype (sep by `:`).

    for (num in names(ref_alt)) {
      gt <- gsub(num, ref_alt[num], gt)
    }

    if (simplify == TRUE) {
      gt <- gsub("([^/|]+)[/|]\\1$", "\\1", gt)  # e.g. "T/T" > "T", "A/T" > "A/T"
    }

    gt_mat[i,] <- gt
  }

  return(gt_mat)
}

#' Calculate MAF and Missing rate of variants in vcf df
#'
#' @param vcf a vcf tibble.
#'
#' @export
vcf_variant_stats <- function(vcf) {
  missing_list <- list()
  maf_list <- list()
  for (i in seq_len(nrow(vcf))) {
    t <- vcf[i, 10:length(vcf)] |>
      unlist() |>
      strsplit("[/|]") |>
      unlist() |>
      table() |>
      sort(decreasing = TRUE)
    missing_list[[length(missing_list) + 1]] <- t[names(t) == "."] / sum(t)
    maf_list[[length(maf_list) + 1]] <- t[names(t) != "."][2] / sum(t[names(t) != "."])
  }

  out <- cbind(
    vcf[1:2],
    Missing = unlist(missing_list),
    MAF = unlist(maf_list)
  )

  return(out)
}
