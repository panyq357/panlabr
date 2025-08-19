make_pos_string <- function(start, end, strand) {
  start_dot_dot_end <- paste(start, end, sep = "..")
  pos_list <- list()
  for (i in seq_along(start_dot_dot_end)) {
    if (strand[i] != "-") {
      pos_list[[i]] <- start_dot_dot_end[i]
    } else {
      pos_list[[i]] <- sprintf("complement(%s)", start_dot_dot_end[i])
    }
  }
  pos_string <- sprintf("join(%s)", paste(pos_list, collapse = ","))
  return(pos_string)
}


make_feature_string <- function(type, start, end, strand, gene, standard_name) {
  lines <- list(
    sprintf("     %-16s%s", type, make_pos_string(start, end, strand)),
    sprintf('                     /gene="%s"', gene),
    sprintf('                     /standard_name="%s"', standard_name)
  )
  return(paste(lines, collapse = "\n"))
}

#' @export
setGeneric("make_genbank_string",
  function(genome, gtf, chromosome, start, end) {
    standardGeneric("make_genbank_string")
  }
)

setMethod("make_genbank_string",
  signature(
    genome = "DNAStringSet",
    gtf = "GRanges",
    chromosome = "character",
    start = "numeric",
    end = "numeric"
  ),
  function(genome, gtf, chromosome, start, end) {
    # 1. subset genome and gtf, and shift gtf pos.
    gr <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(start, end))
    local_genome <- genome[gr][[1]]
    local_gtf <- gtf[IRanges::`%over%`(gtf, gr)] |>
      GenomicRanges::shift(1 - start)

    # 2. make feature string,
    #    split gtf by transcripts,
    #    loop over transcripts to make cds string.
    splitted <- split(local_gtf, local_gtf$transcript_id)
    feature_string_list <- list()
    for (gr in splitted) {
      cds <- subset(gr, gr$type == "CDS")
      if (length(cds) > 0) {
        feature_string_list[[length(feature_string_list) + 1]] <- make_feature_string(
          "CDS",
          BiocGenerics::start(cds),
          BiocGenerics::end(cds),
          BiocGenerics::strand(cds) |> as.character(),
          unique(cds$gene_id),
          unique(cds$transcript_id)
        )
      }
      exon <- subset(gr, gr$type == "exon")
      if (length(exon) > 0) {
        feature_string_list[[length(feature_string_list) + 1]] <- make_feature_string(
          "mRNA",
          BiocGenerics::start(exon),
          BiocGenerics::end(exon),
          BiocGenerics::strand(exon) |> as.character(),
          unique(exon$gene_id),
          unique(exon$transcript_id)
        )
      }
    }
    feature_string <- paste(feature_string_list, collapse = "\n")

    # 3. make origin seq string,
    #    split seq to 10 bp pieces,
    #    loop over, add row numbers and newlines.
    start_vect <- seq(1, length(local_genome), by = 10)
    end_vect <- c((start_vect + 9)[-length(start_vect)], length(local_genome))
    ten_bp_list <- lapply(seq_along(start_vect), function(i) {
      Biostrings::subseq(local_genome, start_vect[i], end_vect[i]) |> as.character()
    })
    origin_list <- lapply(seq_along(ten_bp_list), function(i) {
      if (i %% 6 == 1) {
        x <- sprintf("\n%9d %s", (i %/% 6) * 60 + (i %% 6), ten_bp_list[[i]])
      } else {
        x <- sprintf(" %s", ten_bp_list[[i]])
      }
      return(x)
    })
    origin_string <- paste(c(list("ORIGIN"), origin_list), collapse = "")

    name <- sprintf("%s:%d-%d", chromosome, start, end)

    final_string <- list(
      sprintf("LOCUS       %s %d bp", name, length(local_genome)),
      sprintf("DEFINITION  %s", name),
      "FEATURES             Location/Qualifiers",
      feature_string,
      origin_string,
      "//\n"
    ) |> paste(collapse = "\n")

    return(final_string)
  }
)

setMethod("make_genbank_string",
  signature(
    genome = "character",
    gtf = "character",
    chromosome = "character",
    start = "numeric",
    end = "numeric"
  ),
  function(genome, gtf, chromosome, start, end) {
    make_genbank_string(
      Biostrings::readDNAStringSet(genome),
      rtracklayer::import(gtf),
      chromosome,
      start,
      end
    )
  }
)
