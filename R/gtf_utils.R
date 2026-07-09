#' Convert GFF GRanges to GTF GRanges.
#'
#' @param gff a GRanges object read by [rtracklayer::import()]
#' @param gene_type string for gene feature type. Default is `"gene"`.
#' @param transcript_type string for transcript feature type. Default is `"mRNA"`.
#'
#' @export
#'
#' @importFrom GenomicRanges mcols GRangesList
#'
gff_to_gtf <- function(gff, gene_type="gene", transcript_type="mRNA") {

  gene <- gff[gff$type == gene_type, ]
  gene$gene_id <- as.character(gene$ID)
  GenomicRanges::mcols(gene) <- GenomicRanges::mcols(gene)[c("source", "type", "score", "phase", "gene_id")]

  transcript <- gff[gff$type == transcript_type, ]
  transcript$transcript_id <- as.character(transcript$ID)
  transcript$gene_id <- as.character(transcript$Parent)
  GenomicRanges::mcols(transcript) <- GenomicRanges::mcols(transcript)[c("source", "type", "score", "phase", "transcript_id", "gene_id")]

  sub_tx_feature_names <- setdiff(levels(gff$type), c("gene", "mRNA"))
  sub_tx_feature_list <- sub_tx_feature_names |>
    lapply(function(feature_name) {
      feature <- gff[gff$type == feature_name, ]
      feature$transcript_id <- as.character(feature$Parent)
      feature$gene_id <- transcript$gene_id[match(feature$transcript_id, transcript$transcript_id)]
      GenomicRanges::mcols(feature) <- GenomicRanges::mcols(feature)[c("source", "type", "score", "phase", "transcript_id", "gene_id")]
      return(feature)
    }) |>
    setNames(sub_tx_feature_names) |>
    GRangesList()


  gtf <- unlist(c(GRangesList(list(gene, transcript)), sub_tx_feature_list))
  names(gtf) <- NULL

  return(gtf)
}

#' Sort GTF GRanges by position and type.
#' @param gtf a GTF GRanges object.
#'
#' @export
#'
#' @importFrom GenomicRanges width start seqnames
sort_gtf <- function(gtf) {

  gtf <- gtf[order(gtf$type, decreasing=TRUE), ]
  gtf <- gtf[order(width(gtf), decreasing=TRUE), ]
  gtf <- gtf[order(start(gtf), decreasing=FALSE), ]
  gtf <- gtf[order(seqnames(gtf), decreasing=FALSE), ]

  return(gtf)
}


#' GTF type strings.
#'
#' @export
#'
ensembl_gtf_types <- c("gene", "transcript", "five_prime_utr", "exon", "CDS", "three_prime_utr")


#' Add biotype for VEP annotation.
#'
#' @param gtf a GTF GRanges.
#' @param cds_type CDS type string. Default is `"CDS"`.
#' @param coding_biotype biotype string for coding genes and transcripts. Default is `"protein_coding"`
#' @param noncoding_biotype biotype string for non-coding genes and transcripts. Default is `"ncRNA"`
#'
#' @export
#'
add_biotype_by_cds <- function(gtf, cds_type="CDS", coding_biotype="protein_coding", noncoding_biotype="ncRNA") {

  cds <- gtf[gtf$type == cds_type, ]

  gtf$gene_biotype <- noncoding_biotype
  gtf$gene_biotype[gtf$gene_id %in% cds$gene_id] <- coding_biotype
  gtf$transcript_biotype <- NA
  gtf$transcript_biotype[gtf$type != "gene"] <- noncoding_biotype
  gtf$transcript_biotype[gtf$type != "gene" & gtf$transcript_id %in% cds$transcript_id] <- coding_biotype

  return(gtf)
}
