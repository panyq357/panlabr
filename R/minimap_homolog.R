minimap_wrapper <- function(query_fa, target_genome, out) {
  system2("minimap2", args = c("-ax", "splice", target_genome, query_fa), stdout = out)
}

sam_to_bam_wrapper <- function(sam_file, bam_file) {
  system2("samtools", args = c("view", "-b", "-o", bam_file, sam_file))
}

#' @export
setGeneric("minimap_homolog",
  function(target, query, chromosome, start, end) {
    standardGeneric("minimap_homolog")
  }
)

setMethod("minimap_homolog",
  signature(
    target = "character",
    query = "character",
    chromosome = "missing",
    start = "missing",
    end = "missing"
  ),
  function(target, query) {

    # 1. set up some tmp files.
    tmp_files <- list(
      sam = tempfile(fileext = ".sam"),
      bam = tempfile(fileext = ".bam")
    )

    # 2. minimap to target_genome.
    minimap_wrapper(query, target, tmp_files$sam)
    sam_to_bam_wrapper(tmp_files$sam, tmp_files$bam)
    aln <- GenomicAlignments::readGAlignments(tmp_files$bam)

    # 3. clean tmp files.
    lapply(tmp_files, unlink)

    return(aln)
  }
)

setMethod("minimap_homolog",
  signature(
    target = "character",
    query = "DNAStringSet",
    chromosome = "missing",
    start = "missing",
    end = "missing"
  ),
  function(target, query) {
    tmp_fasta <- tempfile(fileext = ".fasta")
    Biostrings::writeXStringSet(query, tmp_fasta)
    homolog_dna <- minimap_homolog(target, tmp_fasta)
    unlink(tmp_fasta)
    return(homolog_dna)
  }
)


setMethod("minimap_homolog",
  signature(
    target = "character",
    query = "character",
    chromosome = "character",
    start = "numeric",
    end = "numeric"
  ),
  function(target, query, chromosome, start, end) {
    tmp_fasta <- tempfile(fileext = ".fasta")
    gr <- GenomicRanges::GRanges(seqnames = chromosome, IRanges::IRanges(start, end))
    Biostrings::readDNAStringSet(query)[gr] |>
      Biostrings::writeXStringSet(tmp_fasta)
    homolog_dna <- minimap_homolog(target, tmp_fasta)
    unlink(tmp_fasta)
    return(homolog_dna)
  }
)


#' @export
aln_to_dna <- function(aln, target_genome) {
  target_gr <- GenomicRanges::GRanges(aln)
  homolog_dna <- Biostrings::readDNAStringSet(target_genome)[target_gr] |>
    setNames(paste(as.character(target_gr), aln@cigar, sep = ":"))
  return(homolog_dna)
}

