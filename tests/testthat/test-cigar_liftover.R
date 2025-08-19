test_that("cigar_liftover", {
  aln <- minimap_homolog(
    "/home/panyq/repos/index-scripts/os/lagacy/rawdata/IRGSPb5.fa.masked.gz",
    "/home/panyq/repos/index-scripts/os/rap-db/results/genome/os.rap-db.genome.fa",
    "chr05",
    5365122,
    5366701
  )

  irgsp1_0_vcf <- read_vcf(
    "/mnt/g/Genome/3K-RG/results/vcf/NB_final_snp.vcf.gz",
    "5",
    5365122,
    5366701
  )

  build5_vcf <- read_vcf(
    "/mnt/g/Genome/NP.2023.GeSong/rawdata/1578_samples_snp.vcf.gz",
    as.character(aln@seqnames),
    aln@start,
    BiocGenerics::end(aln)
  )

  build5_vcf$NewPOS <- cigar_liftover(build5_vcf$POS, 5365122, aln)

  common_variants <- intersect(irgsp1_0_vcf$POS, build5_vcf$NewPOS)

  expect_equal(
    irgsp1_0_vcf$REF[match(common_variants, irgsp1_0_vcf$POS)],
    build5_vcf$REF[match(common_variants, build5_vcf$NewPOS)]
  )
})
