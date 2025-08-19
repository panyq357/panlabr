test_that("minimap_homolog", {

  target_path <- "/home/panyq/repos/index-scripts/os/lagacy/results/genome/os.build5.genome.fa"
  query_path <- test_path("fixtures", "build5-chr05-5365128-5366707.fasta")
  query_dna <- Biostrings::readDNAStringSet(query_path)

  homolog_dna <- minimap_homolog(target_path, query_path) |> aln_to_dna(target_path)

  expect_equal(as.character(query_dna), as.character(homolog_dna))
  expect_equal(as.character(names(query_dna)), as.character(names(homolog_dna)))

  homolog_dna <- minimap_homolog(target_path, query_dna) |> aln_to_dna(target_path)

  expect_equal(as.character(query_dna), as.character(homolog_dna))
  expect_equal(as.character(names(query_dna)), as.character(names(homolog_dna)))

  homolog_dna <- minimap_homolog(
    target_path,
    "/home/panyq/repos/index-scripts/os/rap-db/results/genome/os.rap-db.genome.fa",
    "chr05",
    5365122,
    5366701
  ) |> aln_to_dna(target_path)

  expect_equal(as.character(query_dna), as.character(homolog_dna))
  expect_equal(as.character(names(query_dna)), as.character(names(homolog_dna)))
})
