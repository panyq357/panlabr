test_that("genbank", {
  expect_equal(
    make_pos_string(c(1, 201, 301), c(100, 250, 450), c("-", "-", "-")),
    "join(complement(1..100),complement(201..250),complement(301..450))"
  )

  expect_equal(
    make_genbank_string(
      genome = "/home/panyq/repos/index-scripts/os/rap-db/results/genome/os.rap-db.genome.fa",
      gtf = "/home/panyq/repos/index-scripts/os/rap-db/results/gtf/os.rap-db.make.gtf",
      chromosome = "chr02",
      start = 28863274,
      end = 28866997
    ),
    readLines(
      test_path("fixtures", "rap-db-chr02-28863274-28866997.gb")
    ) |> paste(collapse = "\n")
  )
})
