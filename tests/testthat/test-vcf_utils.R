test_that("bcftools_variant_stats", {
  original <- read_vcf(test_path("fixtures", "3K-RG-5-5365122-5366701.vcf.gz"))
  a <- bcftools_variant_stats(original)
  b <- vcf_variant_stats_R(original)
  expect_equal(mean(a$MAF - b$MAF), 0, tolerance = 1e6)
  expect_equal(mean(a$Missing - b$Missing), 0, tolerance = 1e6)

  imputed <- read_vcf(test_path("fixtures", "3K-RG-5-5365122-5366701.beagle_imputed.vcf.gz"))
  a <- bcftools_variant_stats(imputed)
  b <- vcf_variant_stats_R(imputed)
  expect_equal(mean(a$MAF - b$MAF), 0, tolerance = 1e6)
  expect_equal(mean(a$Missing - b$Missing), 0, tolerance = 1e6)
})
