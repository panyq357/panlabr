test_that("beagle_wrapper", {
  original_path <- test_path("fixtures", "3K-RG-5-5365122-5366701.vcf.gz")
  original <- read_vcf(original_path)
  imputed <- read_vcf(test_path("fixtures", "3K-RG-5-5365122-5366701.beagle_imputed.vcf.gz"))

  tmp_log <- file.path(tempdir(), "test-beagle_wrapper.log")
  tmp_vcf_gz <- file.path(tempdir(), "test-beagle_wrapper.vcf.gz")
  tmp_out <- file.path(tempdir(), "test-beagle_wrapper")

  beagle_wrapper(original_path, tmp_out)
  new_imputed <- read_vcf(tmp_vcf_gz)

  unlink(tmp_log)
  unlink(tmp_vcf_gz)

  expect_true(all(imputed[-(1:9)] == new_imputed[-(1:9)]))
})

test_that("vcf_to_hap_mat", {

  original <- read_vcf(test_path("fixtures", "3K-RG-5-5365122-5366701.vcf.gz"))
  expect_error(vcf_to_hap_mat(original))

  hap_mat <- readRDS(test_path("fixtures", "3K-RG-5-5365122-5366701.beagle_imputed.hap_mat.rds"))
  imputed <- read_vcf(test_path("fixtures", "3K-RG-5-5365122-5366701.beagle_imputed.vcf.gz"))
  new_hap_mat <- vcf_to_hap_mat(imputed)
  expect_true(all(hap_mat == new_hap_mat))
})

test_that("hap_wide_remove_consensus", {
  hap_mat <- readRDS(test_path("fixtures", "3K-RG-5-5365122-5366701.beagle_imputed.hap_mat.rds"))
  hap_long <- hap_mat_to_hap_long(hap_mat, group_info_k9)
  hap_wide <- hap_long_to_hap_wide(hap_long)
  expect_error(hap_wide_remove_consensus(subset(hap_wide, Count > 10)))

  hap_wide_subset <- hap_wide[hap_wide$Count > 10, ]
  expect_equal(
    hap_wide_remove_consensus(hap_wide_subset)$Seq,
    c("GGTGC", "GATGC", "GGCAG", "GGCAC", "GGTAC", "AGCAC")
  )
})
