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

test_that("vcf_to_field_matrix", {
  vcf <- data.frame(
    FORMAT = rep("AD:DP", 3),
    A1 = c("2,0:2", "3,1,0:4", "0,6:6"),
    A2 = c("3,0:3", "3,0,0:3", "0,6:6"),
    A3 = c("4,0:4", "0,0,5:5", "0,6:6")
  )

  ad_string_mat <- vcf_to_field_matrix(vcf, "AD")

  ad_3d_array <- ad_string_mat_to_ad_3d_array(ad_string_mat)

  expect_equal(unname(ad_string_mat[1, 1]), "2,0")
  expect_equal(ad_3d_array[2, 3, 3], 5)

  vcf$FORMAT[3] <- "AD:PL"
  expect_error(vcf_to_field_matrix(vcf))
})
