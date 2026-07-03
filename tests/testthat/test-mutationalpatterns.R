test_that("sig_contribution returns zero RelFreq (not NaN) when mut_mat is all zeros", {
  zero_mat <- matrix(0L, nrow = 96, ncol = 1, dimnames = list(paste0("ctx", seq_len(96)), "sample"))
  # signatures arg unused when sum(mut_mat) == 0
  fake_sigs <- matrix(0, nrow = 96, ncol = 3, dimnames = list(paste0("ctx", seq_len(96)), c("A", "B", "C")))
  result <- sig_contribution(mut_mat = zero_mat, signatures = fake_sigs)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("Rank", "Signature", "Contribution", "RelFreq"))
  expect_true(all(!is.nan(result$RelFreq)))
  expect_equal(result$Signature, "No Signatures found!")
  expect_equal(result$RelFreq, 0)
})

test_that("sig_plot_dbs returns placeholder ggplot list when counts are all zero", {
  dbs_contexts <- rownames(MutationalPatterns::get_known_signatures(muttype = "dbs"))
  zero_mat <- matrix(0L, nrow = length(dbs_contexts), ncol = 1,
                     dimnames = list(dbs_contexts, "sample"))
  result <- sig_plot_dbs(dbs_counts = zero_mat)

  expect_type(result, "list")
  expect_named(result, c("p_dbs_main", "p_dbs_cont"))
  expect_s3_class(result$p_dbs_main, "ggplot")
  expect_s3_class(result$p_dbs_cont, "ggplot")
})

test_that("sig_plot_mbs returns placeholder ggplot list when counts are all zero", {
  zero_counts <- c("3" = 0L, "4" = 0L, "5+" = 0L)
  result <- sig_plot_mbs(mbs_counts = zero_counts)

  expect_type(result, "list")
  expect_named(result, "p_mbs")
  expect_s3_class(result$p_mbs, "ggplot")
})

test_that("sig_mbs_table returns empty tibble when all counts are zero", {
  zero_counts <- c("3" = 0L, "4" = 0L, "5+" = 0L)
  result <- sig_mbs_table(mbs_counts = zero_counts)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("Type", "Count"))
  expect_equal(nrow(result), 0L)
})

test_that("sig_mbs_table returns correct counts ordered descending", {
  counts <- c("3" = 5L, "4" = 12L, "5+" = 2L)
  result <- sig_mbs_table(mbs_counts = counts)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("Type", "Count"))
  expect_equal(nrow(result), 3L)
  expect_equal(result$Count, c(12L, 5L, 2L))
  expect_equal(result$Type, c("4", "3", "5+"))
})

test_that("sig_mbs_table handles matrix input by summing rows", {
  mat <- matrix(c(3L, 7L, 1L), nrow = 3, ncol = 1,
                dimnames = list(c("3", "4", "5+"), "sample"))
  result <- sig_mbs_table(mbs_counts = mat)

  expect_equal(nrow(result), 3L)
  expect_equal(result$Count, c(7L, 3L, 1L))
})
