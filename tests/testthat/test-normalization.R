test_that("counts_per_million normalizes columns", {
  counts <- matrix(
    c(100, 300, 50, 50),
    nrow = 2,
    dimnames = list(c("guide_1", "guide_2"), c("sample_a", "sample_b"))
  )

  result <- counts_per_million(counts)

  expect_equal(result["guide_1", "sample_a"], 250000)
  expect_equal(result["guide_2", "sample_a"], 750000)
  expect_equal(colSums(result), c(sample_a = 1e6, sample_b = 1e6))
})

test_that("counts_per_million rejects zero totals", {
  counts <- matrix(c(0, 0), ncol = 1)

  expect_error(counts_per_million(counts), "total count zero")
})
