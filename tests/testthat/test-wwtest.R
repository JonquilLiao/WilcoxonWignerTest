test_that("p-value is between 0 and 1", {
  mat <- matrix(rnorm(100), 10, 10)
  mat[lower.tri(mat, diag = TRUE)] = 0
  mat = mat + t(mat)

  result <- wwtest(mat)
  expect_true(result$p.value >= 0 & result$p.value <= 1)
})
