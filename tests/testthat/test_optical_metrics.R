test_that("add_vegetation_indices", {
  df <- data.frame(nir = c(110, 150, 20), r = c(25, 50, 30), g = c(10, 60, 10))
  res <- add_vegetation_indices(df, all = TRUE)
  expect_equal(res[1,1], 110)
})
