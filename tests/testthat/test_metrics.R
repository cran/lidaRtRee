test_that("clouds_metrics", {
  data(las_chablais3)
  # extract four point clouds from LAS object
  llas <- list()
  llas[["A"]] <- lidR::clip_circle(las_chablais3, 974350, 6581680, 10)
  llas[["B"]] <- lidR::clip_circle(las_chablais3, 974390, 6581680, 10)
  llas[["C"]] <- lidR::clip_circle(las_chablais3, 974350, 6581640, 10)
  # normalize point clouds
  llas <- lapply(llas, function(x) {
    lidR::normalize_height(x, lidR::tin())
  })
  # compute metrics
  res1 <- clouds_metrics(llas)
  testthat::expect_equal(res1[1,1] , 24.97)
  res2 <- clouds_metrics(llas, func = ~ mean(Z, na.rm = TRUE))
  testthat::expect_equal(res2[1,1], 8.2003752)
})

test_that("aba_metrics", {
  data(las_chablais3)
  # extract four point clouds from LAS object
  llas <- list()
  llas[["A"]] <- lidR::clip_circle(las_chablais3, 974350, 6581680, 10)
  llas[["B"]] <- lidR::clip_circle(las_chablais3, 974390, 6581680, 10)
  llas[["C"]] <- lidR::clip_circle(las_chablais3, 974350, 6581640, 10)
  # normalize point clouds
  llas <- lapply(llas, function(x) {
    lidR::normalize_height(x, lidR::tin())
  })
  # compute metrics
  res <- clouds_metrics(llas, ~ aba_metrics(
    Z, Intensity, ReturnNumber, Classification, 2, c(-Inf, 0, 2, 10, 20, +Inf)
  ))
  testthat::expect_equal(res[1,1], 24.97)
  res <- clouds_metrics(llas, ~ aba_metrics(
    Z, Intensity, ReturnNumber, Classification, 200, c(-Inf, 0, 2, 10, 20, +Inf)
  ))
  testthat::expect_equal(class(res), "data.frame")
})

test_that("std_tree_metrics", {
  # sample 50 height values
  h <- runif(50, 5, 40)
  # simulate tree data.frame
  trees <- data.frame(h = h, s = h, sp = h * 0.95, v = h * h * 0.6, vp = h * h * 0.55)
  set.seed(1234)
  res <- std_tree_metrics(trees, area_ha = 0.1)
  testthat::expect_equal(res[[1]], 21.260642, tolerance = 0.01)
})

test_that("terrain_points_metrics", {
  # sample points
  XYZ <- data.frame(x = runif(200, -10, 10), y = runif(200, -10, 10))
  XYZ$z <- 350 + 0.3 * XYZ$x + 0.1 * XYZ$y + rnorm(200, mean = 0, sd = 0.5)
  # compute terrain statistics
  set.seed(1234)
  res1 <- terrain_points_metrics(XYZ)
  testthat::expect_equal(res1[[1]], 349.7, tolerance = 0.01)
  set.seed(1234)
  res2 <- terrain_points_metrics(XYZ, centre = c(5, 5), r = 5)
  testthat::expect_equal(res2[[1]], 351.4, tolerance = 0.01)
  # with a LAS object
  data(las_chablais3)
  terrain_points <- lidR::filter_ground(las_chablais3)
  set.seed(1234)
  res3 <- terrain_points_metrics(terrain_points)
  testthat::expect_equal(res3[[1]], 1362.9, tolerance = 0.01)
  set.seed(1234)
  res4 <- terrain_points_metrics(terrain_points, centre = c(974360, 6581650), r = 10)
  testthat::expect_equal(res4[[1]], 1367.1, tolerance = 0.01)
})

test_that("clouds_tree_metrics", {
  data(las_chablais3)
  # extract three point clouds of 10 m radius from LAS object
  llas <- list()
  llas[[1]] <- lidR::clip_circle(las_chablais3, 974350, 6581680, 10)
  llas[[2]] <- lidR::clip_circle(las_chablais3, 974390, 6581680, 10)
  llas[[3]] <- lidR::clip_circle(las_chablais3, 974350, 6581640, 10)
  # normalize point clouds
  llas <- lapply(llas, function(x) {
    lidR::normalize_height(x, lidR::tin())
  })
  # compute tree metrics restricted to disks of radius 8 m.
  res1 <- clouds_tree_metrics(llas,
    cbind(c(974350, 974390, 974350), c(6581680, 6581680, 6581640)),
    8,
    res = 0.5
  )
  testthat::expect_equal(res1[1,1], 20.2, tolerance = 0.01)
})