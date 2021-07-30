test_that("circle2Raster", {
  res <- circle2Raster(100, 100, 20, 1, 5)
  expect_equal(class(res)[1], "RasterLayer")
})

test_that("rasters2Cor", {
  # create raster
  r_b <- raster::raster()
  raster::extent(r_b) <- c(0, 40, 0, 40)
  raster::res(r_b) <- 1
  xy <- raster::xyFromCell(r_b, 1:length(r_b))

  # add Gaussian surface and noise
  z <- 3 * exp(-((xy[, 1] - 20)^2 + (xy[, 2] - 20)^2 / 2) / 6)
  r_b <- raster::rasterFromXYZ(cbind(xy, z))

  # create circular mask of radius 5
  z_mask <- (xy[, 1] - 20)^2 + (xy[, 2] - 20)^2 < 5^2
  r_mask <- raster::rasterFromXYZ(cbind(xy, z_mask))

  # create small raster of size 20
  r_s <- raster::crop(r_b, raster::extent(c(10, 30, 10, 30)))

  # add noise to small raster
  raster::values(r_s) <- raster::values(r_s) + rnorm(length(r_s), 0, 0.5)
  r_mask <- raster::crop(r_mask, raster::extent(c(10, 30, 10, 30)))

  # compute correlation on masked area where signal to noise ratio is lower
  set.seed(1234)
  res1 <- rasters2Cor(r_b, r_s, r_mask, small.SC = FALSE)
  expect_equal(res1, 0.92851597, tolerance = 0.1)
  # compute correlation for whole small raster
  set.seed(1234)
  res2 <- rasters2Cor(r_b, r_s, small.SC = FALSE)
  expect_equal(res2, 0.70630613, tolerance = 0.1)
  set.seed(1234)
  res3 <- rasters2Cor(r_b, r_s, small.SC = TRUE)
  expect_equal(res3, 0.51103586, tolerance = 0.1)
})

test_that("rasters_moving_cor", {
  # create raster
  r.b <- raster::raster()
  raster::extent(r.b) <- c(0, 40, 0, 40)
  raster::res(r.b) <- 1
  xy <- raster::xyFromCell(r.b, 1:length(r.b))

  # add Gaussian surfaces
  z1 <- 1.5 * exp(-((xy[, 1] - 22)^2 + (xy[, 2] - 22)^2 / 2) / 5)
  z2 <- exp(-((xy[, 1] - 20)^2 + (xy[, 2] - 22)^2 / 2) / 3)
  z3 <- 1.5 * exp(-((xy[, 1] - 17)^2 + (xy[, 2] - 17)^2 / 2) / 5)
  r.b <- raster::rasterFromXYZ(cbind(xy, z1 + z2 + z3))

  # create small raster
  r.s <- raster::crop(r.b, raster::extent(c(15, 25, 15, 25)))
  # offset raster by (-2, -2)
  raster::extent(r.s) <- c(13, 23, 13, 23)

  # compute correlations for translations inside buffer
  rr <- rasters_moving_cor(r.b, r.s, buffer = 6, step = 1)
  expect_equal(class(rr)[1], "RasterLayer")
})

test_that("raster_local_max", {
  # create raster
  r_b <- raster::raster()
  raster::extent(r_b) <- c(0, 40, 0, 40)
  raster::res(r_b) <- 1
  xy <- raster::xyFromCell(r_b, 1:length(r_b))
  # add Gaussian surfaces
  z1 <- 1.5 * exp(-((xy[, 1] - 22)^2 + (xy[, 2] - 22)^2 / 2) / 5)
  z2 <- exp(-((xy[, 1] - 20)^2 + (xy[, 2] - 22)^2 / 2) / 3)
  z3 <- 1.5 * exp(-((xy[, 1] - 17)^2 + (xy[, 2] - 17)^2 / 2) / 5)
  r_b <- raster::rasterFromXYZ(cbind(xy, z1 + z2 + z3))
  # create small raster
  r_s <- raster::crop(r_b, raster::extent(c(15, 25, 15, 25)))
  # offset raster by (-2, -2)
  raster::extent(r_s) <- c(13, 23, 13, 23)
  rr <- rasters_moving_cor(r_b, r_s, buffer = 6, step = 1)
  loc_max <- raster_local_max(rr)
  expect_equal(loc_max[[1]], 0.99)
})

test_that("coregistration", {
  # tree inventory
  trees <- data.frame(x = c(22.2, 18.3, 18.1), y = c(22.1, 22.7, 18.4),
  z = c(15, 10, 15))
  # mask of inventory area
  # empty raster with extent
  tree.mask <- circle2Raster(20, 20, 9, resolution = 1)
  # fill binary mask
  tree.mask <- raster_xy_mask(rbind(c(20, 20), c(20, 20)), c(9, 9), tree.mask,
  binary = TRUE)
  # simulate chm raster
  chm <- raster::raster()
  raster::extent(chm) <- c(0, 40, 0, 40)
  raster::res(chm) <- 1
  xy <- raster::xyFromCell(chm, 1:length(chm))
  # add Gaussian surfaces to simulate tree crowns
  z1 <- trees$z[1] * exp(-((xy[, 1] - trees$x[1])^2 + (xy[, 2] - trees$y[1])^2 / 2) * trees$z[1] / 50)
  z2 <- trees$z[2] * exp(-((xy[, 1] - trees$x[2])^2 + (xy[, 2] - trees$y[2])^2 / 2) * trees$z[2] / 50)
  z3 <- trees$z[3] * exp(-((xy[, 1] - trees$x[3])^2 + (xy[, 2] - trees$y[3])^2 / 2) * trees$z[3] / 50)
  chm <- raster::rasterFromXYZ(cbind(xy, pmax(z1, z2, z3))) #+rnorm(length(z1),0,1)))
  # translate trees
  trees$x <- trees$x + 1
  trees$y <- trees$y + 2
  coreg <- coregistration(chm, trees, mask = tree.mask, buffer = 5, step = 1, dm = 1, plot = FALSE)
  expect_equal(class(coreg$correlation_raster)[1], "RasterLayer")
  expect_equal(coreg$local_max[[1]], 0.34411264)
  coreg <- coregistration(chm, trees, mask = tree.mask, buffer = 5, step = 1, dm = 1, plot = TRUE)
})