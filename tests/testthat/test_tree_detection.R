test_that("create_disk", {
  res <- create_disk(7)
  expect_equal(res[1,1], FALSE)
})

test_that("dem_filtering", {
  data(chm_chablais3)
  # filtering with median and Gaussian smoothing
  im <- dem_filtering(chm_chablais3, nl_filter = "Median", nl_size = 3, sigmap = 0.8)
  # filtering with median filter and value-dependent Gaussian smoothing
  # (less smoothing for values between 0 and 15)
  im2 <- dem_filtering(chm_chablais3,
    nl_filter = "Median", nl_size = 3,
    sigmap = cbind(c(0.2, 0.8), c(0, 15))
  )
  expect_equal(class(im2)[1], "RasterStack")
})

test_that("maxima_detection", {
  data(chm_chablais3)
  # maxima detection
  maxi <- maxima_detection(chm_chablais3)
  expect_equal(class(maxi)[1], "RasterLayer")
})

test_that("maxima_selection", {
  data(chm_chablais3)
  # maxima detection
  maxi <- maxima_detection(chm_chablais3)
  # several maxima selection settings
  selected_maxi_hmin <- maxima_selection(maxi, chm_chablais3, hmin = 15)
  expect_equal(class(selected_maxi_hmin)[1], "RasterLayer")
  selected_maxi_dm <- maxima_selection(maxi, chm_chablais3, dm = 2.5)
  expect_equal(class(selected_maxi_dm)[1], "RasterLayer")
  selected_maxi <- maxima_selection(maxi, chm_chablais3, dm = 1, dprop = 0.1)
  expect_equal(class(selected_maxi)[1], "RasterLayer")
})

test_that("segmentation", {
  data(chm_chablais3)
  # median filter
  chm_chablais3 <- dem_filtering(chm_chablais3,
    nl_filter = "Median", nl_size = 3,
    sigmap = 0
  )$non_linear_image
  # maxima detection
  maxi <- maxima_detection(chm_chablais3)
  # maxima selection
  selected_maxi <- maxima_selection(maxi, chm_chablais3, dm = 1, dprop = 0.1)
  # segmentation
  seg_maxi <- segmentation(maxi, chm_chablais3)
  expect_equal(class(seg_maxi)[1], "RasterLayer")
  seg_selected_maxi <- segmentation(selected_maxi, chm_chablais3)
  expect_equal(class(seg_selected_maxi)[1], "RasterLayer")
})

test_that("raster_zonal_stats", {
  data(chm_chablais3)
  # median filter
  chm_chablais3 <- dem_filtering(chm_chablais3,
    nl_filter = "Median", nl_size = 3,
    sigmap = 0
  )$non_linear_image
  # maxima detection
  maxi <- maxima_detection(chm_chablais3)
  # segmentation
  seg_maxi <- segmentation(maxi, chm_chablais3)
  # compute image of maximum value in each segment
  max_in_segment <- raster_zonal_stats(seg_maxi, chm_chablais3)
  expect_equal(class(max_in_segment)[1], "RasterLayer")
})


test_that("seg_adjust", {
  data(chm_chablais3)
  # median filter
  chm_chablais3 <- dem_filtering(chm_chablais3,
    nl_filter = "Median", nl_size = 3,
    sigmap = 0
  )$non_linear_image
  # maxima detection and selection
  maxi <- maxima_detection(chm_chablais3)
  selected_maxi <- maxima_selection(maxi, chm_chablais3, dm = 1, dprop = 0.1)
  # segmentation
  seg_selected_maxi <- segmentation(selected_maxi, chm_chablais3)
  # max value in segments
  max_in_segment <- raster_zonal_stats(seg_selected_maxi , chm_chablais3)
  # segmentation modification
  seg_modif1 <- seg_adjust(seg_selected_maxi , max_in_segment,
    chm_chablais3,
    prop = 0.5
  )
  expect_equal(class(seg_modif1)[1], "RasterLayer")
  seg_modif2 <- seg_adjust(seg_selected_maxi , max_in_segment,
    chm_chablais3,
    prop = 0, min.value = 5, min.maxvalue = 10
  )
  expect_equal(class(seg_modif2)[1], "RasterLayer")
})

test_that("tree_segmentation", {
  data(chm_chablais3)
  # tree segmentation
  segments <- tree_segmentation(chm_chablais3)
  expect_equal(class(segments)[1], "RasterStack")
  segments2 <- tree_segmentation(chm_chablais3,
    nl_filter = "Median", nl_size = 3,
    sigma = cbind(c(0.2, 0.8), c(0, 15)), dmin = 0, dprop = 0, hmin = 10, crown_prop = 0.5, crown_hmin = 5
  )
  expect_equal(class(segments2)[1], "RasterStack")
})

test_that("tree_extraction", {
  data(chm_chablais3)
  # tree segmentation
  segments <- tree_segmentation(chm_chablais3)
  # tree extraction
  trees <- tree_extraction(
    segments$filled_dem, segments$local_maxima,
    segments$segments_id
  )
  expect_equal(class(trees)[1], "SpatialPointsDataFrame")
})

test_that("cimg2Raster", {
  data(chm_chablais3)
  # convert rasterLayer to cimg object
  chm_cim <- raster2Cimg(chm_chablais3)
  expect_equal(class(chm_cim)[1], "cimg")
  # apply filtering
  chm_cim_filt <- dem_filtering(
    chm_cim,
    nl_filter = "Closing",
    nl_size = 3,
    sigmap = 0
  )$non_linear_image
  # convert to RasterLayer
  chm_filt <- cimg2Raster(chm_cim_filt, chm_chablais3)
  expect_equal(class(chm_filt)[1], "RasterLayer")
})

test_that("raster2Cimg", {
  data(chm_chablais3)
  chm_cim <- raster2Cimg(chm_chablais3)
  expect_equal(class(chm_cim)[1], "cimg")
})
