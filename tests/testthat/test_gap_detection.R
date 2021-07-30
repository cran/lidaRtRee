test_that("gap_detection", {
  data(chm_chablais3)

  # fill NA values in canopy height model
  chm_chablais3[is.na(chm_chablais3)] <- 0

  # gap detection with distance larger than canopy height / 2
  gaps <- gap_detection(chm_chablais3, ratio = 2, gap_max_height = 1,
  min_gap_surface = 0)

  # gap detection with distance larger than canopy height / 2
  # and reconstruction of border areas
  gaps1 <- gap_detection(chm_chablais3,
    ratio = 2, gap_max_height = 1, min_gap_surface = 0,
    gap_reconstruct = TRUE
  )

  # gap detection without distance criterion
  gaps2 <- gap_detection(chm_chablais3, ratio = NULL, gap_max_height = 1,
  min_gap_surface = 0)

  # gap id and corresponding surface for third detection parameters
  res <- table(raster::values(gaps2$gap_id)) * raster::res(gaps2$gap_id)[1]^2
  expect_equal(res[[1]], 156)
})

test_that("edge_detection", {
  data(chm_chablais3)
  # fill NA values in canopy height model
  chm_chablais3[is.na(chm_chablais3)] <- 0
  # gap detection with distance larger than canopy height / 2
  gaps <- gap_detection(chm_chablais3,
                        ratio = 2, gap_max_height = 1, min_gap_surface = 10,
                        gap_reconstruct = TRUE
  )
  # edge detection
  edges.inside <- edge_detection(!is.na(gaps$gap_id))
  edges.outside <- edge_detection(!is.na(gaps$gap_id), inside = FALSE)
  # edge propotion
  set.seed(1234)
  res1 <- sum(raster::values(edges.inside)) / (nrow(edges.inside) * ncol(edges.inside))
  res2 <- sum(raster::values(edges.outside)) / (nrow(edges.outside) * ncol(edges.outside))
  expect_equal(res1, 0.0415219907)
  expect_equal(res2, 0.042727623)
})