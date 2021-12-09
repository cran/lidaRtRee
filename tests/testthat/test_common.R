test_that("points2DSM", {
  LASfile <- system.file("extdata", "las_chablais3.laz", package="lidaRtRee")
  las_chablais3 <- lidR::readLAS(LASfile)
  # set projection
  lidR::projection(las_chablais3) <- 2154
  # create a digital surface model with first-return points, resolution 0.5 m
  dsm <- points2DSM(.las = lidR::filter_first(las_chablais3), res = 0.5)
  expect_equal(class(dsm)[1], "RasterLayer")
  # test LAS Null
  # res <- points2DSM(.las = NULL, res = 0.5)
  # res
})

test_that("points2DTM", {
  LASfile <- system.file("extdata", "las_chablais3.laz", package="lidaRtRee")
  las_chablais3 <- lidR::readLAS(LASfile)
  # set projection
  # lidR::projection(las_chablais3) <- 2154
  # create digital terrain model with points classified as ground
  dtm <- points2DTM(las_chablais3)
  expect_equal(class(dtm)[1], "RasterLayer")
  # test LAS Null
  # res <- points2DTM(.las = NULL, res = 0.5)
  # res
})

test_that("polar2Projected", {
  # create data.frame of trees with polar coordinates and diameters
  trees <- data.frame(
    x = rep(c(0, 10), each = 2),
    y = rep(c(0, 10), each = 2),
    z = rep(c(0, 2), each = 2),
    azimuth = rep(c(0, pi / 3)),
    dist = rep(c(2, 4)),
    slope = rep(c(0, pi / 6)),
    diameter.cm = c(15, 20, 25, 30)
  )
  # compute projected coordinates
  res <- polar2Projected(trees$x, trees$y, trees$z, trees$azimuth, trees$dist,
    trees$slope,
    declination = 0.03, convergence = 0.02, trees$diameter.cm / 100
  )
  mat <- data.frame(x = c(0.1037068,3.1718105,10.1062057,13.2163071), 
                    y = c(2.0724068,1.6255579,12.1223443,11.6483625), 
                    z = c(0,2,2,4), 
                    d = c(2.075000,3.5641016,2.125000,3.6141016))
  expect_equal(res, mat)
})

test_that("species_color", {
  tab.species <- species_color()
  expect_equal(tab.species[1,1], "Abies alba")
})

test_that("plot_tree_inventory", {
  # display tree inventory with CHM background
  data("chm_chablais3")
  set.seed(1234)
  raster::plot(chm_chablais3, col = gray(seq(0, 1, 1 / 255)))
  p <- plot_tree_inventory(
    tree_inventory_chablais3[, c("x", "y")],
    height = tree_inventory_chablais3$h,
    species = tree_inventory_chablais3$s,
    add = TRUE
  )
  expect_equal(class(p), "list")
  # height species and diam NULL
  raster::plot(chm_chablais3, col = gray(seq(0, 1, 1 / 255)))
  p <- plot_tree_inventory(
    tree_inventory_chablais3[, c("x", "y")],
    height = NULL,
    diam = NULL,
    species = NULL,
    add = TRUE
  )
  expect_equal(p, NULL)
  # add = FALSE
  raster::plot(chm_chablais3, col = gray(seq(0, 1, 1 / 255)))
  p <- plot_tree_inventory(
    tree_inventory_chablais3[, c("x", "y")],
    height = NULL,
    diam = NULL,
    species = NULL,
    add = FALSE
  )
  expect_equal(p, NULL)
})

test_that("raster_xy_mask", {
  # create raster
  r <- raster::raster()
  raster::extent(r) <- c(0, 40, 0, 40)
  raster::res(r) <- 1
  # xy positions
  xy <- data.frame(
    x = c(10, 20, 31.25, 15),
    y = c(10, 20, 31.25, 25)
  )
  # compute mask
  mask1 <- raster_xy_mask(xy, c(5, 8, 5, 5), r)
  mask2 <- raster_xy_mask(xy, c(5, 8, 5, 5), r, binary = FALSE)
  expect_equal(class(mask1)[1], "RasterLayer")
  expect_equal(class(mask2)[1], "RasterLayer")
})

test_that("raster_chull_mask", {
  # create raster
  r <- raster::raster()
  raster::extent(r) <- c(0, 40, 0, 40)
  raster::res(r) <- 1
  # xy positions
  xy <- data.frame(
    c(10, 20, 31.25, 15),
    c(10, 20, 31.25, 25)
  )
  # compute mask
  mask1 <- raster_chull_mask(xy, r)
  expect_equal(class(mask1)[1], "RasterLayer")
})

test_that("ellipses4Crown", {
  # compute coordinates of ellipses
  ellipses1 <- ellipses4Crown(c(0, 10), c(0, 10), c(2, 2), c(3, 4), c(2.5, 3), c(2, 3),
    id = c("A", "B")
  )
  expect_equal(ellipses1[["A"]][1], 2.5)
  # tilted ellipse
  ellipses2 <- ellipses4Crown(c(0, 10), c(0, 10), c(2, 2), c(3, 4), c(2.5, 3), c(2, 3),
    angle.offset = pi / 6
  )
  expect_equal(ellipses2[[2]][1], 12.5980762)
})

test_that("pointList2SPDF", {
  # Compute coordinates of polygons
  ellipses <- ellipses4Crown(c(0, 10), c(0, 10), c(2, 2), c(3, 4), c(2.5, 3), c(2, 3),
    id = c("A", "B")
  )
  # Convert to Spatial object
  ellipses1 <- pointList2SPDF(ellipses)
  expect_equal(class(ellipses1)[1], "SpatialPolygons")
  # Convert to Spatial object with data.frame
  ellipses2 <- pointList2SPDF(ellipses, df = data.frame(info = 1:2))
  expect_equal(class(ellipses2)[1], "SpatialPolygonsDataFrame")
})
