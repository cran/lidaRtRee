test_that("raster_metrics", {
  data(chm_chablais3)
  # raster metrics from raster
  metrics1 <- raster_metrics(chm_chablais3, res = 10)
  expect_equal(class(metrics1)[1], "RasterBrick")
  # raster metrics from data.frame
  n <- 1000
  df <- data.frame(
    x = runif(n, 0, 100), 
    y = runif(n, 0, 100), 
    z1 = runif(n, 0, 1),
    z2 = runif(n, 10, 20)
  )
  # compute raster metrics
  metrics2 <- raster_metrics(df,
    res = 10,
    fun = function(x) {
      data.frame(max.z = max(x$z1), max.sum = max(x$z1 + x$z2))
    },
    output = "data.frame"
  )
  expect_equal(metrics2[1,1], 5)
})
