test_that("aba_build_model", {
  data(quatre_montagnes)
  set.seed(1234)
  model_aba <- aba_build_model(quatre_montagnes$G_m2_ha,
    quatre_montagnes[, 9:76],
    transform = "boxcox",
    nmax = 3
  )
  expect_equal(model_aba$model$coefficients[[1]], 1.7867643, tolerance = 0.0000001)
  expect_equal(model_aba$model$coefficients[[2]], 0.00053183605, tolerance = 0.0000001)
  expect_equal(model_aba$model$coefficients[[3]], 0.5452766, tolerance = 0.0000001)
  expect_equal(model_aba$model$coefficients[[4]], 0.024203673, tolerance = 0.0000001)
  expect_equal(model_aba$stats[[1]], 96)
  expect_equal(model_aba$stats[[3]], 0.68897522, tolerance = 0.0000001)
  # message("Negative or NA observations of dependent variable have been removed")
  qmg <- quatre_montagnes$G_m2_ha
  qmg[3] <- NA
  model_aba <- aba_build_model(
    qmg,
    quatre_montagnes[, 9:76],
    transform = "boxcox",
    nmax = 3
  )
  expect_equal(model_aba$stats[[1]], 95)
  # transform = "none"
  model_aba <- aba_build_model(
    quatre_montagnes$G_m2_ha,
    quatre_montagnes[, 9:76],
    transform = "none",
    nmax = 3
  )
  expect_equal(model_aba, NULL)
  # transform = "log"
  model_aba <- aba_build_model(
    quatre_montagnes$G_m2_ha,
    quatre_montagnes[, 9:76],
    transform = "log",
    nmax = 3
  )
  expect_equal(model_aba$stats[[4]], "log")
  # threshold = c(1, 5)
  model_aba <- aba_build_model(
    quatre_montagnes$G_m2_ha,
    quatre_montagnes[, 9:76],
    transform = "boxcox",
    nmax = 3,
    threshold = c(1, 5)
  )
  expect_equal(model_aba$model$rank, 4)
  # xy = data.frame(898205.7, 6455708)
  model_aba <- aba_build_model(
    quatre_montagnes$G_m2_ha,
    quatre_montagnes[, 9:76],
    transform = "boxcox",
    nmax = 3,
    xy = data.frame(898205.7, 6455708)
  )
  expect_equal(model_aba$model$rank, 4)
})

test_that("boxcox_tr", {
  x <- 1:10
  res <- boxcox_tr(x, -2)
  expect_equal(res[1], 0)
  # count_reg > 0
  x <- raster::raster(system.file("external/test.grd", package = "raster"))
  x[1, 1] <- -3
  res <- boxcox_tr(x, 0)
  expect_equal(class(res)[1], "RasterLayer")
  # lambda
  x <- c(1:10, NA)
  res <- boxcox_tr(x, 0)
  expect_equal(res[1], 0)
})

test_that("boxcox_itr", {
  x <- 1:10
  res <- boxcox_itr(x, -2)
  expect_equal(res[1], NaN)
  # count_reg > 0
  x <- raster::raster(system.file("external/test.grd", package = "raster"))
  x[1, 1] <- -3
  res <- boxcox_itr(x, 0)
  expect_equal(class(res)[1], "RasterLayer")
  # lambda
  x <- c(1:10, NA)
  res <- boxcox_itr(x, 0)
  expect_equal(res[1], 2.7182818)
})

test_that("boxcox_itr_bias_cor", {
  x <- 1:10
  set.seed(1234)
  res <- boxcox_itr_bias_cor(x, 0.3, 0)
  expect_equal(res[1], 2.39779016)
  # lambda = 0
  res <- boxcox_itr_bias_cor(x, 0, 0)
  expect_equal(res[1], 2.7182818)
  # raster
  x <- raster::raster(system.file("external/test.grd", package = "raster"))
  x[1, 1] <- -3
  res <- boxcox_itr_bias_cor(x, 0.3, 0)
  expect_equal(class(res)[1], "RasterLayer")
})

test_that("aba_combine_strata", {
  # load Quatre Montagnes dataset
  data(quatre_montagnes)
  # initialize list of models
  model_aba_stratified <- list()
  # calibrate basal area prediction model for each stratum
  for (i in levels(quatre_montagnes$stratum)) {
    subsample <- which(quatre_montagnes$stratum == i)
    model_aba_stratified[[i]] <-
      aba_build_model(quatre_montagnes[subsample, "G_m2_ha"],
        quatre_montagnes[subsample, 9:76],
        transform = "boxcox", nmax = 4,
        xy = quatre_montagnes[subsample, c("X", "Y")]
      )
  }
  # combine models in single object
  set.seed(1234)
  model_aba_stratified <- aba_combine_strata(
    model_aba_stratified,
    quatre_montagnes$plotId
  )
  expect_equal(model_aba_stratified$model$private$coefficients[[1]], 1.56456992)
  expect_equal(model_aba_stratified$model$private$coefficients[[2]], -0.00280186505)
  expect_equal(model_aba_stratified$model$private$coefficients[[3]], 0.03882846)
  expect_equal(model_aba_stratified$model$private$coefficients[[4]], 0.56584811)
  expect_equal(model_aba_stratified$stats[[1,1]], 32)
  expect_equal(model_aba_stratified$stats[[1,3]], 0.62842778)
})

test_that("aba_modeler_plot", {
  # load Quatre Montagnes dataset 
  data(quatre_montagnes)
  # build ABA model for basal area, with all metrics as predictors
  model_aba <- aba_build_model(
    quatre_montagnes$G_m2_ha, 
    quatre_montagnes[, 9:76],
    transform = "boxcox", 
    nmax = 3
  )
  # plot field values VS predictions in cross-validation
  aba_plot(model_aba, main = "Basal area")
  # plot field values VS predictions in cross-validation with text
  aba_plot(model_aba, main = "Basal area", disp_text = TRUE)
  #### stratified -------------------------------------------------------------
  # load Quatre Montagnes dataset
  data(quatre_montagnes)
  # initialize list of models
  model_aba_stratified <- list()
  # calibrate basal area prediction model for each stratum
  for (i in levels(quatre_montagnes$stratum)) {
    subsample <- which(quatre_montagnes$stratum == i)
    model_aba_stratified[[i]] <-
      aba_build_model(quatre_montagnes[subsample, "G_m2_ha"],
                  quatre_montagnes[subsample, 9:76],
                  transform = "boxcox", nmax = 4,
                  xy = quatre_montagnes[subsample, c("X", "Y")]
      )
  }
  # combine models in single object
  set.seed(1234)
  model_aba_stratified <- aba_combine_strata(
    model_aba_stratified,
    quatre_montagnes$plotId
  )
  # plot field values VS predictions in cross-validation
  p <- aba_plot(model_aba_stratified, main = "Basal area stratified")
  expect_equal(p, NULL)
})

test_that("aba_predict", {
  # load data
  data(quatre_montagnes)
  # build model
  aba_model <- aba_build_model(
    quatre_montagnes$G_m2_ha, 
    quatre_montagnes[, 9:76],
    transform = "boxcox"
  )
  # build example raster to apply model
  quatre_montagnes$X <- rep(1:8, 12)
  quatre_montagnes$Y <- rep(1:12, each = 8)
  metrics_map <- raster::rasterFromXYZ(quatre_montagnes[, c(2, 3, 9:76)])
  predict_map <- lidaRtRee::aba_predict(aba_model, metrics_map)
  predict_map_error <- lidaRtRee::aba_predict(aba_model, metrics_map, add_error = TRUE)
  # log
  aba_model <- aba_build_model(
    quatre_montagnes$G_m2_ha, 
    quatre_montagnes[, 9:76],
    transform = "log"
  )
  predict_map <- lidaRtRee::aba_predict(aba_model, metrics_map)
  # log with add_error
  aba_model <- aba_build_model(
    quatre_montagnes$G_m2_ha, 
    quatre_montagnes[, 9:76],
    transform = "log"
  )
  predict_map <- lidaRtRee::aba_predict(aba_model, metrics_map, add_error = TRUE)
  # none
  # aba_model <- aba_build_model(
  #   quatre_montagnes$G_m2_ha, 
  #   quatre_montagnes[, 9:76],
  #   transform = "none",
  #   nmax = 3
  # )
  # predict_map <- lidaRtRee::aba_predict(aba_model, metrics_map)
})

test_that("clean_raster", {
  # load data
  data(quatre_montagnes)
  # build model
  aba_model <- lidaRtRee::aba_build_model(quatre_montagnes$G_m2_ha,
    quatre_montagnes[, 9:76],
    transform = "boxcox"
  )
  # build example raster to apply model
  quatre_montagnes$X <- rep(1:8, 12)
  quatre_montagnes$Y <- rep(1:12, each = 8)
  metrics_map <- raster::rasterFromXYZ(quatre_montagnes[, c(2, 3, 9:76)])
  predict_map <- lidaRtRee::aba_predict(aba_model, metrics_map)
  # create raster mask
  mask <- predict_map
  # set values to 1 or NA
  raster::values(mask) <- rep(c(1, 1, NA), each = 32)
  # apply thresholds and mask
  predict_map_clean <- clean_raster(predict_map, c(40, 70), mask)
  expect_equal(class(predict_map_clean)[1], "RasterLayer")
})

test_that("aba_inference", {
  # load data
  data(quatre_montagnes)
  # build model
  aba_model <- lidaRtRee::aba_build_model(quatre_montagnes$G_m2_ha,
                                      quatre_montagnes[, 9:76],
                                      transform = "boxcox"
  )
  # build example raster to apply model
  quatre_montagnes$X <- rep(1:8, 12)
  quatre_montagnes$Y <- rep(1:12, each = 8)
  metrics_map <- raster::rasterFromXYZ(quatre_montagnes[, c(2, 3, 9:76)])
  predict_map <- lidaRtRee::aba_predict(aba_model, metrics_map)
  # create raster mask
  mask <- predict_map
  # set values to 1 or NA
  raster::values(mask) <- rep(c(1, 1, NA), each = 32)
  # apply thresholds and mask
  predict_map_clean <- clean_raster(predict_map, c(40, 70), mask)
  res <- aba_inference(aba_model, predict_map, "SRS")
  expect_equal(res$mean, 40.174538)
  res <- aba_inference(aba_model, predict_map, "ED")
  expect_equal(res$mean, 40.146629)
  # res <- aba_inference(aba_model, predict_map, "D")
  # expect_equal(res$mean, 40.174538)
  # res <- aba_inference(aba_model, predict_map, "STR")
  # expect_equal(res$mean, 40.174538)
  res <- aba_inference(aba_model, predict_map, "SYNT")
  expect_equal(res$mean, 40.16471)
  # aba_inference(
  #   modell = aba_model, 
  #   r.predictions = predict_map, 
  #   type = "SRS", 
  #   r.mask = mask)
})
