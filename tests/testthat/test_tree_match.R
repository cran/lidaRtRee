test_that("tree_matching", {
  # create reference and detected trees
  ref_trees <- cbind(c(1, 4, 3, 4, 2), c(1, 1, 2, 3, 4), c(15, 18, 20, 10, 11))
  def_trees <- cbind(c(2, 2, 4, 4), c(1, 3, 4, 1), c(16, 19, 9, 15))
  # match trees
  match1 <- tree_matching(ref_trees, def_trees)
  match2 <- tree_matching(ref_trees, def_trees, delta_ground = 2, h_prec = 0)
  expect_equal(match1[1,1], 1)
  expect_equal(match2[1,1], 1)
})

test_that("plot_matched", {
  # create reference and detected trees
  ref_trees <- cbind(c(1, 4, 3, 4, 2), c(1, 1.5, 2, 3, 4), c(15, 18, 20, 10, 11))
  def_trees <- cbind(c(2, 2, 4, 4), c(1, 3, 4, 1), c(16, 19, 9, 15))
  # compute matching
  match1 <- tree_matching(ref_trees, def_trees)
  match2 <- tree_matching(ref_trees, def_trees, delta_ground = 2, h_prec = 0)
  # 2D display of matching results
  p1 <- plot_matched(ref_trees, def_trees, match1, xlab = "X", ylab = "Y")
  p2 <- plot_matched(ref_trees, def_trees, match2, xlab = "X", ylab = "Y")
  expect_equal(class(p1), "list")
  expect_equal(class(p2), "list")
})

test_that("hist_detection", {
  # create reference and detected trees
  ref_trees <- cbind(c(1, 4, 3, 4, 2), c(1, 1, 2, 3, 4), c(15, 18, 20, 10, 11))
  def_trees <- cbind(c(2, 2, 4, 4), c(1, 3, 4, 1), c(16, 19, 9, 15))
  # tree matching with different buffer size
  match1 <- tree_matching(ref_trees, def_trees)
  match2 <- tree_matching(ref_trees, def_trees, delta_ground = 2, h_prec = 0)
  # corresponding number of detections
  h1 <- hist_detection(ref_trees, def_trees, match1)
  h2 <- hist_detection(ref_trees, def_trees, match2)
  expect_equal(class(h1), "list")
  expect_equal(class(h2), "list")
})

test_that("hist_stack", {
  # create reference and detected trees
  ref_trees <- cbind(c(1, 4, 3, 4, 2), c(1, 1, 2, 3, 4), c(15, 18, 20, 10, 11))
  def_trees <- cbind(c(2, 2, 4, 4), c(1, 3, 4, 1), c(16, 19, 9, 15))
  # tree matching with different buffer size
  match1 <- tree_matching(ref_trees, def_trees)
  match2 <- tree_matching(ref_trees, def_trees, delta_ground = 2, h_prec = 0)
  # hist_stack
  dummy <- list()
  lr <- ref_trees
  ld <- def_trees
  matched <- match1
  # true detections
  dummy[[1]] <- lr[matched[, 1], 3]
  # true omissions
  dummy[[2]] <- lr[-matched[, 1], 3]
  # false detections
  dummy[[3]] <- ld[-matched[, 2], 3]
  # corresponding number of detections
  # draw histogram
  h1 <- hist_stack(dummy, breaks = seq(
      from = 0, 
      to = ceiling(max(lr[, 3], ld[, 3]) / 5) * 5, 
      by = 5), 
      col = c("green", "red", "blue"))
  expect_equal(class(h1), "numeric")
  # dummy is null
  h1 <- function() stop(hist_stack(NULL, breaks = seq(
    from = 0, 
    to = ceiling(max(lr[, 3], ld[, 3]) / 5) * 5, 
    by = 5), 
    col = c("green", "red", "blue")))
  expect_error(h1(), "'x' must be a list.")
  # col is null
  h1 <- hist_stack(dummy, breaks = seq(
    from = 0, 
    to = ceiling(max(lr[, 3], ld[, 3]) / 5) * 5, 
    by = 5), 
    col = NULL)
  expect_equal(class(h1), "numeric")
  # x have a non numeric
  dum <- dummy
  dum[[1]][1] <- "A"
  h1 <- function() stop(hist_stack(
    dum, 
    breaks = seq(
      from = 0, 
      to = ceiling(max(lr[, 3], ld[, 3]) / 5) * 5, 
      by = 5), 
    col = c("green", "red", "blue")))
  expect_error(h1(), "'x' must be numeric")
})

test_that("height_regression", {
  # create tree locations and heights
  ref_trees <- cbind(c(1, 4, 3, 4, 2), c(1, 1, 2, 3, 4), c(15, 18, 20, 10, 11))
  def_trees <- cbind(c(2, 2, 4, 4), c(1, 3, 4, 1), c(16, 19, 9, 15))
  # tree matching
  match1 <- tree_matching(ref_trees, def_trees)
  # height regression
  reg <- height_regression(ref_trees, def_trees, match1,
    species = c("ABAL", "ABAL", "FASY", "FASY", "ABAL"),
    asp = 1, xlim = c(0, 21), ylim = c(0, 21)
  )
  summary(reg$lm)
  expect_equal(reg$lm$rank, 2)
  expect_equal(reg$stats$rmse, 2)
})

