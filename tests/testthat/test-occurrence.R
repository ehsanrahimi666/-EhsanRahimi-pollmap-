test_that("poll_thin reduces clustering and respects the cap", {
  set.seed(1)
  occ <- data.frame(Longitude = c(rnorm(300, 8.5, 0.02), runif(40, 6, 10)),
                    Latitude  = c(rnorm(300, 47.4, 0.02), runif(40, 46, 47.7)))
  th <- poll_thin(occ, dist = 5)
  expect_lt(nrow(th), nrow(occ))
  expect_true(all(c("Longitude", "Latitude") %in% names(th)))
  expect_equal(unname(attr(th, "thin")["input"]), nrow(occ))
  th2 <- poll_thin(occ, dist = 5, max_n = 10, seed = 1)
  expect_lte(nrow(th2), 10)
})

test_that("poll_thin validates its inputs", {
  occ <- data.frame(x = 1, y = 2)
  expect_error(poll_thin(occ), "coordinate columns")
  occ2 <- data.frame(Longitude = 8, Latitude = 47)
  expect_error(poll_thin(occ2, dist = -1), "positive")
})

test_that("larger thinning distance retains fewer records", {
  set.seed(2)
  occ <- data.frame(Longitude = runif(400, 7, 9), Latitude = runif(400, 46.5, 47.5))
  expect_gte(nrow(poll_thin(occ, dist = 2)), nrow(poll_thin(occ, dist = 20)))
})

test_that("poll_sdm_eval returns AUC and TSS", {
  skip_on_cran()
  library(terra)
  set.seed(1)
  preds <- rast(nrows = 30, ncols = 30, nlyr = 2,
                xmin = 0, xmax = 300, ymin = 0, ymax = 300)
  names(preds) <- c("v1", "v2")
  values(preds) <- cbind(runif(900), runif(900))
  xy  <- xyFromCell(preds, which(values(preds[["v1"]]) > 0.6))
  occ <- data.frame(Longitude = xy[, 1], Latitude = xy[, 2])
  ev <- poll_sdm_eval(occ, preds, k = 3, block = FALSE, n_background = 200, seed = 1)
  expect_s3_class(ev, "poll_eval")
  expect_true(is.finite(ev$mean_auc))
  expect_true(ev$mean_auc >= 0 && ev$mean_auc <= 1)
  expect_equal(nrow(ev$folds), 3)
})

test_that("internal AUC is correct on a separable case", {
  expect_equal(pollmap:::.auc(c(1,1,0,0), c(0.9,0.8,0.2,0.1)), 1)
  expect_equal(pollmap:::.auc(c(1,1,0,0), c(0.1,0.2,0.8,0.9)), 0)
})
