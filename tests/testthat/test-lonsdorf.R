# Shared fixture: farmland (1) with a forest nesting strip (2) on the left.
make_landscape <- function() {
  r <- terra::rast(nrows = 20, ncols = 20,
                   xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
  terra::values(r) <- 1L
  r[, 1:3] <- 2L
  r
}

make_landcover <- function() {
  poll_landcover(
    data.frame(lucode  = c(1L, 2L),
               nesting = c(0.0, 1.0),
               floral  = c(0.8, 0.1))
  )
}

test_that("poll_landcover validates its inputs", {
  expect_s3_class(make_landcover(), "poll_landcover")
  expect_error(
    poll_landcover(data.frame(lucode = 1L, nesting = 0.5)),
    "missing column"
  )
  expect_error(
    poll_landcover(data.frame(lucode = 1L, nesting = 2, floral = 0.5)),
    "\\[0, 1\\]"
  )
  expect_error(
    poll_landcover(data.frame(lucode = c(1L, 1L),
                              nesting = c(0, 1), floral = c(0.5, 0.5))),
    "duplicated"
  )
})

test_that("poll_guild validates alpha", {
  expect_s3_class(poll_guild(500), "poll_guild")
  expect_error(poll_guild(0), "positive")
})

test_that("poll_lonsdorf returns a well-formed poll_map", {
  m <- poll_lonsdorf(make_landscape(), make_landcover(), poll_guild(300))
  expect_s3_class(m, "poll_map")
  expect_setequal(names(poll_raster(m)), c("supply", "abundance"))
  v <- terra::values(poll_raster(m))
  expect_true(all(v >= 0 & v <= 1, na.rm = TRUE))
  expect_identical(m$meta$model, "lonsdorf")
  expect_identical(m$meta$alpha, 300)
})

test_that("layer selection works", {
  m <- poll_lonsdorf(make_landscape(), make_landcover(),
                     poll_guild(300), layers = "abundance")
  expect_identical(names(poll_raster(m)), "abundance")
})

test_that("abundance is higher near nesting habitat than far from it", {
  ab <- poll_raster(
    poll_lonsdorf(make_landscape(), make_landcover(),
                  poll_guild(300), layers = "abundance")
  )
  m <- terra::as.matrix(ab, wide = TRUE)
  near <- mean(m[, 4:6],   na.rm = TRUE)   # just right of the forest strip
  far  <- mean(m[, 17:19], na.rm = TRUE)   # far side of the field
  expect_gt(near, far)
})

test_that("unknown LULC codes raise a warning", {
  lulc <- make_landscape()
  lulc[1, 1] <- 9L   # code with no parameters
  expect_warning(
    poll_lonsdorf(lulc, make_landcover(), poll_guild(300)),
    "no parameters"
  )
})
