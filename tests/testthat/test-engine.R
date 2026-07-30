test_that("a uniform resource field normalises to itself", {
  r <- terra::rast(nrows = 20, ncols = 20,
                   xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
  terra::values(r) <- 1
  out <- poll_kernel(r, alpha = 300)
  v <- terra::values(out)
  expect_true(all(abs(v - 1) < 1e-6, na.rm = TRUE))
})

test_that("weights decay with distance from a point source", {
  r <- terra::rast(nrows = 21, ncols = 21,
                   xmin = 0, xmax = 2100, ymin = 0, ymax = 2100)
  terra::values(r) <- 0
  r[11, 11] <- 1
  out <- poll_kernel(r, alpha = 300, normalize = FALSE)
  m <- terra::as.matrix(out, wide = TRUE)
  expect_gt(m[11, 11], m[11, 12])   # centre > adjacent
  expect_gt(m[11, 12], m[11, 14])   # nearer > farther
  expect_equal(m[11, 12], m[12, 11], tolerance = 1e-8)  # rotational symmetry
})

test_that("NA cells are excluded from the neighbourhood", {
  r <- terra::rast(nrows = 10, ncols = 10,
                   xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  terra::values(r) <- 1
  r[1, 1] <- NA
  out <- poll_kernel(r, alpha = 200)
  # a fully-1 field with holes still normalises to 1 where data exist
  v <- terra::values(out)
  expect_true(all(abs(v - 1) < 1e-6, na.rm = TRUE))
})

test_that("non-square cells are rejected", {
  # 10 cols over 2000 m => 200 m; 10 rows over 1000 m => 100 m (non-square)
  r <- terra::rast(nrows = 10, ncols = 10,
                   xmin = 0, xmax = 2000, ymin = 0, ymax = 1000)
  terra::values(r) <- 1
  expect_error(poll_kernel(r, alpha = 100), "square")
})

test_that("alpha must be a single positive number", {
  r <- terra::rast(nrows = 5, ncols = 5)
  terra::values(r) <- 1
  expect_error(poll_kernel(r, alpha = -1), "positive")
  expect_error(poll_kernel(r, alpha = c(1, 2)), "positive")
})
