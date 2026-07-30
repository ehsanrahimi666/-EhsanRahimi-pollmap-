# Fixtures ---------------------------------------------------------------------
make_lulc <- function() {
  r <- terra::rast(nrows = 24, ncols = 24,
                   xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
  terra::values(r) <- 211L          # arable everywhere
  r[, 1:4] <- 311L                  # broad-leaved forest strip (nesting)
  r[18:22, 14:20] <- 231L           # pasture block
  r
}

bio_tbl <- function() {
  poll_biophysical(data.frame(
    lucode = c(211L, 311L, 231L),
    nesting_cavity_availability_index = c(0.10, 0.65, 0.10),
    nesting_ground_availability_index = c(0.30, 0.35, 0.45),
    floral_resources_spring_index     = c(0.40, 0.55, 0.45),
    floral_resources_summer_index     = c(0.35, 0.30, 0.45)
  ))
}

guilds_tbl <- function() {
  poll_guilds(data.frame(
    SPECIES = c("Apis", "Osmia"),
    nesting_suitability_cavity_index = c(1.0, 1.0),
    nesting_suitability_ground_index = c(0.0, 0.1),
    foraging_activity_spring_index   = c(1.0, 1.0),
    foraging_activity_summer_index   = c(1.0, 0.2),
    alpha = c(800, 250),
    relative_abundance = c(0.9, 0.5)
  ))
}

lc_single <- function() {
  poll_landcover(data.frame(lucode = c(211L, 311L, 231L),
                            nesting = c(0.2, 0.9, 0.3),
                            floral  = c(0.6, 0.2, 0.5)))
}

# Constructors -----------------------------------------------------------------
test_that("poll_biophysical detects substrates and seasons", {
  b <- bio_tbl()
  expect_s3_class(b, "poll_biophysical")
  expect_setequal(attr(b, "substrates"), c("cavity", "ground"))
  expect_setequal(attr(b, "seasons"), c("spring", "summer"))
  expect_error(
    poll_biophysical(data.frame(lucode = 1L,
                                nesting_cavity_availability_index = 2,
                                floral_resources_spring_index = 0.5)),
    "\\[0, 1\\]"
  )
})

test_that("poll_guilds detects substrates and seasons", {
  g <- guilds_tbl()
  expect_s3_class(g, "poll_guilds")
  expect_setequal(attr(g, "substrates"), c("cavity", "ground"))
})

# Kennedy ----------------------------------------------------------------------
test_that("poll_kennedy returns supply and per-season abundance", {
  m <- poll_kennedy(make_lulc(), bio_tbl(), guilds_tbl())
  expect_s3_class(m, "poll_map")
  expect_setequal(names(poll_raster(m)),
                  c("supply", "abundance_spring", "abundance_summer"))
  v <- terra::values(poll_raster(m))
  expect_true(all(v >= 0, na.rm = TRUE))
  expect_true(all(is.finite(v[!is.na(v)])))
})

test_that("poll_kennedy errors on mismatched substrates/seasons", {
  g <- poll_guilds(data.frame(
    SPECIES = "X",
    nesting_suitability_cavity_index = 1,
    foraging_activity_spring_index = 1,
    alpha = 500, relative_abundance = 1))
  expect_error(poll_kennedy(make_lulc(), bio_tbl(), g), "match")
})

# Modified ---------------------------------------------------------------------
test_that("poll_modified returns a poll_map in range and is edge-sensitive", {
  m <- poll_modified(make_lulc(), lc_single(), poll_guild(400))
  expect_s3_class(m, "poll_map")
  expect_setequal(names(poll_raster(m)), c("supply", "abundance"))
  ab <- terra::as.matrix(poll_raster(m)[["abundance"]], wide = TRUE)
  near <- mean(ab[, 5:7],   na.rm = TRUE)   # next to forest strip
  far  <- mean(ab[, 20:22], na.rm = TRUE)
  expect_gt(near, far)
})

test_that("flight-range truncation changes the result", {
  m1 <- poll_raster(poll_modified(make_lulc(), lc_single(), poll_guild(400),
                                  max_dist = 400, layers = "abundance"))
  m2 <- poll_raster(poll_modified(make_lulc(), lc_single(), poll_guild(400),
                                  max_dist = 1200, layers = "abundance"))
  expect_false(isTRUE(all.equal(terra::values(m1), terra::values(m2))))
})

# ESTIMAP ----------------------------------------------------------------------
test_that("poll_estimap returns relative potential and honours activity", {
  lulc <- make_lulc()
  sc <- data.frame(lucode = c(211L, 311L, 231L), availability = c(0.4, 0.5, 0.7))
  base <- poll_raster(poll_estimap(lulc, sc, poll_guild(500)))
  expect_true(all(terra::values(base) >= 0, na.rm = TRUE))

  act <- base * 0 + 0.5                      # uniform 0.5 activity multiplier
  scaled <- poll_raster(poll_estimap(lulc, sc, poll_guild(500), activity = act))
  expect_equal(terra::values(scaled), terra::values(base) * 0.5, tolerance = 1e-8)
})

# MCE --------------------------------------------------------------------------
test_that("poll_mce combines criteria with normalised weights", {
  r <- terra::rast(nrows = 10, ncols = 10)
  a <- terra::setValues(r, seq_len(100))
  b <- terra::setValues(r, rev(seq_len(100)))
  out <- poll_mce(list(a = a, b = b), weights = c(a = 1, b = 1))
  expect_true(inherits(out, "SpatRaster"))
  vals <- terra::values(out)
  expect_true(all(vals >= 0 & vals <= 1, na.rm = TRUE))
  # equal weights on a and reverse(a): standardised sum is ~constant 0.5
  expect_equal(mean(vals), 0.5, tolerance = 0.02)
})

test_that("poll_mce accepts an AHP matrix and reports a consistency ratio", {
  r <- terra::rast(nrows = 8, ncols = 8)
  a <- terra::setValues(r, runif(64)); b <- terra::setValues(r, runif(64))
  # column-major fill: M[a,b]=3 means a is preferred 3x over b
  m <- matrix(c(1, 1/3, 3, 1), 2, 2, dimnames = list(c("a","b"), c("a","b")))
  out <- poll_mce(list(a = a, b = b), weights = m)
  expect_true(is.finite(attr(out, "CR")))
  w <- attr(out, "weights")
  expect_gt(unname(w["a"]), unname(w["b"]))
})
