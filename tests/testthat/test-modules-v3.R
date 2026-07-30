# Fixtures ---------------------------------------------------------------------
mk_lulc2 <- function() {
  r <- terra::rast(nrows = 24, ncols = 24, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
  terra::values(r) <- 211L
  r[, 1:4] <- 311L
  r
}
lc2 <- function() poll_landcover(data.frame(lucode = c(211L, 311L),
                                            nesting = c(0.2, 0.9), floral = c(0.6, 0.2)))

# Demand -----------------------------------------------------------------------
test_that("poll_crop_dependence maps categories and numeric ratios", {
  d1 <- poll_crop_dependence(data.frame(crop = c(1L, 2L, 3L),
                                        dependence = c("modest", "high", "essential")))
  expect_s3_class(d1, "poll_crop_dependence")
  expect_equal(d1$dependence, c(0.25, 0.65, 0.95))
  d2 <- poll_crop_dependence(data.frame(crop = 1L, dependence = 0.5))
  expect_equal(d2$dependence, 0.5)
  expect_error(poll_crop_dependence(data.frame(crop = 1L, dependence = "bogus")), "Unrecognised")
})

test_that("poll_demand maps a crop raster to demand", {
  crop <- terra::rast(nrows = 10, ncols = 10)
  terra::values(crop) <- rep(c(0L, 2310L, 1110L, 2320L), length.out = 100)
  dep <- poll_crop_dependence(data.frame(crop = c(0L, 1110L, 2310L, 2320L),
                                         dependence = c(0, 0, 0.65, 0.65)))
  dm <- poll_demand(crop, dep)
  expect_s3_class(dm, "poll_map")
  v <- terra::values(poll_raster(dm))
  expect_true(all(v %in% c(0, 0.65)))
})

test_that("poll_balance flags deficits", {
  r <- terra::rast(nrows = 8, ncols = 8)
  s <- terra::setValues(r, seq_len(64))            # supply increases
  d <- terra::setValues(r, rev(seq_len(64)))       # demand decreases
  b <- poll_balance(s, d)
  expect_setequal(names(poll_raster(b)), c("balance", "deficit"))
  expect_true(all(terra::values(poll_raster(b)[["deficit"]]) %in% c(0, 1)))
})

# Analysis ---------------------------------------------------------------------
test_that("poll_scenario returns baseline, scenario and delta", {
  sc <- poll_scenario(poll_lonsdorf, mk_lulc2(), lc2(), poll_guild(400),
                      from = 311L, to = 211L)     # remove forest -> arable
  expect_s3_class(sc, "poll_scenario")
  expect_true(inherits(sc$delta, "SpatRaster"))
  # removing the only nesting habitat should reduce supply somewhere
  expect_lt(min(terra::values(sc$delta[["delta_supply"]]), na.rm = TRUE), 0)
})

test_that("poll_sensitivity summarises a stochastic run", {
  lulc <- mk_lulc2(); g <- poll_guild(400)
  run <- function() {
    lc <- poll_landcover(data.frame(lucode = c(211L, 311L),
                                    nesting = c(runif(1, 0, 0.3), runif(1, 0.7, 1)),
                                    floral  = c(runif(1, 0.4, 0.8), runif(1, 0, 0.3))))
    poll_lonsdorf(lulc, lc, g, layers = "abundance")
  }
  s <- poll_sensitivity(run, n = 8)
  expect_setequal(names(poll_raster(s)), c("mean", "sd", "cv", "q5", "q95"))
  expect_true(all(terra::values(poll_raster(s)[["sd"]]) >= 0, na.rm = TRUE))
})

test_that("beta_landcover draws a valid land-cover table", {
  priors <- data.frame(lucode = c(1L, 2L),
                       nesting_shape1 = c(2, 5), nesting_shape2 = c(5, 2),
                       floral_shape1  = c(3, 3), floral_shape2  = c(3, 3))
  f <- beta_landcover(priors)
  lc <- f()
  expect_s3_class(lc, "poll_landcover")
  expect_true(all(lc$nesting >= 0 & lc$nesting <= 1))
})

test_that("poll_validate correlates predictions with observations", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
  r <- terra::init(r, "y")                     # cell value = y coordinate
  obs <- data.frame(x = c(5, 95, 45), y = c(5, 95, 45), observed = c(5, 95, 45))
  v <- poll_validate(r, obs, coords = c("x", "y"))
  expect_true(v$spearman > 0)
  expect_equal(v$n, 3)
})

test_that("poll_compare returns consensus, cv and a correlation matrix", {
  r <- terra::rast(nrows = 8, ncols = 8); terra::values(r) <- runif(64)
  a <- r; b <- r * 2; d <- r * 0 + runif(64)
  cmp <- poll_compare(a = a, b = b, d = d)
  expect_true(inherits(cmp$consensus, "SpatRaster"))
  expect_equal(dim(cmp$correlation), c(3, 3))
  expect_equal(cmp$correlation["a", "b"], 1, tolerance = 1e-6)   # a and 2a rank-identical
})

# Landscape --------------------------------------------------------------------
test_that("poll_simulate produces a binary landscape near the target proportion", {
  s <- poll_simulate(nrow = 60, ncol = 60, p_habitat = 0.3, autocorr = 2, seed = 1)
  expect_true(inherits(s, "SpatRaster"))
  expect_setequal(sort(unique(terra::values(s))), c(0, 1))
  expect_equal(mean(terra::values(s)), 0.3, tolerance = 0.05)
})

test_that("poll_fragment builds a landscape series", {
  fr <- poll_fragment(p_seq = c(0.2, 0.4), autocorr_seq = c(0, 2),
                      nrow = 40, ncol = 40, seed = 1)
  expect_s3_class(fr, "poll_fragment")
  expect_equal(length(fr$landscapes), 4)
})

test_that("poll_metrics summarises supply by class", {
  lulc <- mk_lulc2()
  m <- poll_lonsdorf(lulc, lc2(), poll_guild(400), layers = "supply")
  out <- poll_metrics(m, lulc)
  expect_true("class_summary" %in% names(out))
  expect_true(all(c("class", "mean_value", "proportion") %in% names(out$class_summary)))
})

# SDM (glm path, self-contained) ----------------------------------------------
test_that("poll_sdm fits, projects and couples", {
  set.seed(1)
  preds <- terra::rast(nrows = 30, ncols = 30, nlyr = 2,
                       xmin = 0, xmax = 300, ymin = 0, ymax = 300)
  names(preds) <- c("v1", "v2")
  terra::values(preds) <- cbind(runif(900), runif(900))
  # presences biased to high v1
  xy <- terra::xyFromCell(preds, which(terra::values(preds[["v1"]]) > 0.6))
  occ <- data.frame(Longitude = xy[, 1], Latitude = xy[, 2])
  m <- poll_sdm(occ, preds, method = "glm", n_background = 200, seed = 1)
  expect_s3_class(m, "poll_model")
  sui <- poll_sdm_project(m, preds)
  expect_setequal(names(poll_raster(sui)), "suitability")
  vv <- terra::values(poll_raster(sui))
  expect_true(all(vv >= 0 & vv <= 1, na.rm = TRUE))
  supply <- poll_sdm_supply(list(sui, sui))
  expect_s3_class(supply, "poll_map")
})

# Classify & CPF ---------------------------------------------------------------
test_that("poll_classify returns classes and an area table", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  terra::values(r) <- seq_len(100)
  cl <- poll_classify(r, n = 5)
  expect_true(inherits(cl$map, "SpatRaster"))
  expect_true("area_ha" %in% names(cl$area))
})

test_that("poll_cpf runs and truncates at the profitable range", {
  m <- poll_cpf(mk_lulc2(), lc2(), poll_guild(300))
  expect_s3_class(m, "poll_map")
  expect_setequal(names(poll_raster(m)), c("supply", "abundance"))
  expect_true(all(terra::values(poll_raster(m)) >= 0, na.rm = TRUE))
})
