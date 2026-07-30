# Cross-cutting analysis utilities ############################################

#' Land-use-change scenario analysis
#'
#' Re-runs a supply model under a land-cover transition and returns the change in
#' pollination. The transition reclassifies land-cover codes `from` to codes `to`;
#' all other cells are unchanged (Priess et al. 2007; Rouabah et al. 2024).
#'
#' @param fun A pollmap supply-model function (e.g. [poll_lonsdorf()],
#'   [poll_kennedy()]) whose first argument is the land-cover raster.
#' @param lulc A land-cover `SpatRaster`.
#' @param ... Further arguments passed to `fun`.
#' @param from,to Integer vectors of equal length defining the code transition.
#'
#' @return A list of class `poll_scenario` with `baseline` and `scenario`
#'   `poll_map`s and a `delta` `SpatRaster` (scenario minus baseline, per layer).
#' @seealso [poll_compare()]
#' @examples
#' library(terra)
#' lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
#' values(lulc) <- 1L; lulc[, 1:3] <- 2L
#' lc <- poll_landcover(data.frame(lucode = c(1L, 2L),
#'                                 nesting = c(0.1, 0.9), floral = c(0.6, 0.2)))
#' # convert the forest strip (2) to farmland (1) and map the change
#' sc <- poll_scenario(poll_lonsdorf, lulc, lc, poll_guild(400), from = 2L, to = 1L)
#' plot(sc$delta)
#' @export
poll_scenario <- function(fun, lulc, ..., from, to) {
  if (length(from) != length(to)) stop("`from` and `to` must have equal length.", call. = FALSE)
  baseline <- fun(lulc, ...)
  lulc2    <- terra::subst(lulc, from = from, to = to)
  scenario <- fun(lulc2, ...)
  b <- if (inherits(baseline, "poll_map")) baseline$map else baseline
  s <- if (inherits(scenario, "poll_map")) scenario$map else scenario
  delta <- s - b; names(delta) <- paste0("delta_", names(b))
  structure(list(baseline = baseline, scenario = scenario, delta = delta),
            class = "poll_scenario")
}

#' @export
print.poll_scenario <- function(x, ...) {
  cat("<poll_scenario>\n  delta layers: ", paste(names(x$delta), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Build a Beta-prior land-cover sampler
#'
#' Returns a function that draws nesting and floral indices from Beta
#' distributions (e.g. the priors of Koh et al. 2016) and returns a fresh
#' [poll_landcover()] object. Intended as the `run` ingredient of
#' [poll_sensitivity()].
#'
#' @param priors A data frame with `lucode` and the Beta shape columns
#'   `nesting_shape1`, `nesting_shape2`, `floral_shape1`, `floral_shape2`.
#'
#' @return A function of no arguments returning a `poll_landcover`.
#' @references Koh, I. et al. (2016) *PNAS*, 113, 140-145.
#' @seealso [poll_sensitivity()]
#' @examples
#' priors <- data.frame(lucode = c(1L, 2L),
#'   nesting_shape1 = c(2, 8), nesting_shape2 = c(8, 2),
#'   floral_shape1  = c(5, 3), floral_shape2  = c(3, 5))
#' draw <- beta_landcover(priors)
#' draw()   # a freshly sampled poll_landcover
#' @export
beta_landcover <- function(priors) {
  priors <- as.data.frame(priors)
  req <- c("lucode", "nesting_shape1", "nesting_shape2", "floral_shape1", "floral_shape2")
  if (!all(req %in% names(priors))) stop("`priors` needs columns: ", paste(req, collapse = ", "), ".", call. = FALSE)
  function() {
    poll_landcover(data.frame(
      lucode  = priors$lucode,
      nesting = stats::rbeta(nrow(priors), priors$nesting_shape1, priors$nesting_shape2),
      floral  = stats::rbeta(nrow(priors), priors$floral_shape1,  priors$floral_shape2)
    ))
  }
}

#' Monte-Carlo sensitivity and uncertainty mapping
#'
#' Runs a model repeatedly under sampled inputs and summarises the resulting
#' distribution per cell, addressing the uncertainty in expert-derived resource
#' scores highlighted by Rouabah et al. (2024).
#'
#' @param run A function of no arguments (or of the iteration index) returning a
#'   `poll_map` or single-layer `SpatRaster` for one draw. Typically closes over a
#'   sampler such as [beta_landcover()].
#' @param n Number of draws.
#' @param layer Optional layer name to summarise when `run` returns multiple
#'   layers (defaults to the first).
#' @param probs Quantiles to report (default 5th and 95th percentiles).
#'
#' @return A `poll_map` with `mean`, `sd`, `cv` and quantile layers.
#' @seealso [beta_landcover()]
#' @examples
#' library(terra)
#' lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
#' values(lulc) <- 1L; lulc[, 1:3] <- 2L
#' priors <- data.frame(lucode = c(1L, 2L),
#'   nesting_shape1 = c(2, 8), nesting_shape2 = c(8, 2),
#'   floral_shape1  = c(5, 3), floral_shape2  = c(3, 5))
#' draw <- beta_landcover(priors)
#' unc <- poll_sensitivity(function() poll_lonsdorf(lulc, draw(), poll_guild(400)),
#'                         n = 10)
#' plot(unc)
#' @export
poll_sensitivity <- function(run, n = 50, layer = NULL, probs = c(0.05, 0.95)) {
  if (!is.function(run)) stop("`run` must be a function returning a poll_map/SpatRaster.", call. = FALSE)
  get_layer <- function(res) {
    r <- if (inherits(res, "poll_map")) res$map else res
    if (!is.null(layer)) r[[layer]] else r[[1]]
  }
  first <- get_layer(if (length(formals(run))) run(1) else run())
  stk <- first
  for (i in 2:n) {
    stk <- c(stk, get_layer(if (length(formals(run))) run(i) else run()))
  }
  m  <- terra::app(stk, mean);  names(m) <- "mean"
  sdr <- terra::app(stk, stats::sd); names(sdr) <- "sd"
  cv <- terra::ifel(m == 0, NA, sdr / m); names(cv) <- "cv"
  qs <- terra::quantile(stk, probs = probs)
  names(qs) <- paste0("q", probs * 100)
  new_poll_map(c(m, sdr, cv, qs),
               meta = list(model = "sensitivity", n = n, call = match.call()))
}

#' Validate predictions against observations
#'
#' Compares predicted pollination indices with field observations at survey
#' locations using rank correlation and error summaries (Cunningham et al. 2018).
#'
#' @param x A predicted `poll_map` or `SpatRaster`.
#' @param observed A data frame with coordinate columns and an observed-value
#'   column, or a `SpatVector` of points.
#' @param coords Length-2 character vector naming the x/y columns (data-frame
#'   input).
#' @param value Name of the observed-value column.
#' @param layer Optional predicted-layer name (defaults to the first).
#'
#' @return A list with the Spearman correlation, RMSE on standardised values, the
#'   number of points used, and a data frame of paired predicted/observed values.
#' @references Cunningham, C. et al. (2018) *Int. J. Biodiversity Science,
#'   Ecosystem Services & Management*, 14, 60-70.
#' @examples
#' library(terra)
#' r <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
#' r <- init(r, "y")
#' obs <- data.frame(x = c(10, 50, 90), y = c(10, 50, 90), observed = c(10, 50, 90))
#' poll_validate(r, obs, coords = c("x", "y"))
#' @export
poll_validate <- function(x, observed, coords = c("x", "y"),
                          value = "observed", layer = NULL) {
  r <- if (inherits(x, "poll_map")) x$map else x
  if (!is.null(layer)) r <- r[[layer]] else r <- r[[1]]
  if (inherits(observed, "SpatVector")) {
    pts <- observed; obs <- as.data.frame(observed)[[value]]
  } else {
    observed <- as.data.frame(observed)
    if (!all(c(coords, value) %in% names(observed)))
      stop("`observed` must contain columns ", paste(c(coords, value), collapse = ", "), ".", call. = FALSE)
    pts <- terra::vect(observed, geom = coords, crs = terra::crs(r))
    obs <- observed[[value]]
  }
  pred <- terra::extract(r, pts)[, 2]
  ok <- is.finite(pred) & is.finite(obs)
  pred <- pred[ok]; obs <- obs[ok]
  std <- function(v) { rg <- range(v); if (diff(rg) > 0) (v - rg[1]) / diff(rg) else v * 0 }
  rmse <- sqrt(mean((std(pred) - std(obs))^2))
  list(spearman = stats::cor(pred, obs, method = "spearman"),
       rmse = rmse, n = length(pred),
       data = data.frame(predicted = pred, observed = obs))
}

#' Compare multiple models
#'
#' Runs several models on identical inputs (supplied as their result maps),
#' standardises them, and summarises agreement: pairwise rank correlations, a
#' per-cell consensus (mean) surface and a disagreement (coefficient of variation)
#' surface (Rahimi et al. 2021).
#'
#' @param ... Two or more `poll_map`s or single-layer `SpatRaster`s, or a single
#'   named list of them.
#' @param layer Optional layer name to extract from each (defaults to the first).
#' @param standardize Logical; min-max standardise each map (default `TRUE`).
#'
#' @return A list with a `consensus` `SpatRaster`, a `cv` `SpatRaster` and a
#'   Spearman `correlation` matrix.
#' @seealso [poll_scenario()]
#' @examples
#' library(terra)
#' r <- rast(nrows = 20, ncols = 20); values(r) <- runif(400)
#' cmp <- poll_compare(model_a = r, model_b = r * 2)
#' cmp$correlation
#' @export
poll_compare <- function(..., layer = NULL, standardize = TRUE) {
  maps <- list(...)
  if (length(maps) == 1 && is.list(maps[[1]]) && !inherits(maps[[1]], "poll_map")) maps <- maps[[1]]
  if (length(maps) < 2) stop("Provide at least two maps to compare.", call. = FALSE)
  nms <- names(maps); if (is.null(nms)) nms <- paste0("model", seq_along(maps))
  pull <- function(m) {
    r <- if (inherits(m, "poll_map")) m$map else m
    r <- if (!is.null(layer)) r[[layer]] else r[[1]]
    if (standardize) { mm <- terra::minmax(r); if (mm[2] > mm[1]) (r - mm[1]) / (mm[2] - mm[1]) else r * 0 } else r
  }
  stk <- terra::rast(lapply(maps, pull)); names(stk) <- nms
  consensus <- terra::app(stk, mean); names(consensus) <- "consensus"
  sdr <- terra::app(stk, stats::sd)
  cv <- terra::ifel(consensus == 0, NA, sdr / consensus); names(cv) <- "cv"
  vals <- terra::values(stk); vals <- vals[stats::complete.cases(vals), , drop = FALSE]
  # a spatially constant map has undefined rank correlation (returns NA); the
  # zero-variance warning from cor() is expected here and is suppressed.
  cmat <- suppressWarnings(stats::cor(vals, method = "spearman"))
  list(consensus = consensus, cv = cv, correlation = cmat)
}
