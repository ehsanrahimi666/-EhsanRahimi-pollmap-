# Landscape simulation and structure ##########################################

#' Simulate a neutral landscape
#'
#' Generates a binary categorical landscape (1 = habitat/forest, 0 = matrix) for
#' testing and for isolating the effect of landscape structure on pollination
#' (Rahimi et al. 2021). Clustering is produced natively by thresholding a
#' spatially smoothed random field, so no external dependency is required; set
#' `autocorr = 0` for a spatially random pattern.
#'
#' @param nrow,ncol Landscape dimensions in cells.
#' @param p_habitat Target proportion of habitat cells.
#' @param autocorr Spatial autocorrelation (Gaussian smoothing range, in cells);
#'   larger values give more aggregated habitat. `0` gives a random pattern.
#' @param resolution Cell size in map units.
#' @param crs Coordinate reference system.
#' @param seed Optional random seed.
#'
#' @return A single-layer categorical `SpatRaster` named `habitat`.
#' @references Rahimi, E., Barghjelveh, S. & Dong, P. (2021) *Ecological
#'   Processes*, 10, 22.
#' @seealso [poll_fragment()], [poll_metrics()]
#' @export
poll_simulate <- function(nrow = 100, ncol = 100, p_habitat = 0.3,
                          autocorr = 3, resolution = 30,
                          crs = "EPSG:2056", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol * resolution,
                   ymin = 0, ymax = nrow * resolution, crs = crs)
  noise <- terra::setValues(r, stats::rnorm(terra::ncell(r)))
  if (autocorr > 0) {
    w <- terra::focalMat(r, autocorr * resolution, "Gauss")
    noise <- terra::focal(noise, w, "sum", na.rm = TRUE)
  }
  thr <- stats::quantile(terra::values(noise), 1 - p_habitat, na.rm = TRUE)
  hab <- terra::ifel(noise >= thr, 1L, 0L)
  names(hab) <- "habitat"
  hab
}

#' Generate a fragmentation gradient
#'
#' Produces a series of neutral landscapes spanning combinations of habitat amount
#' and aggregation, for fragmentation experiments (Mitchell et al. 2015; Rahimi &
#' Jung 2024).
#'
#' @param p_seq Vector of habitat proportions.
#' @param autocorr_seq Vector of aggregation levels (see [poll_simulate()]).
#' @param nrow,ncol,resolution Landscape geometry passed to [poll_simulate()].
#' @param seed Optional random seed.
#'
#' @return A list of class `poll_fragment` with a `params` data frame and a
#'   `landscapes` list of `SpatRaster`s.
#' @references Mitchell, M.G.E. et al. (2015) *Environmental Research Letters*,
#'   10, 094014.
#' @export
poll_fragment <- function(p_seq = seq(0.1, 0.7, by = 0.2),
                          autocorr_seq = c(0, 3),
                          nrow = 100, ncol = 100, resolution = 30, seed = NULL) {
  grid <- expand.grid(p_habitat = p_seq, autocorr = autocorr_seq)
  landscapes <- Map(function(p, a)
    poll_simulate(nrow, ncol, p_habitat = p, autocorr = a,
                  resolution = resolution, seed = seed),
    grid$p_habitat, grid$autocorr)
  structure(list(params = grid, landscapes = landscapes), class = "poll_fragment")
}

#' @export
print.poll_fragment <- function(x, ...) {
  cat("<poll_fragment>: ", nrow(x$params), " landscapes\n", sep = "")
  invisible(x)
}

#' Relate pollination to landscape structure
#'
#' Summarises a pollination surface by land-cover class (mean value and class
#' proportion) and, if the \pkg{landscapemetrics} package is available, computes
#' landscape metrics for the land-cover raster (Hesselbarth et al. 2019).
#'
#' @param x A `poll_map` or `SpatRaster` of pollination values.
#' @param lulc A categorical land-cover `SpatRaster` of the same geometry.
#' @param metrics Optional character vector of landscape-metric names
#'   (`landscapemetrics` `what` argument).
#' @param level Metric level: `"class"` (default) or `"landscape"`.
#'
#' @return A list with a `class_summary` data frame (always) and, if available, a
#'   `landscapemetrics` data frame.
#' @references Hesselbarth, M.H.K. et al. (2019) *Ecography*, 42, 1648-1657.
#' @export
poll_metrics <- function(x, lulc, metrics = NULL, level = "class") {
  supply <- if (inherits(x, "poll_map")) x$map[[1]] else x[[1]]
  z <- terra::zonal(supply, lulc, fun = "mean", na.rm = TRUE)
  names(z) <- c("class", "mean_value")
  fr <- terra::freq(lulc)
  fr$proportion <- fr$count / sum(fr$count)
  summary_tbl <- merge(z, fr[, c("value", "proportion")],
                       by.x = "class", by.y = "value", all.x = TRUE)
  out <- list(class_summary = summary_tbl)
  if (requireNamespace("landscapemetrics", quietly = TRUE)) {
    out$landscapemetrics <- landscapemetrics::calculate_lsm(lulc, level = level, what = metrics)
  }
  out
}
