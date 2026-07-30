# Classification and central-place-foraging variant ###########################

#' Classify a pollination surface into ordinal classes
#'
#' Reclassifies a pollination layer into `n` ordinal classes and returns the
#' classified raster together with the area in each class (Rahimi et al. 2021,
#' Iran case study).
#'
#' @param x A `poll_map` or `SpatRaster`.
#' @param n Number of classes (default 5).
#' @param labels Optional class labels (length `n`); defaults to
#'   very low / low / moderate / high / very high for `n = 5`.
#' @param style `"quantile"` (equal-count, default) or `"equal"` (equal-interval)
#'   breaks.
#' @param layer Optional layer name (defaults to the first).
#'
#' @return A list with the classified `map` (a categorical `SpatRaster`) and an
#'   `area` data frame (cells and hectares per class).
#' @export
poll_classify <- function(x, n = 5, labels = NULL,
                          style = c("quantile", "equal"), layer = NULL) {
  style <- match.arg(style)
  r <- if (inherits(x, "poll_map")) x$map else x
  r <- if (!is.null(layer)) r[[layer]] else r[[1]]
  if (is.null(labels)) {
    labels <- if (n == 5) c("very low", "low", "moderate", "high", "very high")
              else paste0("class", seq_len(n))
  }
  if (length(labels) != n) stop("`labels` must have length `n`.", call. = FALSE)

  vals <- terra::values(r); vals <- vals[is.finite(vals)]
  brks <- if (style == "quantile")
    stats::quantile(vals, probs = seq(0, 1, length.out = n + 1))
  else seq(min(vals), max(vals), length.out = n + 1)
  brks <- unique(brks)
  if (length(brks) < 2) stop("Not enough distinct values to classify.", call. = FALSE)

  cl <- terra::classify(r, rcl = brks, include.lowest = TRUE, brackets = TRUE)
  k  <- length(brks) - 1
  levels(cl) <- data.frame(id = seq_len(k), class = labels[seq_len(k)])
  names(cl) <- "class"

  cell_area <- prod(terra::res(r))
  ft <- terra::freq(cl)
  ft$area_ha <- ft$count * cell_area / 10000
  list(map = cl, area = ft)
}

# Linear net-gain foraging kernel (central-place foraging).
.linear_kernel <- function(x, r_max, normalize = TRUE) {
  cell <- .cell_size(x)
  rc <- max(1L, ceiling(r_max / cell))
  offs <- seq.int(-rc, rc)
  dist_m <- outer(offs, offs, function(a, b) sqrt(a^2 + b^2) * cell)
  w <- pmax(0, 1 - dist_m / r_max)          # net gain declines to 0 at r_max
  num <- terra::focal(x, w = w, fun = "sum", na.rm = TRUE, fillvalue = NA)
  if (!normalize) return(num)
  valid <- terra::ifel(is.na(x), 0, 1)
  den <- terra::focal(valid, w = w, fun = "sum", na.rm = TRUE, fillvalue = NA)
  out <- num / den
  terra::ifel(is.nan(out), NA, out)
}

#' Central-place-foraging pollination model
#'
#' A model in which foragers exploit resources only within a profitable range,
#' inspired by the foraging-theory visitation model of Olsson et al. (2015). The
#' exponential kernel of the Lonsdorf family is replaced by a linear net-gain
#' kernel, \eqn{w(d)=\max(0,\,1-d/r_{\max})}, so that visitation declines with
#' distance and stops beyond the profitable range \eqn{r_{\max}}. This is a
#' tractable approximation of central-place foraging rather than the full
#' marginal-value optimisation; the exact optimisation is planned for a future
#' release.
#'
#' @inheritParams poll_lonsdorf
#' @param r_max Profitable foraging range (map units); defaults to twice
#'   `guild$alpha`.
#'
#' @return A `poll_map` with `supply` and/or `abundance` layers.
#' @references Olsson, O. et al. (2015) *Ecological Modelling*, 316, 133-143.
#'   \doi{10.1016/j.ecolmodel.2015.08.009}
#' @seealso [poll_lonsdorf()], [poll_modified()]
#' @export
poll_cpf <- function(lulc, landcover, guild,
                     r_max = 2 * guild$alpha, layers = c("supply", "abundance")) {
  if (!inherits(lulc, "SpatRaster")) stop("`lulc` must be a SpatRaster.", call. = FALSE)
  if (!inherits(landcover, "poll_landcover")) stop("`landcover` must come from poll_landcover().", call. = FALSE)
  if (!inherits(guild, "poll_guild")) stop("`guild` must come from poll_guild().", call. = FALSE)
  layers <- match.arg(layers, c("supply", "abundance"), several.ok = TRUE)

  nesting <- .reclass(lulc, landcover$lucode, landcover$nesting)
  floral  <- .reclass(lulc, landcover$lucode, landcover$floral)

  accessible <- .linear_kernel(floral, r_max = r_max, normalize = TRUE)
  supply     <- nesting * accessible; names(supply) <- "supply"
  abundance  <- .linear_kernel(supply, r_max = r_max, normalize = TRUE); names(abundance) <- "abundance"

  out <- terra::subset(c(supply, abundance), layers)
  new_poll_map(out, meta = list(model = "cpf", alpha = guild$alpha, r_max = r_max,
                                guild = guild$name, call = match.call()))
}
