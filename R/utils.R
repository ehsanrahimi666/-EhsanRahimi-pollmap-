# Internal helpers -------------------------------------------------------------

# NULL-coalescing operator.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Reclassify a land-cover SpatRaster to a continuous index using a lookup.
# Unmatched codes become NA.
.reclass <- function(x, from, to) {
  terra::subst(x, from = from, to = to, others = NA)
}

# Assert a single positive finite number.
.check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("`alpha` must be a single positive number (mean foraging distance, in map units).",
         call. = FALSE)
  }
  invisible(alpha)
}

# Assert square cells and return the (single) cell size.
.cell_size <- function(x) {
  r <- terra::res(x)
  if (abs(r[1] - r[2]) > sqrt(.Machine$double.eps)) {
    stop("pollmap requires square raster cells; got resolution ",
         paste(signif(r, 4), collapse = " x "), ".", call. = FALSE)
  }
  r[1]
}
