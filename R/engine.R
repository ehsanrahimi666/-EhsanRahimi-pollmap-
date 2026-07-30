# Distance-decay engine --------------------------------------------------------

#' Distance-weighted foraging convolution
#'
#' Computes, for every cell, a distance-weighted neighbourhood sum of a resource
#' surface using an exponential decay kernel \eqn{\exp(-D/\alpha)}. This is the
#' shared computational core of every distance-decay pollination model in the
#' package. With \code{normalize = TRUE} the weighted sum is divided by the sum
#' of weights over the cells that actually contributed (i.e. non-\code{NA} cells
#' and cells inside the raster extent), giving the normalised accessible-resource
#' index used by the Lonsdorf family:
#'
#' \deqn{K(x) = \frac{\sum_{x'} v(x')\, \exp(-D(x,x')/\alpha)}{\sum_{x'} \exp(-D(x,x')/\alpha)}.}
#'
#' The kernel is truncated at \code{max_dist} for efficiency; beyond roughly
#' \eqn{3\alpha} the weights are negligible.
#'
#' @param x A single-layer \code{terra} \code{SpatRaster} of resource values
#'   (for example floral resources, or a pollinator-supply surface). Cells may be
#'   \code{NA}; they are excluded from both numerator and denominator.
#' @param alpha Mean foraging distance in the raster's linear map units. A single
#'   positive number.
#' @param mask Optional single-layer \code{SpatRaster} of the same geometry whose
#'   non-zero, non-\code{NA} cells mark the valid domain over which weights are
#'   accumulated. Defaults to the non-\code{NA} cells of \code{x}.
#' @param max_dist Radius (in map units) at which the kernel is truncated.
#'   Defaults to \code{3 * alpha}.
#' @param normalize Logical; if \code{TRUE} (default) divide by the sum of
#'   weights to return a normalised index, otherwise return the raw weighted sum.
#'
#' @return A single-layer \code{SpatRaster}.
#'
#' @details
#' The current implementation uses \code{\link[terra]{focal}} with an explicit
#' weight matrix, which is exact within the truncation radius and backed by
#' compiled code. A fast-Fourier-transform path for very large foraging ranges is
#' planned for a future release and will be selectable without changing this
#' interface.
#'
#' @references
#' Lonsdorf, E., Kremen, C., Ricketts, T., Winfree, R., Williams, N. &
#' Greenleaf, S. (2009) Modelling pollination services across agricultural
#' landscapes. \emph{Annals of Botany}, 103, 1589-1600.
#' \doi{10.1093/aob/mcp069}
#'
#' @examples
#' library(terra)
#' r <- rast(nrows = 21, ncols = 21, xmin = 0, xmax = 2100, ymin = 0, ymax = 2100)
#' values(r) <- 0
#' r[11, 11] <- 1
#' k <- poll_kernel(r, alpha = 300, normalize = FALSE)
#' # values decay with distance from the central source
#' @export
poll_kernel <- function(x, alpha, mask = NULL,
                        max_dist = 3 * alpha, normalize = TRUE) {
  if (!inherits(x, "SpatRaster")) {
    stop("`x` must be a terra SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(x) != 1L) {
    stop("`x` must have a single layer.", call. = FALSE)
  }
  .check_alpha(alpha)
  cell <- .cell_size(x)

  # Exponential-decay weight matrix, truncated at max_dist.
  r_cells <- max(1L, ceiling(max_dist / cell))
  offs    <- seq.int(-r_cells, r_cells)
  dist_m  <- outer(offs, offs, function(a, b) sqrt(a^2 + b^2) * cell)
  w       <- exp(-dist_m / alpha)
  w[dist_m > max_dist] <- 0

  # Numerator: distance-weighted sum of resource values.
  num <- terra::focal(x, w = w, fun = "sum", na.rm = TRUE, fillvalue = NA)

  if (!normalize) {
    names(num) <- names(x)
    return(num)
  }

  # Denominator: sum of weights over the valid domain.
  if (is.null(mask)) {
    valid <- terra::ifel(is.na(x), 0, 1)
  } else {
    if (!inherits(mask, "SpatRaster")) {
      stop("`mask` must be a terra SpatRaster.", call. = FALSE)
    }
    valid <- terra::ifel(is.na(mask) | mask == 0, 0, 1)
  }
  den <- terra::focal(valid, w = w, fun = "sum", na.rm = TRUE, fillvalue = NA)

  out <- num / den                       # cells with zero weight become NaN
  out <- terra::ifel(is.nan(out), NA, out)
  names(out) <- names(x)
  out
}
