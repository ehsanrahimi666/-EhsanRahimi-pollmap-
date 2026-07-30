# Land-use supply models -------------------------------------------------------

#' Pollination supply and abundance (Lonsdorf 2009 model)
#'
#' Implements the original Lonsdorf \emph{et al.} (2009) distance-decay model in
#' its single-season, single-guild form. The model proceeds in two steps:
#'
#' \enumerate{
#'   \item \strong{Pollinator supply} (the nesting-source index): the nesting
#'     suitability of each cell multiplied by the normalised, distance-weighted
#'     floral resources it can access,
#'     \deqn{G_i = N_i \cdot \frac{\sum_j F_j \exp(-D_{ij}/\alpha)}{\sum_j \exp(-D_{ij}/\alpha)}.}
#'   \item \strong{Pollinator abundance} (the visitation index): the supply
#'     projected forward onto the landscape with the same decay kernel,
#'     \deqn{P_j = \frac{\sum_i G_i \exp(-D_{ij}/\alpha)}{\sum_i \exp(-D_{ij}/\alpha)}.}
#' }
#'
#' Both outputs are relative indices in \code{[0, 1]}. This original formulation
#' spreads foragers equally with distance and does not weight destinations by
#' floral quality; the quality-weighted, central-place variant (following Olsson
#' \emph{et al.} 2015 and used in current InVEST) will be provided by
#' \code{poll_kennedy()} and \code{poll_modified()}.
#'
#' @param lulc A land-cover \code{terra} \code{SpatRaster} of integer codes, in a
#'   projected coordinate system with square cells (linear units, e.g. metres).
#' @param landcover A \code{\link{poll_landcover}} table mapping each code in
#'   \code{lulc} to nesting and floral indices.
#' @param guild A \code{\link{poll_guild}} giving the foraging distance
#'   \code{alpha}.
#' @param layers Which output layer(s) to return: any of \code{"supply"} and
#'   \code{"abundance"} (both by default).
#'
#' @return A \code{poll_map} object wrapping a \code{SpatRaster} with the
#'   requested layer(s) and provenance metadata.
#'
#' @references
#' Lonsdorf, E., Kremen, C., Ricketts, T., Winfree, R., Williams, N. &
#' Greenleaf, S. (2009) Modelling pollination services across agricultural
#' landscapes. \emph{Annals of Botany}, 103, 1589-1600.
#' \doi{10.1093/aob/mcp069}
#'
#' Rahimi, E., Barghjelveh, S. & Dong, P. (2021) Using the Lonsdorf model for
#' estimating habitat loss and fragmentation effects on pollination service.
#' \emph{Ecological Processes}, 10, 22. \doi{10.1186/s13717-021-00291-8}
#'
#' @seealso [poll_kernel()], [poll_landcover()], [poll_guild()]
#'
#' @examples
#' library(terra)
#' # A 20 x 20 landscape: farmland (code 1) with a forest nesting strip (code 2)
#' lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
#' values(lulc) <- 1L
#' lulc[, 1:3] <- 2L
#'
#' lc <- poll_landcover(
#'   data.frame(lucode  = c(1L, 2L),
#'              nesting = c(0.0, 1.0),
#'              floral  = c(0.8, 0.1))
#' )
#' g <- poll_guild(alpha = 300, name = "wild_bee")
#'
#' m <- poll_lonsdorf(lulc, lc, g)
#' m
#' # terra::plot(poll_raster(m))
#' @examples
#' library(terra)
#' lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
#' values(lulc) <- 1L; lulc[, 1:3] <- 2L    # farmland with a forest strip
#' lc <- poll_landcover(data.frame(lucode = c(1L, 2L),
#'                                 nesting = c(0.1, 0.9), floral = c(0.6, 0.2)))
#' m <- poll_lonsdorf(lulc, lc, poll_guild(alpha = 400))
#' plot(m)
#' @export
poll_lonsdorf <- function(lulc, landcover, guild,
                          layers = c("supply", "abundance")) {
  # --- input validation ---
  if (!inherits(lulc, "SpatRaster")) {
    stop("`lulc` must be a terra SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(lulc) != 1L) {
    stop("`lulc` must have a single layer.", call. = FALSE)
  }
  if (!inherits(landcover, "poll_landcover")) {
    stop("`landcover` must be built with poll_landcover().", call. = FALSE)
  }
  if (!inherits(guild, "poll_guild")) {
    stop("`guild` must be built with poll_guild().", call. = FALSE)
  }
  layers <- match.arg(layers, c("supply", "abundance"), several.ok = TRUE)
  alpha  <- guild$alpha

  # Warn about land-cover codes present in the raster but absent from the table.
  present <- terra::unique(lulc)[[1]]
  present <- present[!is.na(present)]
  unknown <- setdiff(present, landcover$lucode)
  if (length(unknown)) {
    warning("LULC code(s) with no parameters (mapped to NA): ",
            paste(unknown, collapse = ", "), ".", call. = FALSE)
  }

  # --- map land cover to resource surfaces ---
  nesting <- .reclass(lulc, landcover$lucode, landcover$nesting)
  floral  <- .reclass(lulc, landcover$lucode, landcover$floral)

  # Step 1: supply = nesting x accessible floral resources (G_i).
  accessible <- poll_kernel(floral, alpha = alpha, normalize = TRUE)
  supply     <- nesting * accessible
  names(supply) <- "supply"

  # Step 2: abundance = forward projection of supply (P_j).
  abundance     <- poll_kernel(supply, alpha = alpha, normalize = TRUE)
  names(abundance) <- "abundance"

  both <- c(supply, abundance)
  out  <- terra::subset(both, layers)

  new_poll_map(
    out,
    meta = list(
      model = "lonsdorf",
      alpha = alpha,
      guild = guild$name,
      layers = layers,
      call  = match.call()
    )
  )
}
