# Input and output classes ----------------------------------------------------

#' Define a pollinator guild
#'
#' Creates a lightweight description of a pollinator species or guild. In this
#' first release a guild is characterised solely by its mean foraging distance
#' \code{alpha}; seasonal activity, nesting-substrate preferences and relative
#' abundance (as used by the full InVEST/Kennedy formulation) will be added to
#' this object in a later version.
#'
#' @param alpha Mean foraging distance in the same linear map units as the
#'   land-cover raster (typically metres). A single positive number.
#' @param name Optional label for the guild.
#'
#' @return An object of class \code{poll_guild}.
#' @references
#' Greenleaf, S.S., Williams, N.M., Winfree, R. & Kremen, C. (2007) Bee foraging
#' ranges and their relationship to body size. \emph{Oecologia}, 153, 589-596.
#' \doi{10.1007/s00442-007-0752-9}
#' @seealso [poll_landcover()], [poll_lonsdorf()]
#' @examples
#' poll_guild(alpha = 500, name = "Apis")
#' @export
poll_guild <- function(alpha, name = "pollinator") {
  .check_alpha(alpha)
  structure(
    list(name = as.character(name)[1], alpha = as.numeric(alpha)),
    class = "poll_guild"
  )
}

#' Build a land-cover parameter table
#'
#' Maps each land-cover class to its nesting-suitability and floral-resource
#' indices. Both indices are dimensionless and must lie in \code{[0, 1]}, where 1
#' denotes maximal suitability. The schema follows the InVEST "biophysical
#' table" so that existing InVEST parameter sets can be reused; the single-season,
#' single-substrate form used here corresponds to the original Lonsdorf (2009)
#' model.
#'
#' @param data A data frame (or object coercible to one) containing at least the
#'   three named columns below.
#' @param lucode,nesting,floral Column names in \code{data} giving, respectively,
#'   the integer land-cover code, the nesting-suitability index, and the
#'   floral-resource index.
#'
#' @return An object of class \code{poll_landcover} (a validated data frame).
#' @references
#' Lonsdorf, E., Kremen, C., Ricketts, T., Winfree, R., Williams, N. &
#' Greenleaf, S. (2009) Modelling pollination services across agricultural
#' landscapes. \emph{Annals of Botany}, 103, 1589-1600.
#' \doi{10.1093/aob/mcp069}
#' @seealso [poll_guild()], [poll_lonsdorf()]
#' @examples
#' poll_landcover(
#'   data.frame(
#'     lucode  = c(1L, 2L, 3L),
#'     nesting = c(0.0, 1.0, 0.3),
#'     floral  = c(0.8, 0.1, 0.5)
#'   )
#' )
#' @export
poll_landcover <- function(data, lucode = "lucode",
                           nesting = "nesting", floral = "floral") {
  data <- as.data.frame(data)
  req  <- c(lucode, nesting, floral)
  miss <- setdiff(req, names(data))
  if (length(miss)) {
    stop("`data` is missing column(s): ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }

  tbl <- data.frame(
    lucode  = suppressWarnings(as.integer(data[[lucode]])),
    nesting = suppressWarnings(as.numeric(data[[nesting]])),
    floral  = suppressWarnings(as.numeric(data[[floral]]))
  )

  if (anyNA(tbl$lucode)) {
    stop("`", lucode, "` must be coercible to integer with no missing values.",
         call. = FALSE)
  }
  if (anyDuplicated(tbl$lucode)) {
    stop("`", lucode, "` contains duplicated codes.", call. = FALSE)
  }
  in01 <- function(v) all(v >= 0 & v <= 1, na.rm = TRUE)
  if (!in01(tbl$nesting) || !in01(tbl$floral)) {
    stop("`nesting` and `floral` indices must lie in [0, 1].", call. = FALSE)
  }

  structure(tbl, class = c("poll_landcover", "data.frame"))
}

# Constructor for the result object. Wraps a SpatRaster with provenance so that
# downstream utilities (scenarios, model comparison) can inspect how it was made.
new_poll_map <- function(raster, meta = list()) {
  if (!inherits(raster, "SpatRaster")) {
    stop("internal: `raster` must be a terra SpatRaster.", call. = FALSE)
  }
  structure(list(map = raster, meta = meta), class = "poll_map")
}

#' Extract the raster from a poll_map
#'
#' @param x A \code{poll_map} object.
#' @return The underlying \code{terra} \code{SpatRaster}.
#' @examples
#' lulc <- terra::rast(nrows = 10, ncols = 10,
#'                     xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
#' terra::values(lulc) <- 1L
#' lulc[, 1:2] <- 2L
#' lc <- poll_landcover(data.frame(lucode = c(1L, 2L),
#'                                 nesting = c(0, 1), floral = c(0.8, 0.1)))
#' m <- poll_lonsdorf(lulc, lc, poll_guild(alpha = 200))
#' poll_raster(m)
#' @export
poll_raster <- function(x) {
  if (!inherits(x, "poll_map")) {
    stop("`x` must be a <poll_map> object.", call. = FALSE)
  }
  x$map
}
