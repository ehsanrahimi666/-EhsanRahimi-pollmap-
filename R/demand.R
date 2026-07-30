# Pollination demand and supply-demand balance ################################

# Default Klein et al. (2007) dependence ratios by category.
.klein_ratios <- c(none = 0, independent = 0, little = 0.05, modest = 0.25,
                   moderate = 0.25, high = 0.65, great = 0.65, essential = 0.95)

#' Build a crop pollinator-dependence table
#'
#' Maps crop classes to a pollinator-dependence ratio in \[0, 1\]. Dependence may
#' be supplied directly as a numeric ratio or as a Klein et al. (2007) category
#' (`none`, `little`, `modest`, `high`, `essential`, and common synonyms), which
#' is converted to a ratio (0, 0.05, 0.25, 0.65, 0.95).
#'
#' @param data A data frame with a crop-code column and a dependence column.
#' @param crop Name of the crop-code column (integer codes matching a crop
#'   raster).
#' @param dependence Name of the dependence column (numeric ratio or category
#'   string).
#' @param ratios Named vector mapping categories to ratios; defaults to the
#'   Klein et al. (2007) scheme.
#'
#' @return An object of class `poll_crop_dependence` (a validated data frame with
#'   columns `crop` and `dependence`).
#' @references Klein, A.-M. et al. (2007) *Proceedings of the Royal Society B*,
#'   274, 303-313. \doi{10.1098/rspb.2006.3721}
#' @seealso [poll_demand()]
#' @examples
#' dep <- poll_crop_dependence(data.frame(crop = c(1L, 2L, 3L),
#'          dependence = c("modest", "high", "essential")))
#' dep
#' @export
poll_crop_dependence <- function(data, crop = "crop", dependence = "dependence",
                                 ratios = .klein_ratios) {
  data <- as.data.frame(data)
  if (!all(c(crop, dependence) %in% names(data)))
    stop("`data` must contain columns `", crop, "` and `", dependence, "`.", call. = FALSE)
  dep <- data[[dependence]]
  if (is.numeric(dep)) {
    ratio <- dep
  } else {
    key <- tolower(trimws(as.character(dep)))
    ratio <- unname(ratios[key])
    if (anyNA(ratio))
      stop("Unrecognised dependence categories: ",
           paste(unique(dep[is.na(ratio)]), collapse = ", "),
           ". Known: ", paste(names(ratios), collapse = ", "), ".", call. = FALSE)
  }
  if (any(ratio < 0 | ratio > 1, na.rm = TRUE))
    stop("Dependence ratios must lie in [0, 1].", call. = FALSE)
  out <- data.frame(crop = as.integer(data[[crop]]), dependence = as.numeric(ratio))
  structure(out, class = c("poll_crop_dependence", "data.frame"))
}

#' Map pollination demand
#'
#' Converts a crop-type raster into a relative pollination-demand surface by
#' assigning each crop its pollinator-dependence ratio (Schulp et al. 2014;
#' Picanco et al. 2017). Non-crop or unmatched cells become 0.
#'
#' @param crop A crop-type `SpatRaster` of integer codes.
#' @param dependence A [poll_crop_dependence()] table.
#'
#' @return A `poll_map` with a single `demand` layer in \[0, 1\].
#' @references Schulp, C.J.E., Lautenbach, S. & Verburg, P.H. (2014)
#'   *Ecological Indicators*, 36, 131-141.
#' @seealso [poll_crop_dependence()], [poll_balance()]
#' @examples
#' library(terra)
#' crop <- rast(nrows = 20, ncols = 20)
#' values(crop) <- rep(c(0L, 2310L, 1110L), length.out = 400)
#' dep <- poll_crop_dependence(data.frame(crop = c(0L, 1110L, 2310L),
#'                                        dependence = c(0, 0, 0.65)))
#' d <- poll_demand(crop, dep)
#' plot(d)
#' @export
poll_demand <- function(crop, dependence) {
  if (!inherits(crop, "SpatRaster")) stop("`crop` must be a SpatRaster.", call. = FALSE)
  if (!inherits(dependence, "poll_crop_dependence"))
    stop("`dependence` must come from poll_crop_dependence().", call. = FALSE)
  d <- terra::subst(crop, from = dependence$crop, to = dependence$dependence, others = 0)
  names(d) <- "demand"
  new_poll_map(d, meta = list(model = "demand", call = match.call()))
}

#' Match pollination supply and demand
#'
#' Couples a supply surface (from any land-use or SDM model) with a demand surface
#' to locate matches, mismatches and deficits (Schulp et al. 2014; Fernandes et
#' al. 2020). Both are placed on a common scale, then compared by difference or
#' ratio.
#'
#' @param supply A supply `poll_map` or `SpatRaster` (single layer).
#' @param demand A demand `poll_map` or `SpatRaster` (single layer).
#' @param method `"difference"` (standardised supply minus demand, default) or
#'   `"ratio"` (supply / demand).
#' @param standardize Logical; min-max standardise supply and demand before
#'   comparison (default `TRUE`).
#'
#' @return A `poll_map` with a `balance` layer and a `deficit` layer (1 where
#'   demand exceeds supply, else 0).
#' @references Fernandes, J. et al. (2020) *Ecosystems and People*, 16, 212-229.
#' @seealso [poll_demand()]
#' @examples
#' library(terra)
#' r <- rast(nrows = 20, ncols = 20)
#' supply <- setValues(r, runif(400)); demand <- setValues(r, runif(400))
#' b <- poll_balance(supply, demand)
#' plot(b)
#' @export
poll_balance <- function(supply, demand,
                         method = c("difference", "ratio"), standardize = TRUE) {
  method <- match.arg(method)
  s <- if (inherits(supply, "poll_map")) supply$map[[1]] else supply[[1]]
  d <- if (inherits(demand, "poll_map")) demand$map[[1]] else demand[[1]]
  std <- function(r) {
    mm <- terra::minmax(r); if (mm[2] > mm[1]) (r - mm[1]) / (mm[2] - mm[1]) else r * 0
  }
  if (standardize) { s <- std(s); d <- std(d) }
  bal <- if (method == "difference") s - d else s / (d + 1e-9)
  names(bal) <- "balance"
  deficit <- terra::ifel(d > s, 1, 0); names(deficit) <- "deficit"
  new_poll_map(c(bal, deficit),
               meta = list(model = "balance", method = method, call = match.call()))
}
