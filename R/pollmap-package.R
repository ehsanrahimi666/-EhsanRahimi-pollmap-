#' pollmap: Spatial Modelling of Pollination Supply, Demand and Species
#' Distributions
#'
#' A raster-native framework for mapping wild-bee pollination in agricultural
#' landscapes. It brings the distance-decay family of land-use models, a
#' species-distribution pathway, and pollinator-dependence demand mapping into a
#' single, reproducible workflow built on \pkg{terra}.
#'
#' The core objects are:
#' \describe{
#'   \item{\code{\link{poll_guild}}}{a pollinator guild (foraging range).}
#'   \item{\code{\link{poll_landcover}}}{a land-cover parameter table mapping
#'     land-cover codes to nesting and floral-resource indices.}
#'   \item{\code{poll_map}}{the result of a model run: a \code{SpatRaster}
#'     of index layers plus provenance metadata.}
#' }
#'
#' The computational core is \code{\link{poll_kernel}}, a normalised
#' distance-decay (foraging) convolution reused by every land-use model.
#'
#' @keywords internal
"_PACKAGE"

# Quiet R CMD check regarding the `.` used in some terra pipelines (none yet).
NULL
