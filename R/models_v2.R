# Multi-guild inputs and the seasonal / quality-weighted / ESTIMAP / MCE models #

# ---- Input constructors (InVEST-schema) -------------------------------------

#' Build a multi-substrate, multi-season biophysical table
#'
#' Validates a land-cover parameter table in the InVEST "biophysical" schema, with
#' one nesting-availability column per substrate
#' (`nesting_<substrate>_availability_index`) and one floral-resource column per
#' season (`floral_resources_<season>_index`). Substrates and seasons are detected
#' automatically from the column names. All indices must lie in \[0, 1\].
#'
#' @param data A data frame (or coercible) with an integer `lucode` column and the
#'   nesting/floral index columns described above.
#' @param lucode Name of the land-cover-code column.
#'
#' @return An object of class `poll_biophysical` with `substrates` and `seasons`
#'   attributes.
#' @references Kennedy, C.M. et al. (2013) *Ecology Letters*, 16, 584-599.
#'   \doi{10.1111/ele.12082}
#' @seealso [poll_guilds()], [poll_kennedy()]
#' @examples
#' bio <- poll_biophysical(read.csv(system.file("extdata",
#'           "clc_biophysical_table.csv", package = "pollmap")))
#' bio
#' @export
poll_biophysical <- function(data, lucode = "lucode") {
  data <- as.data.frame(data)
  if (!lucode %in% names(data)) {
    stop("`data` must contain the code column `", lucode, "`.", call. = FALSE)
  }
  nest_cols <- grep("^nesting_.+_availability_index$", names(data), value = TRUE)
  flor_cols <- grep("^floral_resources_.+_index$",     names(data), value = TRUE)
  if (!length(nest_cols)) stop("No `nesting_<substrate>_availability_index` columns found.", call. = FALSE)
  if (!length(flor_cols)) stop("No `floral_resources_<season>_index` columns found.", call. = FALSE)

  substrates <- sub("^nesting_(.+)_availability_index$", "\\1", nest_cols)
  seasons    <- sub("^floral_resources_(.+)_index$",     "\\1", flor_cols)

  idx <- unlist(data[, c(nest_cols, flor_cols)])
  if (any(idx < 0 | idx > 1, na.rm = TRUE)) {
    stop("All nesting and floral indices must lie in [0, 1].", call. = FALSE)
  }
  data[[lucode]] <- as.integer(data[[lucode]])
  if (anyNA(data[[lucode]]) || anyDuplicated(data[[lucode]])) {
    stop("`", lucode, "` must be unique integers with no missing values.", call. = FALSE)
  }
  names(data)[names(data) == lucode] <- "lucode"
  structure(data, class = c("poll_biophysical", "data.frame"),
            substrates = substrates, seasons = seasons)
}

#' Build a multi-species guild table
#'
#' Validates a pollinator guild table in the InVEST schema: one row per species,
#' with nesting-suitability columns (`nesting_suitability_<substrate>_index`),
#' foraging-activity columns (`foraging_activity_<season>_index`), an `alpha`
#' foraging distance (map units) and a `relative_abundance`.
#'
#' @param data A data frame (or coercible) with a species column plus the columns
#'   above.
#' @param species Name of the species column.
#'
#' @return An object of class `poll_guilds` with `substrates` and `seasons`
#'   attributes.
#' @references Greenleaf, S.S. et al. (2007) *Oecologia*, 153, 589-596.
#'   \doi{10.1007/s00442-007-0752-9}
#' @seealso [poll_biophysical()], [poll_kennedy()]
#' @examples
#' guilds <- poll_guilds(read.csv(system.file("extdata",
#'              "swiss_guild_table.csv", package = "pollmap")))
#' guilds
#' @export
poll_guilds <- function(data, species = "SPECIES") {
  data <- as.data.frame(data)
  if (!species %in% names(data)) {
    stop("`data` must contain the species column `", species, "`.", call. = FALSE)
  }
  for (col in c("alpha", "relative_abundance")) {
    if (!col %in% names(data)) stop("Missing required column `", col, "`.", call. = FALSE)
  }
  ns_cols <- grep("^nesting_suitability_.+_index$", names(data), value = TRUE)
  fa_cols <- grep("^foraging_activity_.+_index$",   names(data), value = TRUE)
  if (!length(ns_cols) || !length(fa_cols)) {
    stop("Guild table needs `nesting_suitability_<substrate>_index` and ",
         "`foraging_activity_<season>_index` columns.", call. = FALSE)
  }
  if (any(data$alpha <= 0, na.rm = TRUE)) stop("`alpha` must be positive.", call. = FALSE)

  substrates <- sub("^nesting_suitability_(.+)_index$", "\\1", ns_cols)
  seasons    <- sub("^foraging_activity_(.+)_index$",   "\\1", fa_cols)
  names(data)[names(data) == species] <- "SPECIES"
  structure(data, class = c("poll_guilds", "data.frame"),
            substrates = substrates, seasons = seasons)
}

#' @export
print.poll_biophysical <- function(x, ...) {
  cat("<poll_biophysical>: ", nrow(x), " classes | substrates: ",
      paste(attr(x, "substrates"), collapse = ", "), " | seasons: ",
      paste(attr(x, "seasons"), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
print.poll_guilds <- function(x, ...) {
  cat("<poll_guilds>: ", nrow(x), " species (", paste(x$SPECIES, collapse = ", "),
      ")\n", sep = "")
  invisible(x)
}

# ---- Kennedy / InVEST seasonal, guild-based model ---------------------------

#' Seasonal, guild-based pollination model (Kennedy / InVEST)
#'
#' Implements the multi-substrate, multi-season, multi-guild model used in InVEST
#' (Kennedy et al. 2013) with the central-place-foraging quality weighting of
#' Olsson et al. (2015). For each species \eqn{s}, habitat nesting suitability,
#' accessible floral resources and pollinator supply are
#' \deqn{HN(x,s)=\max_n N(l(x),n)\,ns(s,n),\quad
#'       FR(x,s)=\mathrm{kernel}\big(\textstyle\sum_j F_j\,fa(s,j)\big),\quad
#'       PS(x,s)=FR(x,s)\,HN(x,s)\,sa(s),}
#' and seasonal on-site abundance is the supply projected forward and weighted by
#' local floral quality,
#' \deqn{PA(x,s,j)=\frac{F_j(x)\,fa(s,j)}{FR(x,s)}\;\mathrm{kernel}\big(PS(\cdot,s)\big).}
#' Outputs are aggregated over species. Yield and managed-bee functions are not
#' computed (this is an ecological model).
#'
#' @param lulc A land-cover `SpatRaster` (projected, square cells, metres).
#' @param biophysical A [poll_biophysical()] table.
#' @param guilds A [poll_guilds()] table whose substrates and seasons match
#'   `biophysical`.
#'
#' @return A `poll_map` with a `supply` layer (summed over species) and one
#'   `abundance_<season>` layer per season.
#' @references
#' Kennedy, C.M. et al. (2013) *Ecology Letters*, 16, 584-599.
#' Olsson, O. et al. (2015) *Ecological Modelling*, 316, 133-143.
#' @seealso [poll_biophysical()], [poll_guilds()], [poll_lonsdorf()]
#' @examples
#' library(terra)
#' lulc <- rast(nrows = 40, ncols = 40, xmin = 0, xmax = 8000, ymin = 0, ymax = 8000)
#' values(lulc) <- 211L; lulc[, 1:3] <- 311L
#' bio <- poll_biophysical(data.frame(lucode = c(211L, 311L),
#'   nesting_cavity_availability_index = c(0.10, 0.65),
#'   nesting_ground_availability_index = c(0.30, 0.35),
#'   floral_resources_spring_index     = c(0.40, 0.55),
#'   floral_resources_summer_index     = c(0.35, 0.30)))
#' guilds <- poll_guilds(data.frame(SPECIES = c("Apis", "Osmia"),
#'   nesting_suitability_cavity_index = c(1, 1),
#'   nesting_suitability_ground_index = c(0, 0.1),
#'   foraging_activity_spring_index   = c(1, 1),
#'   foraging_activity_summer_index   = c(1, 0.2),
#'   alpha = c(800, 250), relative_abundance = c(0.9, 0.5)))
#' m <- poll_kennedy(lulc, bio, guilds)
#' plot(m)
#' @export
poll_kennedy <- function(lulc, biophysical, guilds) {
  if (!inherits(lulc, "SpatRaster")) stop("`lulc` must be a SpatRaster.", call. = FALSE)
  if (!inherits(biophysical, "poll_biophysical")) stop("`biophysical` must come from poll_biophysical().", call. = FALSE)
  if (!inherits(guilds, "poll_guilds")) stop("`guilds` must come from poll_guilds().", call. = FALSE)

  substrates <- attr(biophysical, "substrates")
  seasons    <- attr(biophysical, "seasons")
  if (!all(substrates %in% attr(guilds, "substrates")))
    stop("Guild substrates do not match the biophysical table.", call. = FALSE)
  if (!all(seasons %in% attr(guilds, "seasons")))
    stop("Guild seasons do not match the biophysical table.", call. = FALSE)

  Nsub <- lapply(substrates, function(n)
    .reclass(lulc, biophysical$lucode, biophysical[[paste0("nesting_", n, "_availability_index")]]))
  Fsea <- lapply(seasons, function(j)
    .reclass(lulc, biophysical$lucode, biophysical[[paste0("floral_resources_", j, "_index")]]))
  names(Nsub) <- substrates; names(Fsea) <- seasons

  supply_total <- NULL
  ab_season <- setNames(vector("list", length(seasons)), seasons)

  for (si in seq_len(nrow(guilds))) {
    a  <- guilds$alpha[si]
    sa <- guilds$relative_abundance[si]
    ns <- setNames(as.numeric(guilds[si, paste0("nesting_suitability_", substrates, "_index")]), substrates)
    fa <- setNames(as.numeric(guilds[si, paste0("foraging_activity_",   seasons,   "_index")]), seasons)

    HN <- terra::app(terra::rast(lapply(substrates, function(n) Nsub[[n]] * ns[[n]])), fun = max)
    floral_s <- terra::app(terra::rast(lapply(seasons, function(j) Fsea[[j]] * fa[[j]])), fun = sum)

    FR <- poll_kernel(floral_s, alpha = a, normalize = TRUE)
    PS <- FR * HN * sa
    supply_total <- if (is.null(supply_total)) PS else supply_total + PS

    fwd <- poll_kernel(PS, alpha = a, normalize = TRUE)
    for (j in seasons) {
      qual <- (Fsea[[j]] * fa[[j]]) / FR
      qual <- terra::ifel(is.nan(qual) | is.infinite(qual), 0, qual)
      pa   <- fwd * qual
      ab_season[[j]] <- if (is.null(ab_season[[j]])) pa else ab_season[[j]] + pa
    }
  }

  names(supply_total) <- "supply"
  layers <- supply_total
  for (j in seasons) {
    l <- ab_season[[j]]; names(l) <- paste0("abundance_", j)
    layers <- c(layers, l)
  }
  new_poll_map(layers, meta = list(model = "kennedy", species = guilds$SPECIES,
                                    substrates = substrates, seasons = seasons,
                                    alpha = guilds$alpha, call = match.call()))
}

# ---- PollMap modified (quality-weighted, flight-range) model ----------------

#' Quality-weighted, flight-range pollination model (PollMap)
#'
#' The modified Lonsdorf model of Rahimi et al. (2021): as [poll_lonsdorf()], but
#' (i) the foraging neighbourhood is truncated at the flight range `max_dist`
#' (default one foraging distance) and (ii) visitation to each destination is
#' weighted by its floral availability relative to its accessible neighbourhood,
#' \deqn{PSF_i=\mathrm{kernel}(G)\cdot \frac{F_i}{\mathrm{kernel}(F)},}
#' so that habitat edges adjacent to floral-rich fields receive the highest scores.
#'
#' @inheritParams poll_lonsdorf
#' @param max_dist Flight-range truncation (map units); defaults to `guild$alpha`.
#'
#' @return A `poll_map` with `supply` and/or `abundance` layers.
#' @references Rahimi, E. et al. (2021) PollMap: a software for crop pollination
#'   mapping. *Journal of Ecology and Environment*, 45, 22.
#'   Olsson, O. et al. (2015) *Ecological Modelling*, 316, 133-143.
#' @seealso [poll_lonsdorf()], [poll_kennedy()]
#' @examples
#' library(terra)
#' lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
#' values(lulc) <- 1L; lulc[, 1:3] <- 2L
#' lc <- poll_landcover(data.frame(lucode = c(1L, 2L),
#'                                 nesting = c(0.1, 0.9), floral = c(0.6, 0.2)))
#' m <- poll_modified(lulc, lc, poll_guild(400))
#' plot(m)
#' @export
poll_modified <- function(lulc, landcover, guild,
                          max_dist = guild$alpha,
                          layers = c("supply", "abundance")) {
  if (!inherits(lulc, "SpatRaster")) stop("`lulc` must be a SpatRaster.", call. = FALSE)
  if (!inherits(landcover, "poll_landcover")) stop("`landcover` must come from poll_landcover().", call. = FALSE)
  if (!inherits(guild, "poll_guild")) stop("`guild` must come from poll_guild().", call. = FALSE)
  layers <- match.arg(layers, c("supply", "abundance"), several.ok = TRUE)
  alpha  <- guild$alpha

  nesting <- .reclass(lulc, landcover$lucode, landcover$nesting)
  floral  <- .reclass(lulc, landcover$lucode, landcover$floral)

  accessible <- poll_kernel(floral, alpha = alpha, max_dist = max_dist, normalize = TRUE)
  supply     <- nesting * accessible
  names(supply) <- "supply"

  fwd  <- poll_kernel(supply, alpha = alpha, max_dist = max_dist, normalize = TRUE)
  qual <- floral / accessible
  qual <- terra::ifel(is.nan(qual) | is.infinite(qual), 0, qual)
  abundance <- fwd * qual
  names(abundance) <- "abundance"

  out <- terra::subset(c(supply, abundance), layers)
  new_poll_map(out, meta = list(model = "modified", alpha = alpha, guild = guild$name,
                                max_dist = max_dist, call = match.call()))
}

# ---- ESTIMAP relative-potential model ---------------------------------------

#' ESTIMAP relative pollination potential
#'
#' Implements the resource-plus-activity index of Zulian et al. (2013): a relative
#' availability score per land-cover class is smoothed over the foraging
#' neighbourhood and scaled by an activity multiplier combining abiotic and
#' landscape modifiers,
#' \deqn{PP(x)=C(x)\cdot \mathrm{kernel}(RA)(x).}
#'
#' @param lulc A land-cover `SpatRaster`.
#' @param scores A data frame with columns `lucode` and `availability` (in
#'   \[0, 1\]) giving the relative pollinator availability of each class.
#' @param guild A [poll_guild()] giving `alpha`.
#' @param activity Optional activity multiplier: a single `SpatRaster` in \[0, 1\],
#'   or a list of such layers combined by multiplication (e.g. standardised solar
#'   radiation, temperature).
#' @param edge Optional additive enhancement layer (e.g. proximity to forest edge
#'   or riparian zones), added to `availability` before smoothing.
#'
#' @return A `poll_map` of relative pollination potential.
#' @references Zulian, G., Maes, J. & Paracchini, M.L. (2013) *Land*, 2, 472-492.
#'   \doi{10.3390/land2030472}
#' @seealso [poll_mce()]
#' @examples
#' library(terra)
#' lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
#' values(lulc) <- 1L; lulc[, 1:3] <- 2L
#' scores <- data.frame(lucode = c(1L, 2L), availability = c(0.4, 0.7))
#' m <- poll_estimap(lulc, scores, poll_guild(500))
#' plot(m)
#' @export
poll_estimap <- function(lulc, scores, guild, activity = NULL, edge = NULL) {
  if (!inherits(lulc, "SpatRaster")) stop("`lulc` must be a SpatRaster.", call. = FALSE)
  scores <- as.data.frame(scores)
  if (!all(c("lucode", "availability") %in% names(scores)))
    stop("`scores` must have columns `lucode` and `availability`.", call. = FALSE)
  if (any(scores$availability < 0 | scores$availability > 1, na.rm = TRUE))
    stop("`availability` must lie in [0, 1].", call. = FALSE)
  if (!inherits(guild, "poll_guild")) stop("`guild` must come from poll_guild().", call. = FALSE)

  RA <- .reclass(lulc, as.integer(scores$lucode), scores$availability)
  if (!is.null(edge)) RA <- RA + edge
  neigh <- poll_kernel(RA, alpha = guild$alpha, normalize = TRUE)

  C <- NULL
  if (!is.null(activity)) {
    if (inherits(activity, "SpatRaster")) {
      C <- if (terra::nlyr(activity) > 1) terra::app(activity, prod) else activity
    } else if (is.list(activity)) {
      C <- Reduce(`*`, activity)
    } else stop("`activity` must be a SpatRaster or a list of SpatRasters.", call. = FALSE)
  }
  pp <- if (is.null(C)) neigh else neigh * C
  names(pp) <- "potential"
  new_poll_map(pp, meta = list(model = "estimap", alpha = guild$alpha, call = match.call()))
}

# ---- Multi-criteria evaluation ----------------------------------------------

# Analytic-hierarchy-process weights from a pairwise-comparison matrix.
.ahp_weights <- function(mat) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  ev <- eigen(mat)
  k  <- which.max(Re(ev$values))
  w  <- Re(ev$vectors[, k]); w <- w / sum(w)
  lambda <- Re(ev$values[k])
  CI <- (lambda - n) / (n - 1)
  RI <- c(0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49)
  ri <- if (n <= length(RI)) RI[n] else 1.49
  CR <- if (is.na(ri) || ri == 0) 0 else CI / ri
  list(weights = stats::setNames(w, rownames(mat)), CR = CR, lambda = lambda)
}

#' Multi-criteria evaluation by weighted linear combination
#'
#' Integrates several criterion layers into one suitability surface. Each
#' criterion is (optionally) standardised to \[0, 1\] and combined as a weighted
#' sum, \eqn{S(x)=\sum_k w_k x_k(x)}, with optional Boolean constraint masks. This
#' is the overlay device used by ESTIMAP and PollMap (Zulian et al. 2013; Rahimi
#' et al. 2021), following weighted-linear-combination practice (Malczewski 2000)
#' with weights optionally from the analytic hierarchy process (Saaty 1977).
#'
#' @param criteria A named list of single-layer `SpatRaster`s (or a multi-layer
#'   `SpatRaster`).
#' @param weights Either a named numeric vector (normalised to sum to 1) or a
#'   square pairwise-comparison matrix, in which case AHP priorities and the
#'   consistency ratio are computed.
#' @param directions Optional named vector of `"benefit"` (higher is better,
#'   default) or `"cost"` (higher is worse) per criterion.
#' @param standardize Logical; min-max standardise each criterion to \[0, 1\]
#'   (default `TRUE`).
#' @param constraints Optional list of Boolean (`0`/`1`) `SpatRaster`s multiplied
#'   into the result.
#'
#' @return A single-layer `SpatRaster` of the integrated suitability, with the
#'   weights and (if applicable) consistency ratio attached as attributes.
#' @references Saaty, T.L. (1977) *Journal of Mathematical Psychology*, 15,
#'   234-281. Malczewski, J. (2000) *Transactions in GIS*, 4, 5-22.
#' @seealso [poll_estimap()]
#' @examples
#' library(terra)
#' r <- rast(nrows = 20, ncols = 20)
#' access   <- setValues(r, runif(400))
#' resource <- setValues(r, runif(400))
#' s <- poll_mce(list(access = access, resource = resource),
#'               weights = c(access = 0.4, resource = 0.6))
#' plot(s)
#' @export
poll_mce <- function(criteria, weights, directions = NULL,
                     standardize = TRUE, constraints = NULL) {
  if (inherits(criteria, "SpatRaster")) {
    criteria <- stats::setNames(lapply(seq_len(terra::nlyr(criteria)),
                                       function(i) criteria[[i]]), names(criteria))
  }
  if (!is.list(criteria) || is.null(names(criteria)))
    stop("`criteria` must be a named list of SpatRasters or a named SpatRaster.", call. = FALSE)
  nms <- names(criteria)

  CR <- NA_real_
  if (is.matrix(weights) || inherits(weights, "matrix")) {
    aw <- .ahp_weights(weights); w <- aw$weights; CR <- aw$CR
    if (is.null(names(w))) names(w) <- nms
  } else {
    if (is.null(names(weights))) names(weights) <- nms
    w <- weights[nms] / sum(weights[nms])
  }
  if (!setequal(names(w), nms)) stop("Names of `weights` must match `criteria`.", call. = FALSE)

  terms <- lapply(nms, function(nm) {
    r <- criteria[[nm]]
    if (standardize) {
      mm <- terra::minmax(r); lo <- mm[1]; hi <- mm[2]
      s <- if (hi > lo) (r - lo) / (hi - lo) else r * 0
      if (!is.null(directions) && !is.na(directions[nm]) && directions[nm] == "cost") s <- 1 - s
      r <- s
    }
    r * w[[nm]]
  })
  out <- Reduce(`+`, terms)
  if (!is.null(constraints)) out <- out * Reduce(`*`, constraints)
  names(out) <- "suitability"
  attr(out, "weights") <- w
  attr(out, "CR") <- CR
  out
}
