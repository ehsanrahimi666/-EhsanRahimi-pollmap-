# Occurrence handling: thinning and model evaluation #########################

#' Spatially thin occurrence records
#'
#' Reduces spatial clustering in occurrence data by retaining at most one record
#' per grid cell of a chosen size, optionally followed by a random cap.
#' Opportunistic records accumulate around towns, roads and popular recording
#' sites; fitting a distribution model to unthinned data allows that sampling
#' effort to be modelled as habitat preference, and the resulting spatial
#' autocorrelation inflates apparent predictive performance. Thinning is a
#' standard remedy (Boria et al. 2014).
#'
#' Thinning is performed on a projected grid, so `crs` should be a projected
#' (metric) coordinate reference system appropriate to the study region.
#'
#' @param occ A data frame of occurrences.
#' @param coords Length-2 character vector naming the longitude and latitude
#'   columns.
#' @param dist Minimum spacing in kilometres; one record is kept per `dist` by
#'   `dist` grid cell.
#' @param crs_from Coordinate reference system of `occ` (default `"EPSG:4326"`).
#' @param crs_to Projected CRS in which thinning is carried out.
#' @param max_n Optional cap on the number of retained records; a random subset
#'   is taken if more survive thinning.
#' @param seed Optional random seed, for reproducibility of the cap.
#'
#' @return A data frame with the same columns as `occ`, containing the retained
#'   rows, with an attribute `thin` recording the input and output counts.
#' @references Boria, R.A., Olson, L.E., Goodman, S.M. & Anderson, R.P. (2014)
#'   Spatial filtering to reduce sampling bias can improve the performance of
#'   ecological niche models. *Ecological Modelling*, 275, 73-77.
#' @seealso [poll_sdm()], [poll_sdm_eval()]
#' @examples
#' set.seed(1)
#' # 200 records clustered in one place plus 40 spread out
#' occ <- data.frame(
#'   Longitude = c(rnorm(200, 8.5, 0.02), runif(40, 6, 10)),
#'   Latitude  = c(rnorm(200, 47.4, 0.02), runif(40, 46, 47.7)))
#' thinned <- poll_thin(occ, dist = 5)
#' nrow(occ)
#' nrow(thinned)
#' @export
poll_thin <- function(occ, coords = c("Longitude", "Latitude"), dist = 5,
                      crs_from = "EPSG:4326", crs_to = "EPSG:2056",
                      max_n = NULL, seed = NULL) {
  occ <- as.data.frame(occ)
  if (!all(coords %in% names(occ)))
    stop("`occ` must contain the coordinate columns ",
         paste(coords, collapse = " and "), ".", call. = FALSE)
  if (!is.numeric(dist) || length(dist) != 1L || dist <= 0)
    stop("`dist` must be a single positive number (kilometres).", call. = FALSE)
  n_in <- nrow(occ)
  if (n_in == 0L) return(occ)
  if (!is.null(seed)) set.seed(seed)

  keep <- stats::complete.cases(occ[, coords])
  occ <- occ[keep, , drop = FALSE]
  if (nrow(occ) == 0L) return(occ)

  v <- terra::vect(occ, geom = coords, crs = crs_from)
  v <- terra::project(v, crs_to)
  xy <- terra::crds(v)
  cell <- paste(floor(xy[, 1] / (dist * 1000)),
                floor(xy[, 2] / (dist * 1000)), sep = "_")
  occ <- occ[!duplicated(cell), , drop = FALSE]

  if (!is.null(max_n) && nrow(occ) > max_n)
    occ <- occ[sort(sample.int(nrow(occ), max_n)), , drop = FALSE]

  attr(occ, "thin") <- c(input = n_in, retained = nrow(occ), dist_km = dist)
  occ
}

# Mann-Whitney (rank-based) AUC; no external dependency.
.auc <- function(obs, pred) {
  pos <- pred[obs == 1]; neg <- pred[obs == 0]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

# Maximum of Youden's J over candidate thresholds.
.tss <- function(obs, pred) {
  if (!length(unique(obs)) == 2L) return(NA_real_)
  thr <- sort(unique(stats::quantile(pred, probs = seq(0, 1, 0.01), na.rm = TRUE)))
  best <- -Inf
  for (t in thr) {
    p <- as.integer(pred >= t)
    tp <- sum(p == 1 & obs == 1); fn <- sum(p == 0 & obs == 1)
    tn <- sum(p == 0 & obs == 0); fp <- sum(p == 1 & obs == 0)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    if (!is.na(sens) && !is.na(spec)) best <- max(best, sens + spec - 1)
  }
  if (is.finite(best)) best else NA_real_
}

#' Cross-validated evaluation of a species-distribution model
#'
#' Assesses the predictive performance of the [poll_sdm()] pathway by k-fold
#' cross-validation, reporting the area under the receiver-operating
#' characteristic curve (AUC) and the true skill statistic (TSS). Folds may be
#' assigned at random or as spatial blocks; spatial blocking gives a less
#' optimistic and generally more honest estimate when records are spatially
#' autocorrelated (Roberts et al. 2017).
#'
#' @param occ A data frame of occurrences.
#' @param predictors A multi-layer `SpatRaster` of environmental predictors.
#' @param method Algorithm(s) passed to [poll_sdm()].
#' @param k Number of folds.
#' @param block Logical; if `TRUE` (default) folds are contiguous spatial blocks
#'   rather than random subsets.
#' @param coords Length-2 names of the coordinate columns.
#' @param n_background Number of background points sampled per fold.
#' @param seed Optional random seed.
#'
#' @return A list of class `poll_eval` with a per-fold data frame and the mean
#'   AUC and TSS.
#' @references
#' Roberts, D.R. et al. (2017) Cross-validation strategies for data with
#' temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40,
#' 913-929.
#' Allouche, O., Tsoar, A. & Kadmon, R. (2006) Assessing the accuracy of species
#' distribution models: prevalence, kappa and the true skill statistic (TSS).
#' *Journal of Applied Ecology*, 43, 1223-1232.
#' @seealso [poll_sdm()], [poll_thin()]
#' @examples
#' library(terra)
#' set.seed(1)
#' preds <- rast(nrows = 30, ncols = 30, nlyr = 2,
#'               xmin = 0, xmax = 300, ymin = 0, ymax = 300)
#' names(preds) <- c("v1", "v2")
#' values(preds) <- cbind(runif(900), runif(900))
#' xy  <- xyFromCell(preds, which(values(preds[["v1"]]) > 0.6))
#' occ <- data.frame(Longitude = xy[, 1], Latitude = xy[, 2])
#' ev <- poll_sdm_eval(occ, preds, k = 3, block = FALSE,
#'                     n_background = 200, seed = 1)
#' ev$mean_auc
#' @export
poll_sdm_eval <- function(occ, predictors, method = "glm", k = 5, block = TRUE,
                          coords = c("Longitude", "Latitude"),
                          n_background = 1000, seed = NULL) {
  if (!inherits(predictors, "SpatRaster"))
    stop("`predictors` must be a SpatRaster.", call. = FALSE)
  occ <- as.data.frame(occ)
  if (!all(coords %in% names(occ)))
    stop("`occ` must contain the coordinate columns ",
         paste(coords, collapse = " and "), ".", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(occ)
  if (n < 2 * k) stop("Too few occurrences for ", k, "-fold cross-validation.",
                      call. = FALSE)

  if (block) {
    # contiguous blocks along the dominant axis of the point cloud; the SAME
    # coordinate breaks are applied to the background pool, so that presences
    # and background in a test fold come from the same region. Evaluating
    # block-restricted presences against study-wide background would compare
    # different areas and bias the estimate downwards.
    rx <- diff(range(occ[[coords[1]]])); ry <- diff(range(occ[[coords[2]]]))
    ax <- if (rx >= ry) coords[1] else coords[2]
    brks <- stats::quantile(occ[[ax]], probs = seq(0, 1, length.out = k + 1),
                            na.rm = TRUE)
    brks <- unique(brks)
    if (length(brks) < 2) stop("Cannot form spatial blocks from these coordinates.",
                               call. = FALSE)
    brks[1] <- -Inf; brks[length(brks)] <- Inf
    fold <- cut(occ[[ax]], breaks = brks, labels = FALSE, include.lowest = TRUE)
    k <- max(fold, na.rm = TRUE)
  } else {
    ax <- NULL; brks <- NULL
    fold <- sample(rep_len(seq_len(k), n))
  }

  # one background pool, assigned to folds by the same rule
  bg_pts <- terra::spatSample(predictors[[1]], min(5 * n_background, 10000),
                              method = "random", as.points = TRUE,
                              na.rm = TRUE, values = FALSE)
  bg_xy <- terra::crds(bg_pts)
  if (block) {
    bg_ax <- if (identical(ax, coords[1])) bg_xy[, 1] else bg_xy[, 2]
    bg_fold <- cut(bg_ax, breaks = brks, labels = FALSE, include.lowest = TRUE)
  } else {
    bg_fold <- sample(rep_len(seq_len(k), nrow(bg_xy)))
  }

  res <- data.frame(fold = seq_len(k), n_test = NA_integer_,
                    auc = NA_real_, tss = NA_real_)
  for (i in seq_len(k)) {
    train <- occ[fold != i, , drop = FALSE]
    test  <- occ[fold == i, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 1) next
    fit <- tryCatch(
      poll_sdm(train, predictors, method = method, coords = coords,
               n_background = n_background),
      error = function(e) NULL)
    if (is.null(fit)) next
    sui <- poll_sdm_project(fit, predictors)$map

    pres_pts <- terra::vect(test, geom = coords, crs = terra::crs(predictors))
    bg_i <- bg_pts[which(bg_fold == i)]
    if (length(bg_i) < 10) next
    p_pres <- terra::extract(sui, pres_pts, ID = FALSE)[, 1]
    p_bg   <- terra::extract(sui, bg_i,     ID = FALSE)[, 1]
    obs  <- c(rep(1L, length(p_pres)), rep(0L, length(p_bg)))
    pred <- c(p_pres, p_bg)
    ok <- is.finite(pred)
    obs <- obs[ok]; pred <- pred[ok]

    res$n_test[i] <- nrow(test)
    res$auc[i] <- .auc(obs, pred)
    res$tss[i] <- .tss(obs, pred)
  }

  out <- list(folds = res,
              mean_auc = mean(res$auc, na.rm = TRUE),
              mean_tss = mean(res$tss, na.rm = TRUE),
              k = k, block = block, method = method)
  structure(out, class = "poll_eval")
}

#' @export
print.poll_eval <- function(x, ...) {
  cat("<poll_eval>: ", x$k, "-fold ", if (x$block) "spatial-block" else "random",
      " CV | method: ", paste(x$method, collapse = ", "), "\n", sep = "")
  cat("  mean AUC = ", round(x$mean_auc, 3),
      " | mean TSS = ", round(x$mean_tss, 3), "\n", sep = "")
  invisible(x)
}
