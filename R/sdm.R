# Species-distribution-model pathway ##########################################

new_poll_model <- function(fits, method, predictors, meta = list()) {
  structure(list(fits = fits, method = method, predictors = predictors, meta = meta),
            class = "poll_model")
}

#' @export
print.poll_model <- function(x, ...) {
  cat("<poll_model>: method(s) ", paste(x$method, collapse = ", "),
      " | predictors: ", paste(x$predictors, collapse = ", "), "\n", sep = "")
  if (!is.null(x$meta$n_pres))
    cat("  presences: ", x$meta$n_pres, " | background: ", x$meta$n_bg, "\n", sep = "")
  invisible(x)
}

#' Fit a species-distribution model for a pollinator
#'
#' Fits one or more correlative models relating occurrence to environmental
#' predictors and returns them as an ensemble (Elith et al. 2011; Perennes et al.
#' 2021). Presence records are complemented with randomly sampled background
#' points. The self-contained algorithms are `"glm"` (always available), `"gam"`
#' (needs \pkg{mgcv}) and `"rf"` (needs \pkg{randomForest}); requesting several
#' fits an ensemble. For a fuller workflow (spatial-block cross-validation,
#' additional algorithms, ensemble weighting) the \pkg{flexsdm} package can be
#' used to build the fits externally and passed through unchanged.
#'
#' @param occ A data frame of occurrences with coordinate columns.
#' @param predictors A multi-layer `SpatRaster` of environmental predictors.
#' @param method Character vector of algorithms to fit (`"glm"`, `"gam"`, `"rf"`).
#' @param coords Length-2 names of the coordinate columns.
#' @param n_background Number of background points to sample.
#' @param seed Optional random seed.
#'
#' @return An object of class `poll_model`.
#' @references Perennes, M. et al. (2021) *Ecological Modelling*, 444, 109484.
#' @seealso [poll_sdm_project()], [poll_sdm_supply()]
#' @examples
#' library(terra)
#' preds <- rast(nrows = 30, ncols = 30, nlyr = 2,
#'               xmin = 0, xmax = 300, ymin = 0, ymax = 300)
#' names(preds) <- c("v1", "v2"); values(preds) <- cbind(runif(900), runif(900))
#' xy  <- xyFromCell(preds, which(values(preds[["v1"]]) > 0.6))
#' occ <- data.frame(Longitude = xy[, 1], Latitude = xy[, 2])
#' fit <- poll_sdm(occ, preds, method = "glm", n_background = 200)
#' fit
#' @export
poll_sdm <- function(occ, predictors, method = "glm",
                     coords = c("Longitude", "Latitude"),
                     n_background = 1000, seed = NULL) {
  if (!inherits(predictors, "SpatRaster")) stop("`predictors` must be a SpatRaster.", call. = FALSE)
  method <- match.arg(method, c("glm", "gam", "rf"), several.ok = TRUE)
  if (!is.null(seed)) set.seed(seed)
  occ <- as.data.frame(occ)
  if (!all(coords %in% names(occ))) stop("`occ` must contain columns ", paste(coords, collapse = ", "), ".", call. = FALSE)

  pres <- terra::vect(occ, geom = coords, crs = terra::crs(predictors))
  bg   <- terra::spatSample(predictors[[1]], n_background, method = "random",
                            as.points = TRUE, na.rm = TRUE, values = FALSE)
  pe <- terra::extract(predictors, pres, ID = FALSE)
  be <- terra::extract(predictors, bg,   ID = FALSE)
  df <- rbind(data.frame(presence = 1, pe), data.frame(presence = 0, be))
  df <- df[stats::complete.cases(df), ]
  vars <- names(predictors)
  form <- stats::as.formula(paste("presence ~", paste(vars, collapse = " + ")))

  fits <- list()
  for (m in method) {
    fits[[m]] <- switch(m,
      glm = stats::glm(form, data = df, family = stats::binomial()),
      gam = {
        if (!requireNamespace("mgcv", quietly = TRUE)) stop("method 'gam' needs the mgcv package.", call. = FALSE)
        sform <- stats::as.formula(paste("presence ~", paste(sprintf("s(%s)", vars), collapse = " + ")))
        mgcv::gam(sform, data = df, family = stats::binomial())
      },
      rf = {
        if (!requireNamespace("randomForest", quietly = TRUE)) stop("method 'rf' needs the randomForest package.", call. = FALSE)
        df$presence <- factor(df$presence)
        randomForest::randomForest(form, data = df)
      })
  }
  new_poll_model(fits, method, vars,
                 meta = list(n_pres = sum(df$presence == 1), n_bg = sum(df$presence == 0)))
}

#' Project a fitted SDM to a suitability surface
#'
#' Predicts habitat suitability across a predictor raster from a [poll_sdm()]
#' model, averaging over the ensemble members (Polce et al. 2013).
#'
#' @param model A `poll_model` from [poll_sdm()].
#' @param predictors A `SpatRaster` of predictors (current or scenario) with the
#'   layers used to fit the model.
#'
#' @return A `poll_map` with a single `suitability` layer in \[0, 1\].
#' @references Polce, C. et al. (2013) *PLoS ONE*, 8, e76308.
#' @seealso [poll_sdm()], [poll_sdm_supply()]
#' @examples
#' library(terra)
#' preds <- rast(nrows = 30, ncols = 30, nlyr = 2,
#'               xmin = 0, xmax = 300, ymin = 0, ymax = 300)
#' names(preds) <- c("v1", "v2"); values(preds) <- cbind(runif(900), runif(900))
#' xy  <- xyFromCell(preds, which(values(preds[["v1"]]) > 0.6))
#' occ <- data.frame(Longitude = xy[, 1], Latitude = xy[, 2])
#' fit <- poll_sdm(occ, preds, method = "glm", n_background = 200)
#' sui <- poll_sdm_project(fit, preds)
#' plot(sui)
#' @export
poll_sdm_project <- function(model, predictors) {
  if (!inherits(model, "poll_model")) stop("`model` must be a poll_model.", call. = FALSE)
  preds <- lapply(names(model$fits), function(m) {
    f <- model$fits[[m]]
    if (m == "rf") {
      terra::predict(predictors, f, type = "prob", index = 2, na.rm = TRUE)
    } else {
      terra::predict(predictors, f, type = "response", na.rm = TRUE)
    }
  })
  sui <- if (length(preds) > 1) terra::app(terra::rast(preds), mean) else preds[[1]]
  names(sui) <- "suitability"
  new_poll_map(sui, meta = list(model = "sdm", method = model$method, call = match.call()))
}

#' @rdname poll_sdm_project
#' @param object A `poll_model`.
#' @param ... Unused.
#' @export
predict.poll_model <- function(object, predictors, ...) {
  poll_sdm_project(object, predictors)
}

#' Couple SDM suitabilities into a community supply surface
#'
#' Aggregates per-species suitability surfaces into a community pollinator-supply
#' index, optionally weighted by relative abundance or pollination effectiveness
#' (Polce et al. 2013; Perennes et al. 2021).
#'
#' @param suitabilities A list of `poll_map`s or single-layer `SpatRaster`s (one
#'   per species), or a multi-layer `SpatRaster`.
#' @param weights Optional numeric weights (normalised to sum to 1); equal by
#'   default.
#'
#' @return A `poll_map` with a single `supply` layer.
#' @references Polce, C. et al. (2013) *PLoS ONE*, 8, e76308.
#' @seealso [poll_sdm_project()], [poll_balance()]
#' @examples
#' library(terra)
#' r <- rast(nrows = 20, ncols = 20)
#' s1 <- setValues(r, runif(400)); s2 <- setValues(r, runif(400))
#' supply <- poll_sdm_supply(list(s1, s2), weights = c(0.7, 0.3))
#' plot(supply)
#' @export
poll_sdm_supply <- function(suitabilities, weights = NULL) {
  if (inherits(suitabilities, "SpatRaster")) {
    rl <- lapply(seq_len(terra::nlyr(suitabilities)), function(i) suitabilities[[i]])
  } else {
    rl <- lapply(suitabilities, function(s) if (inherits(s, "poll_map")) s$map[[1]] else s[[1]])
  }
  if (is.null(weights)) weights <- rep(1, length(rl))
  if (length(weights) != length(rl)) stop("`weights` length must match the number of suitabilities.", call. = FALSE)
  weights <- weights / sum(weights)
  supply <- Reduce(`+`, Map(function(r, w) r * w, rl, as.list(weights)))
  names(supply) <- "supply"
  new_poll_map(supply, meta = list(model = "sdm_supply", call = match.call()))
}
