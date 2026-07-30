# S3 methods -------------------------------------------------------------------

#' @export
print.poll_guild <- function(x, ...) {
  cat("<poll_guild>\n")
  cat("  name:  ", x$name, "\n", sep = "")
  cat("  alpha: ", x$alpha, " (map units)\n", sep = "")
  invisible(x)
}

#' @export
print.poll_landcover <- function(x, ...) {
  cat("<poll_landcover>: ", nrow(x), " land-cover classes\n", sep = "")
  print(utils::head(as.data.frame(x), 10L), row.names = FALSE)
  if (nrow(x) > 10L) cat("  ... ", nrow(x) - 10L, " more\n", sep = "")
  invisible(x)
}

#' Print a poll_map
#' @param x A \code{poll_map} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.poll_map <- function(x, ...) {
  r <- x$map
  res <- terra::res(r)
  cat("<poll_map>\n")
  cat("  model:  ", x$meta$model %||% "NA", "\n", sep = "")
  cat("  layers: ", paste(names(r), collapse = ", "), "\n", sep = "")
  cat("  alpha:  ", x$meta$alpha %||% NA, " (map units)\n", sep = "")
  cat(sprintf("  grid:   %d rows x %d cols, resolution %s\n",
              terra::nrow(r), terra::ncol(r),
              paste(signif(res, 4), collapse = " x ")))
  invisible(x)
}

#' Summarise a poll_map
#'
#' @param object A \code{poll_map} object.
#' @param ... Ignored.
#' @return A data frame of per-layer minimum, mean and maximum (invisibly).
#' @export
summary.poll_map <- function(object, ...) {
  r <- object$map
  stats_df <- terra::global(r, fun = c("min", "mean", "max"), na.rm = TRUE)
  stats_df <- data.frame(layer = rownames(stats_df), stats_df, row.names = NULL)
  print(stats_df)
  invisible(stats_df)
}

#' Plot a poll_map
#'
#' @param x A \code{poll_map} object.
#' @param ... Passed to \code{\link[terra]{plot}}.
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#' @export
plot.poll_map <- function(x, ...) {
  terra::plot(x$map, ...)
  invisible(NULL)
}
