#' Write gzipped JSON
#'
#' Serializes an object to JSON and writes to a gzipped file.
#'
#' @param x An object to be serialized to JSON.
#' @param path File on disk (should end in '.gz').
#' @param ... Additional arguments passed to [jsonlite::write_json()]
#'
#' @export
write_jsongz <- function(x, path, ...) {
  assertthat::assert_that(endsWith(path, ".gz"))
  fs::dir_create(dirname(path))
  gz <- gzfile(path, open = "w")
  jsonlite::write_json(x = x, path = gz, ...)
  close(gz)
}

#' Print current timestamp for logging
#'
#' @return Current timestamp as character.
#' @export
date_log <- function() {
  as.character(paste0("[", as.POSIXct(Sys.time()), "]"))
}

#' Does R Package Exist
#'
#' Checks if the specified R package exists on the local system.
#'
#' @param p The R package to check for.
#' @return TRUE if package exists, FALSE otherwise.
#'
#' @export
pkg_exists <- function(p) {
  assertthat::assert_that(is.character(p))
  nzchar(system.file(package = p))
}

#' @noRd
dummy1 <- function() {
  # Solves R CMD check: Namespaces in Imports field not imported from
  argparse::ArgumentParser
}
