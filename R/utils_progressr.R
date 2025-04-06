#' Set the best available progressr handler
#' @keywords internal
set_fastfocal_progress_handler <- function() {
  if (!requireNamespace("progressr", quietly = TRUE)) return()
  
  available_handlers <- progressr::handlers()
  
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      "handler_rstudio" %in% names(available_handlers)) {
    progressr::handlers("rstudio")
  } else {
    progressr::handlers("txtprogressbar")
  }
}
