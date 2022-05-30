clear_libraries <- function() {
  
  pacs <- names(sessionInfo()$otherPkgs)
  
  if (!(is.null(pacs))) {
    suppressWarnings(invisible(lapply(
      paste("package:", pacs, sep = ""),
      detach,
      character.only = TRUE,
      unload = TRUE
    )))
  }
}
