clean_slate <- function() {
  # Clear libraries
  invisible(lapply(
    paste("package:", names(sessionInfo()$otherPkgs), sep = ""),
    detach,
    character.only = TRUE,
    unload = TRUE
  ))

  # Clear objects
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

  # Restart
  .rs.restartR()
}
