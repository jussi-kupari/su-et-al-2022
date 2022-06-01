clear_environment <- function(except = NULL, pattern = ls()) {
    except = except
    pattern = pattern
    formula = c(c(except), ls(pattern = pattern, envir = .GlobalEnv))
    rm(list = setdiff(ls(envir = .GlobalEnv), formula), envir = .GlobalEnv)
  }
