extract_top_GOs <- function(GO_list) {
  GO_list %>%
    map( ~ {
      .x[1,] <- names(.x)
      return(.x)
    })
}
