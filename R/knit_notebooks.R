
# Notebook directory
notebook_dir <- paste0(here::here(), "/notebooks/")

# Create knitting function
knit_notebook <- function(file) {
  ezknitr::ezknit(
    file = paste0(notebook_dir, file),
    wd = notebook_dir,
    out_dir = notebook_dir,
    verbose = TRUE,
    keep_md = TRUE,
    keep_html = FALSE
  )
}

# Style all Notebooks
styler::style_dir(notebook_dir)

# Knit all Notebooks
fs::dir_ls("notebooks") |>
  stringr::str_subset(".Rmd$") |>
  stringr::str_remove("notebooks/") |>
  purrr::walk(knit_notebook)


