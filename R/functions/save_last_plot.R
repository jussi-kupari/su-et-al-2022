save_last_plot <- function(filename, type = "svg", w_mm = 200, h_mm = 200) {
  ggsave(
    filename = filename,
    plot = last_plot(),
    device = type,
    # path = here("ouputs", "fig_drafts"),
    scale = 1,
    width = w_mm,
    height = h_mm,
    units = "mm",
    dpi = 300
  )
}
