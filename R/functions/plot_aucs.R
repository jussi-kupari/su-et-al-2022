plot_aucs <- function(x, lower_limit = 0.5, upper_limit = 1) {
  x %>%
    ggplot(aes(time, auc)) +
    geom_segment(aes(xend = time, yend = 0.5)) +
    geom_point() +
    geom_hline(yintercept = 0.7, color = "red") +
    cowplot::theme_cowplot() +
    coord_flip() +
    ylim(lower_limit, upper_limit) +
    labs(y = "", x = "", title = unique(x$cell_type))
}
