get_regulated_genes <- function(df, by, direction, p_adj = 10^-10, log2FC = 0.25) {
  
  stopifnot("by %not_in% c('celltype', 'timepoint')"  =
              by %in% c("celltype", "timepoint")
  )
  
  stopifnot("direction %not_in% c('up', 'down', 'both')"  =
              direction %in% c("up", "down", "both")
  )
  
  
  if (by == "celltype") {
    second_split <- "timepoint"
  } else {
    second_split <- "celltype"
  }
  
  data <-
    df %>% split(.[[by]]) %>%
    map(~ .x %<>% mutate(abs_log2FC = abs(avg_log2FC))) %>%
    map(~ .x  %<>% filter(p_val_adj < p_adj, abs_log2FC > (abs(log2FC))))
  
  
  get_up <-
    . %>%
    filter(avg_log2FC > 0) %>%
    split(.[[second_split]]) %>%
    map(pull, gene)
  
  get_down <-
    . %>%
    filter(avg_log2FC < 0) %>%
    split(.[[second_split]]) %>%
    map(pull, gene)
  
  get_both <-
    . %>%
    filter(abs_log2FC > 0) %>%
    split(.[[second_split]]) %>%
    map(pull, gene)
  
  
  if (direction == "up") {
    data %>% map(get_up)
    
  } else if (direction == "down") {
    data %>% map(get_down)
    
  } else {
    data %>% map(get_both)
    
  }
}
