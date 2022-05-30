
get_upset_gene_lists <- function(genes_by_timepoint) {
  
  upset_gene_lists <-
    overlapGroups(genes_by_timepoint)
  
  gene_lists <-
    upset_gene_lists %>% 
    map(~ enframe(.x) %>% select(gene = name)) %>% 
    set_names(str_replace_all(names(upset_gene_lists), ":", "_"))
  
  gene_lists
}