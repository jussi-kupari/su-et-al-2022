get_enriched <- function(genes, dbs_vector) {
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  
  if (is.null(dbs)) {
    websiteLive <- FALSE
  }
  
  dbs <- dbs_vector
  
  if (websiteLive) {
    enriched <- enrichr(genes, dbs)
  }
  
  return(enriched)
}
