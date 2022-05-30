
set_active_assay <- function(object, assay) {
  DefaultAssay(object) <- assay

  return(object)
}
