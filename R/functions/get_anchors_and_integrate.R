
get_anchors_and_integrate <- function(object, dims) {
  object <- IntegrateData(
    anchorset = FindIntegrationAnchors(
      object.list = SplitObject(object, split.by = "orig.ident"), dims = dims
    ), dims = dims
  )

  return(object)
}
