predict_labels <- function(query, reference, dims, k.weight = 50) {
  predictions <- TransferData(
    anchorset = FindTransferAnchors(
      reference = reference,
      query = query,
      dims = dims
    ),
    k.weight = k.weight,
    refdata = reference@active.ident,
    dims = dims
  )

  query %<>%
    AddMetaData(metadata = predictions)

  return(query)
}
