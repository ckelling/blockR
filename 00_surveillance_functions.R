#adapted from Surveillance R package
#original function described here: https://surveillance.r-forge.r-project.org/pkgdown/reference/polyAtBorder.html
#taking out unionSpatialPolygons step

polyAtBorder_cust <- function (SpP, W=W,
                          snap = sqrt(.Machine$double.eps),
                          method = "sf", ...)
{
  SpP <- as(SpP, "SpatialPolygons")
  if (length(W@polygons) > 1)
    warning("unionSpatialPolygons() produced >1 Polygons-components")
  Wcoords <- unique(do.call("rbind",
                            lapply(W@polygons[[1]]@Polygons, coordinates)))
  atBorder <- sapply(SpP@polygons, function (x) {
    coords <- unique(do.call("rbind", lapply(x@Polygons, coordinates)))
    res <- FALSE
    for (i in seq_len(nrow(coords))) {
      if (any(spDistsN1(Wcoords, coords[i,], longlat=FALSE) < snap)) {
        res <- TRUE
        break
      }
    }
    res
  })
  names(atBorder) <- row.names(SpP)
  atBorder
}
