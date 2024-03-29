% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_detection.R
\name{tree_extraction}
\alias{tree_extraction}
\title{Tree extraction}
\usage{
tree_extraction(
  r_dem_nl,
  r_maxi = NULL,
  r_dem_w = NULL,
  r_mask = NULL,
  crown = FALSE
)
}
\arguments{
\item{r_dem_nl}{SpatRaster object. Output raster of \code{\link{tree_segmentation}}. Otherwise a raster of canopy height model,
preferably filtered to avoid effect of holes on volume and surface computation can be provided. In this case arguments `r_maxi`, `r_dem_w` 
have to be provided}

\item{r_maxi}{SpatRaster object. raster with positive values at local maxima (in case `r_dem_nl` does not contain it)}

\item{r_dem_w}{SpatRaster object. segmented raster  (in case `r_dem_nl` does not contain it)}

\item{r_mask}{SpatRaster object. only segments which maxima are inside the mask are extracted. Values should be NA outside the mask, 1 inside.}

\item{crown}{boolean. Should the 2D crown geometry be added in wkt format 
to the output data.frame ?}
}
\value{
A \code{sf} collection of POINTs with 7 fields: tree id, local maximum stats
 (height, dominance radius), segment stats (surface and volume), coordinates 
 (x and y). In case argument `crown` is `TRUE`, a `crown` field 
 containing the WKT geometry of the 2D crown is also present. Coordinates are 
 written with one decimal to the right of the order of magnitude of 
 the SpatRaster resolution (e.g. if resolution is 1/3 then 2 decimals are written).
}
\description{
creates a data.frame with segment id, height and coordinates of maxima, surface and volume, computed from three images: 
initial, local maxima and segmented, obtained with \code{\link{tree_segmentation}}. The 2D polygon associated to each crown 
can be added as a WKT field
}
\examples{
data(chm_chablais3)
chm_chablais3 <- terra::rast(chm_chablais3)

# tree segmentation
segments <- tree_segmentation(chm_chablais3)

# tree extraction
trees <- tree_extraction(segments, crown = TRUE)
# create crown polygons from WKT field
trees_crowns <- sf::st_as_sf(sf::st_drop_geometry(trees), wkt = "crown")

# summary of trees without wkt field
summary(trees[, -which(names(trees) == "crown")])

# plot initial image
terra::plot(chm_chablais3, main = "CHM and extracted trees")

# add treetop positions
plot(trees["h"], add = TRUE, cex = trees$h/20, col = "black")
# add crowns
plot(sf::st_geometry(trees_crowns), add = TRUE, border = "black", col = NA)

# plot segments
terra::plot(segments$segments_id, main = "Segments")
# add crowns
plot(sf::st_geometry(trees_crowns), add = TRUE, border = "black", col = NA)

}
\seealso{
\code{\link{tree_segmentation}}, \code{\link{tree_detection}}
}
