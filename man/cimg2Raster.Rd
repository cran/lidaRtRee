% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_detection.R
\name{cimg2Raster}
\alias{cimg2Raster}
\title{Cimg to SpatRaster conversion}
\usage{
cimg2Raster(cimg, r = NULL)
}
\arguments{
\item{cimg}{raster object. raster of canopy height model, preferably filtered to avoid effect of holes on volume and surface computation}

\item{r}{SpatRaster object. defines the extent and projection of conversion result}
}
\value{
A SpatRaster object
}
\description{
converts a cimg object to a SpatRaster object
}
\examples{
data(chm_chablais3)
chm_chablais3 <- terra::rast(chm_chablais3)

# convert raster to cimg object
chm_cim <- raster2Cimg(chm_chablais3)

# apply filtering
chm_cim_filt <- dem_filtering(chm_cim,
  nl_filter = "Closing",
  nl_size = 3,
  sigma = 0
)$non_linear_image

# convert to SpatRaster
chm_filt <- cimg2Raster(chm_cim_filt, chm_chablais3)

# plot SpatRaster
terra::plot(chm_chablais3)

# plot cimg object
plot(chm_cim)

# plot filtered cimg object
plot(chm_cim_filt)

# plot filtered SpatRaster
terra::plot(chm_filt)
}
\seealso{
\code{\link{raster2Cimg}}
}
