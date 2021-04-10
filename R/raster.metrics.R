# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
################################
#' Computes metrics by aggregating a raster at lower resolution or summarizing attributes based on XY locations
#' 
#' Compute statistics by aggregating a raster at lower resolution. Aggregation groups are larger cells, new values are computed by applying a user-specified function to original cells contained in the larger cells. Results are provided as a data.frame which also contains the XY coordinates of the larger cells.
#'
#' @param r raster object or data.frame with xy coordinates in two first columns
#' @param res numeric. Resolution of the aggregation raster, should be a multiple of r resolution if a raster is provided
#' @param fun function. Function to compute metrics in each aggregated cell from the values contained in the initial raster (use x$layer to access raster values) / data.frame (use x$colum_name to access values)
#' @param output string. indicates the class of output object "raster" or "data.frame"
#' @return a data.frame with the XY center coordinates of the aggregated cells, and the values computed with the user-specified function or a Raster object
#' @examples
#' data(chmchablais3)
#' 
#' # raster metrics from raster
#' metrics1 <- rasterMetrics(chmchablais3, res=10)
#' metrics1
#' 
#' # raster metrics from data.frame
#' n <- 1000
#' df <- data.frame(x=runif(n, 0,100), y=runif(n, 0, 100), z1=runif(n, 0,1),
#'                  z2=runif(n,10,20))
#' # compute raster metrics
#' metrics2 <- rasterMetrics(df, res=10,
#' fun=function(x){data.frame(max.z=max(x$z1), max.sum=max(x$z1+x$z2))},
#'                 output="data.frame")
#' summary(metrics2)
#' 
#' # display raster metrics
#' raster::plot(metrics1)
#' # display data.frame metrics
#' raster::plot(raster::rasterFromXYZ(metrics2))
#' @export
rasterMetrics <-
  function(r,
           res = 20,
           fun = function(x) {
             data.frame(mean = mean(x$layer), sd = stats::sd(x$layer))
           },
           output = "raster")
    
{
  if (is.element(class(r), c("RasterLayer", "RasterBrick", "RasterStack")))
  {
    # convert to data.frame
    st <- as.data.frame(raster::rasterToPoints(r))
    projinfo <- r@crs
  } else {
    if (class(r) == "SpatialPointsDataFrame")
    {
      st <- cbind(r@coords, r@data)
      projinfo <- r@proj4string
    } else {
      st <- r
      projinfo <- NA
    }
  }
  # compute coordinates of new cell center at metrics resolution
  st$X <- round((st[,1]-res/2)/res)*res+res/2
  st$Y <- round((st[,2]-res/2)/res)*res+res/2
  # compute metrics by grouping factor
  dummy <- lapply(split(st, list(st$X, st$Y), sep="_"), FUN=fun)
  # convert to data.frame
  dummy <- as.data.frame(do.call(rbind, dummy))
  # 
  dummy <- cbind(matrix(as.numeric(unlist(strsplit(row.names(dummy),"_"))),ncol=2,byrow=TRUE), dummy)
  names(dummy)[1:2] <- c("X", "Y")
  # add id column or coordinates
  if(output=="raster")
  {
    raster::rasterFromXYZ(dummy, res=c(res, res), crs=projinfo)
  } else {dummy}
}