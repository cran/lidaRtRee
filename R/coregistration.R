# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
########################### FUNCTIONS FOR COREGISTRATION #############################
#
######################################
#' Raster corresponding to circle extent
#' 
#' Creates an empty raster which extents corresponds to the circle specified by center coordinates, radius and optional buffer size.
#'
#' @param X numeric. easting coordinate of plot center in meters
#' @param Y numeric. northing coordinate of plot center in meters
#' @param radius numeric. plot radius in meters
#' @param resolution numeric. raster resolution in meters
#' @param buffer numeric. buffer to be added to plot radius in meters
#' @return A raster object
#' @examples
#' circle2Raster(100, 100, 20, 1, 5)
#' @export
circle2Raster <- function(X, Y, radius, resolution=0.5, buffer=0.5)
{
  r <- raster::raster()
  r@extent@xmin <- floor((X-radius-buffer)/resolution)*resolution
  r@extent@ymin <- floor((Y-radius-buffer)/resolution)*resolution
  r@extent@xmax <- ceiling((X+radius+buffer)/resolution)*resolution
  r@extent@ymax <- ceiling((Y+radius+buffer)/resolution)*resolution
  raster::res(r) <- resolution
  raster::crs(r) <- NA
  r
}
######################################
#' Correlation between two rasters
#' 
#' computes correlation between two rasters, based on the extent of the smallest one.
#'
#' @param raster.b raster. raster to correlate with largest extent
#' @param raster.s raster. raster to correlate with smallest extent
#' @param mask raster. mask of area to correlate
#' @param small.SC boolean. is the small raster already standardized and centered ?
#' @seealso \code{\link{rastersMovingCor}} to compute correlation between rasters for different translations
#' @return A numeric
#' @examples
#' # create raster
#' r.b <- raster::raster()
#' raster::extent(r.b) <- c(0,40,0,40)
#' raster::res(r.b) <- 1
#' xy <- raster::xyFromCell(r.b,1:length(r.b))
#' 
#' # add Gaussian surface and noise
#' z <- 3*exp(-((xy[,1]-20)^2+(xy[,2]-20)^2/2)/6)
#' r.b <- raster::rasterFromXYZ(cbind(xy,z))
#' 
#' # create circular mask of radius 5
#' z.mask <- (xy[,1]-20)^2 + (xy[,2]-20)^2 < 5^2
#' r.mask <- raster::rasterFromXYZ(cbind(xy, z.mask))
#' 
#' # create small raster of size 20
#' r.s <- raster::crop(r.b, raster::extent(c(10,30,10,30)))
#' 
#' # add noise to small raster
#' raster::values(r.s) <- raster::values(r.s) + rnorm(length(r.s),0,0.5)
#' r.mask <- raster::crop(r.mask, raster::extent(c(10,30,10,30)))
#'
#' # compute correlation on masked area where signal to noise ratio is lower
#' rasters2Cor(r.b, r.s, r.mask, small.SC=FALSE)
#' 
#' # compute correlation for whole small raster
#' rasters2Cor(r.b, r.s, small.SC=FALSE)
#' 
#' # display large raster
#' raster::plot(r.b, main="Large raster")
#' # display small raster
#' raster::plot(r.s, main="Small raster")
#' # display mask
#' raster::plot(r.mask, main="Computation mask")
#' @export
rasters2Cor <- function(raster.b, raster.s, mask=NULL, small.SC=TRUE)
{
  # crop large raster to smaller size
  raster.b <- raster::crop(raster.b,raster.s)
  # apply raster mask if existing
  if (!is.null(mask))
  {
    # apply raster mask to small raster if not already centered / standardized
    if (small.SC==F) {raster.s <- raster.s*mask}
    raster.b <- raster.b*mask
  }
  # center small raster if not already done
  if (small.SC==F)
  {
    # compute standard deviation
    sd.raster.s <- stats::sd(raster::values(raster.s),na.rm=TRUE)
    # compute mean
    m.raster.s <- mean(raster::values(raster.s),na.rm=TRUE)
    # center values of small raster
    raster.s <- raster.s - m.raster.s
  } else {sd.raster.s <- 1}
  # compute standard deviation of large raster
  sd.raster.b <- stats::sd(raster::values(raster.b),na.rm=TRUE)
  # compute mean of large raster
  m.raster.b <- mean(raster::values(raster.b),na.rm=TRUE)
  # center large raster
  raster.b <- raster.b-m.raster.b
  #
  # compute correlation value
  mean(raster::values(raster.b*raster.s),na.rm=TRUE)/(sd.raster.s*sd.raster.b)
}
#
######################################
#' Correlation between rasters for different XY translations
#' 
#' computes correlation between two rasters for different XY translations. The correlation values are computed on the extent of the smallest raster using \code{\link{rasters2Cor}}, after applying an optional mask, and for each translation within a buffer area.  
#'
#' @param raster.b raster. raster to correlate with largest extent
#' @param raster.s raster. raster to correlate with smallest extent
#' @param mask raster. mask of area to correlate, applied to small raster
#' @param buffer numeric. radius of the circular buffer area for possible translations
#' @param step numeric. increment step of translations within buffer area to compute correlation values, should be a multiple of raster resolution
#' @seealso \code{\link{rasterLocalmax}} to extract local maximum of resulting correlation raster, \code{\link{rasters2Cor}}
#' @return A raster. Raster value at coordinates x,y correspond to the correlation between the large raster and the small raster when small raster center has been translated of (x,y)
#' @examples
#' # create raster
#' r.b <- raster::raster()
#' raster::extent(r.b) <- c(0,40,0,40)
#' raster::res(r.b) <- 1
#' xy <- raster::xyFromCell(r.b,1:length(r.b))
#' 
#' # add Gaussian surfaces
#' z1 <- 1.5*exp(-((xy[,1]-22)^2+(xy[,2]-22)^2/2)/5)
#' z2 <- exp(-((xy[,1]-20)^2+(xy[,2]-22)^2/2)/3)
#' z3 <- 1.5*exp(-((xy[,1]-17)^2+(xy[,2]-17)^2/2)/5)
#' r.b <- raster::rasterFromXYZ(cbind(xy,z1+z2+z3))
#' 
#' # create small raster
#' r.s <- raster::crop(r.b, raster::extent(c(15,25,15,25)))
#' # offset raster by (-2, -2)
#' raster::extent(r.s) <- c(13,23,13,23)
#'
#' # compute correlations for translations inside buffer
#' rr <- rastersMovingCor(r.b, r.s, buffer=6, step=1)
#' rr
#' 
#' # display large raster
#' raster::plot(r.b, main="Large raster")
#' # display small raster
#' raster::plot(r.s, main="Small raster")
#' # display correlation
#' raster::plot(rr, xlab="X translation", ylab="Y translation",
#' main="Correlation between rasters")
#' @export
rastersMovingCor <- function(raster.b, raster.s, mask=NULL, buffer=19, step=0.5)
{
  # create data.frame for results
  resultat <- data.frame(xoffset=numeric(),yoffset=numeric(),correlation=numeric())
  # initialize output line
  n <- 0
  # perform computation on small raster
  # apply mask
  if (!is.null(mask))
  {
    raster.s <- raster.s*mask
  }
  # center and standardize small raster
  raster.s <- raster.s-mean(raster::values(raster.s),na.rm=TRUE)
  raster.s <- raster.s/stats::sd(raster::values(raster.s),na.rm=TRUE)
  #
  # compute squared radius
  buff2 <- buffer^2
  # copy small raster before offsetting
  matrice2 <- raster.s
  # copy mask before offsetting
  if (!is.null(mask)) {mask2 <- mask} else {mask2 <- NULL}
  # loop on X offset        
  for (xoffset in seq(from=-buffer,to=buffer,by=step))
  {
    # modify small raster X extent
    matrice2@extent@xmin <- raster.s@extent@xmin+xoffset
    matrice2@extent@xmax <- raster.s@extent@xmax+xoffset
    # modify mask X extent
    if (!is.null(mask))
    {
      mask2@extent@xmin <- mask@extent@xmin+xoffset
      mask2@extent@xmax <- mask@extent@xmax+xoffset
    }
    # loop on Y offset
    for (yoffset in seq(from=-buffer,to=buffer,by=step))
    {
      # check if new position is inside the search buffer
      if (yoffset^2+xoffset^2<=buff2)
      {
        # modify small raster Y extent
        matrice2@extent@ymin <- raster.s@extent@ymin+yoffset
        matrice2@extent@ymax <- raster.s@extent@ymax+yoffset
        # modify mask Y extent
        if (!is.null(mask))
        {
          mask2@extent@ymin <- mask@extent@ymin+yoffset
          mask2@extent@ymax <- mask@extent@ymax+yoffset
        }
        # increment line
        n <- n+1
        # compute correlation and record results
        resultat[n,] <- c(xoffset, yoffset, rasters2Cor(raster.b,matrice2,mask2,T))
      }
    }
  }
  # convert results to raster
  taille <- length(seq(from=-buffer,to=buffer,by=step))
  rr <- raster::raster(nrows=taille,ncols=taille,xmx=buffer+step/2,xmn=-buffer-step/2,ymx=buffer+step/2,ymn=-buffer-step/2)
  rr <- raster::rasterize(resultat[,1:2],rr,resultat[,3],fun=max)
  rr
}
########################################################
#' Statistics of raster local maximum
#' 
#' identifies global maximum and second global maximum from raster (e.g. output from \code{\link{rastersMovingCor}}), and computes related statistics. Local maxima can be excluded based on a minimum distance \code{dm} to nearest local maximum.
#'
#' @param r raster. typically output of \code{\link{rastersMovingCor}}
#' @param dm numeric. minimum distance between two local maxima in meters
#' @param med1 numeric. window radius to compute median value around the maximum position (default: 1m)
#' @param med2 numeric. window radius #2 to compute median value around the maximum position (default: 2m)
#' @param quanta numeric. quantile value to compute for raster values (default: 3rd quartile)
#' @param quantb numeric. quantile #2 value to compute for raster values (default: median)
#' @seealso \code{\link{rastersMovingCor}}, \code{\link{coregistration}} for application to the coregistration of tree inventory data with canopy height models
#' @return A data.frame with value of maximum, position of maximum, position of second maximum, ratio of max value to 2nd max, ratio of max value to median of neighborhood (size1 and size 2), ratio of max value to raster quantiles 1 and 2
#' @examples
#' # create raster
#' r.b <- raster::raster()
#' raster::extent(r.b) <- c(0,40,0,40)
#' raster::res(r.b) <- 1
#' xy <- raster::xyFromCell(r.b,1:length(r.b))
#' 
#' # add Gaussian surfaces
#' z1 <- 1.5*exp(-((xy[,1]-22)^2+(xy[,2]-22)^2/2)/5)
#' z2 <- exp(-((xy[,1]-20)^2+(xy[,2]-22)^2/2)/3)
#' z3 <- 1.5*exp(-((xy[,1]-17)^2+(xy[,2]-17)^2/2)/5)
#' r.b <- raster::rasterFromXYZ(cbind(xy,z1+z2+z3))
#' 
#' # create small raster
#' r.s <- raster::crop(r.b, raster::extent(c(15,25,15,25)))
#' # offset raster by (-2, -2)
#' raster::extent(r.s) <- c(13,23,13,23)
#'
#' rr <- rastersMovingCor(r.b, r.s, buffer=6, step=1)
#' loc.max <- rasterLocalmax(rr)
#' loc.max
#' 
#' # plot raster
#' raster::plot(rr)
#' # add location of two local maxima
#' graphics::points(loc.max[1,c("dx1","dx2")], loc.max[1,c("dy1","dy2")],
#'                  cex=c(1,0.5), pch=3)
#' @export
rasterLocalmax <- function(r, dm=2, med1=1, med2=2, quanta=0.75, quantb=0.5)
{
  # check if raster is valid
  if (length(r)>1 && max(raster::values(r),na.rm=TRUE)>-Inf)
  {
    # replace NA values with -Inf
    dummy <- r
    dummy[is.na(dummy)] <- -Inf
    # maxima detection
    maxi <- maximaDetection(dummy, jitter=FALSE)
    dummy[] <- 1
    # maxima selection, without value criteria
    maxi <- maximaSelection(maxi, dummy, hmin=0, dmin=dm, dprop=0)
    # set back values in maxi to correlation values
    maxi <- (maxi > 0) * r
    # indice of global max
    pixel1 <- raster::which.max(maxi)
    # max value
    max1 <- r[pixel1[1]]
    # coordinates of global max
    coord1 <- raster::xyFromCell(maxi,pixel1[1])
    # erase first max
    maxi[pixel1] <- NA
    # second maximum
    pixel2 <- raster::which.max(maxi)
    max2 <- r[pixel2[1]]
    coord2 <- raster::xyFromCell(maxi,pixel2[1])
    # median around global max
    # half resolution
    hres <- raster::res(r)[1]/2
    # median around maximum, neighborhood 1
    medloc1 <- stats::median(raster::values(raster::crop(r,raster::extent(c(coord1[1]-med1-hres, 
                                                       coord1[1]+med1+hres,
                                                       coord1[2]-med1-hres,
                                                       coord1[2]+med1+hres)))), na.rm=TRUE)
    # median around maximum, neighborhood 2
    medloc2 <- stats::median(raster::values(raster::crop(r,raster::extent(c(coord1[1]-med2-hres, 
                                                         coord1[1]+med2+hres,
                                                         coord1[2]-med2-hres,
                                                         coord1[2]+med2+hres)))), na.rm=TRUE)
    # stats of correlation raster
    quant.a <- stats::quantile(r,quanta,na.rm=TRUE)
    quant.b <- stats::quantile(r,quantb,na.rm=TRUE)
    # store results into data.frame
    data.frame(max1=max1,dx1=coord1[1],dy1=coord1[2],dx2=coord2[1],dy2=coord2[2],ratiomax1max2=max1/max2,rmedloc1=max1/medloc1,rmedloc2=max1/medloc2,rquanta=max1/quant.a,rquantb=max1/quant.b)
  } else { NA }
}
########################################################
#' Tree inventory and canopy height model coregistration
#' 
#' Computes the correlation between the canopy height model and a virtual canopy height model simulated from tree locations, for different translations of tree inventory positions, and outputs the translation corresponding to best estimated co-registration.
#' 
#' @param chm raster. canopy height model
#' @param trees data.frame. the first two columns contain xy coordinates, and the third is the value to correlate to the chm (e.g. tree heights or diameters)
#' @param mask raster. raster mask of tree inventory area
#' @param buffer numeric. radius of the circular buffer area of possible translations
#' @param step numeric. increment step of translations within buffer area to compute correlation values, should be a multiple of raster resolution
#' @param dm numeric. minimum distance between two local maxima in meters
#' @param plot boolean. whether to display the results or not
#' @seealso \code{\link{rastersMovingCor}}, \code{\link{rasterLocalmax}}
#' @return A list with two elements : first the correlation raster returned by \code{\link{rastersMovingCor}}, second a data.frame returned by \code{\link{rasterLocalmax}}
#' @references Monnet, J.-M. and Mermin, E. 2014. Cross-Correlation of Diameter Measures for the Co-Registration of Forest Inventory Plots with Airborne Laser Scanning Data. Forests 2014, 5(9), 2307-2326, \doi{10.3390/f5092307}
#' @examples
#' # tree inventory
#' trees <- data.frame(x=c(22.2, 18.3, 18.1), y =c(22.1, 22.7, 18.4), z=c(15,10,15))
#' 
#' # mask of inventory area
#' # empty raster with extent
#' tree.mask <- circle2Raster(20, 20, 9, resolution=1)
#' # fill binary mask
#' tree.mask <- rasterXYMask(rbind(c(20, 20), c(20, 20)), c(9,9), tree.mask, binary=TRUE)
#' 
#' # simulate chm raster
#' chm <- raster::raster()
#' raster::extent(chm) <- c(0,40,0,40)
#' raster::res(chm) <- 1
#' xy <- raster::xyFromCell(chm,1:length(chm))
#' 
#' # add Gaussian surfaces to simulate tree crowns
#' z1 <- trees$z[1]*exp(-((xy[,1]-trees$x[1])^2+(xy[,2]-trees$y[1])^2/2)*trees$z[1]/50)
#' z2 <- trees$z[2]*exp(-((xy[,1]-trees$x[2])^2+(xy[,2]-trees$y[2])^2/2)*trees$z[2]/50)
#' z3 <- trees$z[3]*exp(-((xy[,1]-trees$x[3])^2+(xy[,2]-trees$y[3])^2/2)*trees$z[3]/50)
#' chm <- raster::rasterFromXYZ(cbind(xy,pmax(z1,z2,z3)))#+rnorm(length(z1),0,1)))
#' 
#' # translate trees
#' trees$x <- trees$x+1
#' trees$y <- trees$y+2
#' 
#' coreg <- coregistration(chm, trees, mask=tree.mask, buffer=5, step=1, dm=1, plot=FALSE)
#' coreg$local.max[,c("dx1", "dy1")]
#' 
#' # plot raster
#' raster::plot(coreg$correlation.raster)
#' abline(h=0, lty=2)
#' abline(v=0, lty=2)
#' # add location of two local maxima
#' graphics::points(coreg$local.max[1,c("dx1","dx2")],
#' coreg$local.max[1,c("dy1","dy2")], cex=c(1,0.5), pch=3, col="red")
#' @export
coregistration <- function(chm, trees, mask= NULL, buffer=19, step=0.5, dm=2, plot=TRUE)
{
  # compute virtual chm from tree inventory
  virtualchm <- raster::rasterize(trees[,1:2], mask, trees[,3],fun=max)
  # fill NA values with 0
  virtualchm[is.na(virtualchm)] <- 0
  #
  # compute correlation raster 
  r.correlation <- rastersMovingCor(chm, virtualchm, mask=mask, buffer=buffer, step=step)
  # analyse correlation image
  resul <- rasterLocalmax(r.correlation)
  # display results
  if (plot)
  {
    # CHM
    raster::plot(chm, asp=1)
    # display initial tree positions
    graphics::points(trees[,1], trees[,2], cex=trees[,3]/40)
    # display coregistered tree positions
    graphics::points(trees[,1]+resul$dx1, trees[,2]+resul$dy1, cex=trees[,3]/40, col="red")
    graphics::legend("topleft", c("Initial", "Coregistered"), pch=1, col=c("black", "red"))
  }
  list(correlation.raster=r.correlation, local.max=resul)
}
# 