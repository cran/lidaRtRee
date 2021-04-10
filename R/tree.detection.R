# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
########################### FUNCTIONS FOR TREE DETECTION #############################
#
######################################
#' Disk-shaped matrix mask
#' 
#' Creates a matrix with TRUE values shaping a centered disk
#'
#' @param width numeric. disk width in pixels, should be an uneven number
#' @return A matrix with 1 for pixels inside the disk, 0 outside
#' @examples
#' createDisk(7)
#' @export
createDisk <- function(width=5)
{
  if (width %% 2 != 1) {stop("Mask width should be uneven")}
  # matrix of row indices
  row.mat <- matrix(rep(1:width,width), nrow=width, ncol=width)
  # matrix of column indices
  col.mat <- matrix(rep(1:width,width), nrow=width, ncol=width, byrow=TRUE)
  radius <- width %/%2
  row.mat <- row.mat - radius -1
  col.mat <- col.mat - radius -1
  mask <- (col.mat^2 + row.mat^2) <= radius ^ 2
  mask
}
######################################
#' Image pre-processing (non-linear filtering and Gaussian smoothing)
#'
#' applies two filters to an image:
#' \enumerate{
#' \item A non-linear filter: closing (\code{\link[imager]{mclosing}}) with disk kernel, or median (\code{\link[imager]{medianblur}}) with square kernel
#' \item A 2D Gaussian smoother (The \code{\link[imager]{deriche}} filter is applied on both dimensions). Value-dependent smoothing is possible
#' }
#'
#' @param dem cimg object (e.g. obtained with \code{\link[imager]{as.cimg}}) or Raster Layer object (e.g. obtained with \code{\link[raster]{raster}})
#' @param nlFilter string. type of non-linear filter to apply: "None", "Closing" or "Median"
#' @param nlSize numeric. kernel width in pixel for non-linear filtering
#' @param sigmap numeric or matrix. if a single number is provided, sigmap is the standard deviation of the Gaussian filter in pixel, 0 corresponds to no smoothing. In case of matrix, the first column corresponds to the standard deviation of the filter, and the second to thresholds for image values (e.g. a filter of standard deviation specified in line \code{i} is applied to pixels in image which values are between thresholds indicated in lines \code{i} and \code{i+1}). Threshold values should be ordered in increasing order.
#' @param padding boolean. Whether image should be padded by duplicating edge values before filtering to avoid border effects
#' @examples
#' data(chmchablais3)
#' 
#' # filtering with median and Gaussian smoothing
#' im <- demFiltering(chmchablais3, nlFilter="Median", nlSize=3, sigmap=0.8)
#' 
#' # filtering with median filter and value-dependent Gaussian smoothing
#' # (less smoothing for values between 0 and 15)
#' im2 <- demFiltering(chmchablais3, nlFilter="Median", nlSize=3,
#'                     sigmap=cbind(c(0.2,0.8), c(0,15)))
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Initial image")
#' 
#' # plot image after median filter
#' raster::plot(im$non.linear.image, main="Median filter")
#' 
#' # plot image after median and Gaussian filters
#' raster::plot(im$smoothed.image, main="Smoothed image")
#' 
#' # plot image after median and value-dependent Gaussian filters
#' raster::plot(im2$smoothed.image, main="Value-dependent smoothing")
#' @seealso \code{\link{maximaDetection}}, filters of imager package: \code{\link[imager]{mclosing}}, \code{\link[imager]{medianblur}}, \code{\link[imager]{deriche}}
#' @return A list of two cimg or a RasterStack objects: image after non-linear filter and image after both filters
#' @export
demFiltering <- function(dem, nlFilter="Closing", nlSize=5, sigmap=0.3, padding=TRUE)
{
  # convert rasterLayer to cimg object if necessary
  if (class(dem)[1] =="RasterLayer")
  {
    dem.c <- raster2Cimg(dem)
  } else {
    dem.c <- dem
  }
  #
  if (padding)
  {
    # padding number of cells is maximum of half width of non linear filter or ceiling value of three times sigmap
    if (!is.null(dim(sigmap)))
    {dummy <- max(sigmap[,1])} else {dummy <- sigmap}
    border.size <- max((nlSize-1)/2+1, ceiling(dummy*3))
    # convert to matrix
    dem.c <- as.matrix(dem.c)
    # duplicate columns
    dem.c <- dem.c[,c(rep.int(1,border.size), 1:ncol(dem.c), rep.int(ncol(dem.c), border.size))]
    # duplicate rows
    dem.c <- dem.c[c(rep.int(1,border.size), 1:nrow(dem.c), rep.int(nrow(dem.c), border.size)),]
    dem.c <- imager::as.cimg(dem.c)
  }
  #
  # non-linear filtering
  dem.nl <- dem.c
  if (nlFilter=="Closing")
  {
    # closing with structuring element (disk)
    strel <- imager::as.cimg(createDisk(nlSize))
    dem.nl <- imager::mclosing(dem.nl,strel)
  }
  if (nlFilter=="Median")
  {
    # median on square window
    # strel <- imager::as.cimg(createDisk(nlSize))
    dem.nl <- imager::medianblur(dem.nl,nlSize)
  }
  #
  # linear filtering
  # gaussian smoothing
  # if several values of sigma are provided, value-dependent smoothin is performed
  if (length(sigmap)>1)
  {
    dem.gs <- dem.nl
    # for each value of sigma
    for (i in 1:nrow(sigmap))
    {
      # perform 2D smoothing
      dummy <- imager::deriche(dem.nl, sigmap[i,1], axis="x")
      dummy <- imager::deriche(dummy, sigmap[i,1], axis="y")
      # identify pixels which values are superior to the threshold
      temp <- dem.gs>=sigmap[i,2]
      # update only those values in the output
      dem.gs[temp] <- dummy[temp]
    }
  } else { # if only one value of sigma is provided
    if (sigmap>0)
    {
      dem.gs <- imager::deriche(dem.nl, sigmap, axis="x")
      dem.gs <- imager::deriche(dem.gs, sigmap, axis="y")
    } else { # if 0, no smoothing is performed
      dem.gs <- dem.nl
    }
  }
  # remove padding
  if (padding)
  {
    dem.nl <- as.matrix(dem.nl)
    dem.nl <- imager::as.cimg(dem.nl[(border.size+1):(nrow(dem.nl)-border.size), (border.size+1):(ncol(dem.nl)-border.size)])
    dem.gs <- as.matrix(dem.gs)
    dem.gs <- imager::as.cimg(dem.gs[(border.size+1):(nrow(dem.gs)-border.size), (border.size+1):(ncol(dem.gs)-border.size)])
  }
  # convert cimg objects to rasterLayer if necessary
  if (class(dem)[1] =="RasterLayer")
  {
    output <- raster::addLayer(cimg2Raster(dem.nl, dem), cimg2Raster(dem.gs, dem))
    names(output) <- c("non.linear.image", "smoothed.image")
  } else {
    output <- list(non.linear.image=dem.nl, smoothed.image=dem.gs)
  }
  output
}
###################################
#' Local maxima extraction on image
#' 
#' Variable window size maxima detection is performed on the image to extract local maxima position and calculate the window size where they are global maxima. Gaussian white noise is added to the image to avoid adjacent maxima due to neighbor pixels with identical value.
#'
#' @param dem cimg object (e.g. as created by \code{\link[imager]{cimg}}) or RasterLayer object (e.g. obtained with \code{\link[raster]{raster}})
#' @param dem.res numeric. image resolution, in case \code{dem} is a rasterLayer object, \code{dem.res} is extracted from the object by \code{\link[raster]{res}}
#' @param max.width numeric. maximum kernel width in pixel to check for local maximum
#' @param jitter boolean. indicates if noise should be added to image values to avoid the adjacent maxima due to the adjacent pixels with equal values
#' @return A cimg object or RasterLayer object which values are the radius (n) in meter of the square window (width 2n+1) where the center pixel is global maximum
#' @examples
#' data(chmchablais3)
#' 
#' # maxima detection
#' maxi <- maximaDetection(chmchablais3)
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Initial image")
#' 
#' # plot maxima image
#' raster::plot(maxi, main="Local maxima")
#' @seealso \code{\link{demFiltering}}, \code{\link{maximaSelection}}
#' @export
maximaDetection <- function(dem, dem.res=1, max.width=21, jitter=TRUE)
{
  # convert rasterLayer to cimg object if necessary
  if (class(dem)[1] =="RasterLayer")
  {
    dem.gs <- raster2Cimg(dem)
    dem.res <- raster::res(dem)[1]
  } else {
    dem.gs <- dem
  }
  #
  # add absolute of gaussian white noise, mean=0, sd=sd(dem.gs)/100000 to non 0 pixels
  if (jitter) {dem.gs <- dem.gs + abs(imager::imnoise(dim=dim(dem.gs), sd=stats::sd(dem.gs)/100000)) * (dem.gs!=0)}
  # 
  max.radius <- max.width %/% 2
  # extraction of maxima on variable window size (from 0 to max.radius)
  maxi <- NULL
  for (i in 1:max.radius)
  {
    # create square structuring element for dilation
    strel <- imager::imfill(2*i+1,2*i+1,val=1)
    # if initialization check where pixel values are equal in original image and in image dilated by structuring element
    if (is.null(maxi))
    {
      maxi <- imager::as.cimg(dem.gs == imager::dilate(dem.gs,strel))*i
      # otherwise perform check and update previous result
    } else {
      maxi <- imager::parmax(list(maxi, imager::as.cimg(dem.gs == imager::dilate(dem.gs,strel))*i))
    }
  }
  # convert window size from pixels to meters
  maxi[maxi>0] <- (maxi[maxi>0] + 1)*dem.res
  # convert cimg objects to rasterLayer if necessary
  if (class(dem)[1] =="RasterLayer")
  {
    maxi <- cimg2Raster(maxi, dem)
  }
  maxi
}
################################
#' Image maxima selection based on values and neighborhood of local maxima
#' 
#' In a maxima image (output of \code{\link{maximaDetection}}), sets values to zero for pixels which
#' \enumerate{
#' \item value in the initial image (from which maxima were detected) are below a threshold 
#' \item values in the maxima image (corresponding to the radius of the neighborhood where they are global maxima) are below a threshold depending on the initial image value.
#' }
#'
#' @param maxi cimg object or RasterLayer object. image with local maxima (typically output from \code{\link{maximaDetection}}, image values correspond to neighborhood radius on which pixels are global maxima in the initial image)
#' @param dem.nl cimg object. initial image from which maxima were detected
#' @param hmin numeric. minimum value in initial image for a maximum to be selected
#' @param dmin numeric. intercept term for selection of maxima depending on neighborhood radius: \code{maxi >= dmin + dem.nl * dprop}
#' @param dprop numeric. proportional term for selection of maxima depending on neighborhood radius: \code{maxi >= dmin + dem.nl * dprop}
#' @return A cimg object or rasterLayer object which values are the radius (n) in meter of the square window (width 2n+1) where the center pixel is global maximum and which fulfill the selection criteria
#' @examples
#' data(chmchablais3)
#' 
#' # maxima detection
#' maxi <- maximaDetection(chmchablais3)
#' 
#' # several maxima selection settings
#' selected.maxi.hmin <- maximaSelection(maxi, chmchablais3, hmin=15)
#' selected.maxi.dm <- maximaSelection(maxi, chmchablais3, dm=2.5)
#' selected.maxi <- maximaSelection(maxi, chmchablais3, dm=1, dprop=0.1)
#' 
#' # corresponding count number of remaining maxima
#' table(raster::values(maxi))
#' table(raster::values(selected.maxi.hmin))
#' table(raster::values(selected.maxi.dm))
#' table(raster::values(selected.maxi))
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Initial image")
#' 
#' # plot maxima images, original and first case
#' raster::plot(maxi, main="Local maxima")
#' raster::plot(selected.maxi, main="Selected maxima")
#' 
#' @seealso \code{\link{maximaDetection}}
#' @export
maximaSelection <- function(maxi,dem.nl,hmin=0,dmin=0,dprop=0)
{
  # convert rasterLayer to cimg object if necessary
  if (class(maxi)[1] =="RasterLayer")
  {
    israster <- TRUE
    dem <- maxi
    maxi <- raster2Cimg(maxi)
    dem.nl <- raster2Cimg(dem.nl)
  } else {
    israster <- FALSE
  }
  # height filter
  maxi[dem.nl < hmin] <- 0
  # distance filter
  maxi[maxi < (dmin + dem.nl * dprop)] <- 0
  # convert to rasterLayer if necessary
  if(israster)
  {
    cimg2Raster(maxi, dem)
  } else {maxi}
}
################################
#' Image segmentation by seed-based watershed algorithm
#'
#' performs a seed-based watershed segmentation (wrapper for imager::watershed)
#'
#' @param maxi cimg or rasterLayer object. image with seed points (e.g. from \code{\link{maximaDetection}} or \code{\link{maximaSelection}})
#' @param dem.nl cimg or rasterLayer object. image for seed propagation (typically initial image used for maxima detection).
#' @return A cimg object or rasterlayer object with segments id
#' @examples
#' data(chmchablais3)
#'
#' # median filter
#' chmchablais3 <- demFiltering(chmchablais3, nlFilter="Median", nlSize=3,
#'                              sigmap=0)$non.linear.image
#' 
#' # maxima detection
#' maxi <- maximaDetection(chmchablais3)
#' 
#' # maxima selection
#' selected.maxi <- maximaSelection(maxi, chmchablais3, dm=1, dprop=0.1)
#' 
#' # segmentation
#' seg.maxi <- segmentation(maxi, chmchablais3)
#' seg.selected.maxi <- segmentation(selected.maxi, chmchablais3)
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Median filter")
#' 
#' # plot segmented image
#' # replace segment with id 0 (not a tree) with NA
#' seg.maxi[seg.maxi==0] <- NA
#' raster::plot(seg.maxi %% 8, main="Segments, no maxima selection", col=rainbow(8))
#' seg.selected.maxi[seg.selected.maxi==0] <- NA
#' raster::plot(seg.selected.maxi %% 8, main="Segments, maxima selection", col=rainbow(8))
#' 
#' @seealso \code{\link{maximaDetection}}, \code{\link{maximaSelection}}, \code{\link{segAdjust}}
#' @export
segmentation <- function(maxi,dem.nl)
{
  # convert rasterLayer to cimg object if necessary
  if (class(maxi)[1] =="RasterLayer")
  {
    israster <- TRUE
    dem <- maxi
    maxi <- raster2Cimg(maxi)
    dem.nl <- raster2Cimg(dem.nl)
  } else {
    israster <- FALSE
  }
  # maxima locations
  a <- which(maxi>0)
  # tree number
  na <- length(a)
  # create seed image
  dummy <- maxi
  # initialise seeds with id
  dummy[a] <- sample(1:na, na)
  # watershed segmentation
  dem.w <- imager::watershed(dummy,dem.nl)
  #
  # convert to raster if necessary
  if (israster)
  {
    cimg2Raster(dem.w, dem)
  } else {
    dem.w
  }
}
################################
#' Image statistic in segment
#'
#' compute zonal statistic of an image
#'
#' @param segms cimg or rasterLayer object. image with segments id (e.g. from \code{\link{segmentation}})
#' @param dem.nl cimg or rasterLayer object. image to compute statistic from
#' @param fun function to compute statistis from values in each segment
#' @return A cimg object or raster object with values of the statistic
#' @examples
#' data(chmchablais3)
#'
#' # median filter
#' chmchablais3 <- demFiltering(chmchablais3, nlFilter="Median", nlSize=3,
#'                              sigmap=0)$non.linear.image
#' 
#' # maxima detection
#' maxi <- maximaDetection(chmchablais3)
#' 
#' # segmentation
#' seg.maxi <- segmentation(maxi, chmchablais3)
#' 
#' # compute image of maximum value in each segment
#' max.in.segment <- rasterZonalStats(seg.maxi, chmchablais3)
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Median filter")
#' 
#' # plot segments and image of max value inside segments
#' seg.maxi[seg.maxi==0] <- NA
#' raster::plot(seg.maxi %% 8, main="Segments", col=rainbow(8))
#' raster::plot(max.in.segment, main="Max value in segment")
#' @seealso \code{\link{segmentation}}
#' @export
rasterZonalStats <- function(segms,dem.nl,fun=max)
{
  # convert rasterLayer to cimg object if necessary
  if (class(segms)[1] =="RasterLayer")
  {
    israster <- TRUE
    dem <- segms
    segms <- raster2Cimg(segms)
    dem.nl <- raster2Cimg(dem.nl)
  } else {
    israster <- FALSE
  }
  dummy <- stats::aggregate(as.numeric(dem.nl), by=list(as.numeric(segms)), FUN=fun)
  dem.wh <- dummy$x[match(as.numeric(segms), dummy$Group.1)]
  #
  dem.wh <- imager::as.cimg(matrix(dem.wh,nrow(segms)))
  #
  # convert to raster if necessary
  if (israster)
  {
    cimg2Raster(dem.wh, dem)
  } else {
    dem.wh
  }
}
################################
#' Modification of segments based on values
#' 
#' in a segmented image, removes from segments the pixels which values in a reference image is below a certain percentage of the highest value inside the segment. Removed pixels are attributed 0 value.
#'
#' @param dem.w cimg or rasterLayer object. image with segments id, without 0 values
#' @param dem.wh cimg or rasterLayer object. image with max value inside segment
#' @param dem.nl cimg or rasterLayer object. image with initial values
#' @param prop numeric. proportional threshold for removal of pixels which initial values are lower than the max height of the segment (\code{dem.nl < prop x dem.wh})
#' @param min.value numeric. threshold for removel of pixels which initial values are lower (\code{dem.nl < min.value}) 
#' @param min.maxvalue numeric. threshold for complete removal of segments which maximum value height is smaller to the threshold (\code{dem.wh < min.maxvalue}) 
#' @return A cimg or rasterLayer object: image with modified segments. 
#' @examples
#' data(chmchablais3)
#' 
#' # median filter
#' chmchablais3 <- demFiltering(chmchablais3, nlFilter="Median", nlSize=3,
#'                              sigmap=0)$non.linear.image
#' 
#' # maxima detection and selection
#' maxi <- maximaDetection(chmchablais3)
#' selected.maxi <- maximaSelection(maxi, chmchablais3, dm=1, dprop=0.1)
#' 
#' # segmentation
#' seg.selected.maxi <- segmentation(selected.maxi, chmchablais3)
#' 
#' # max value in segments
#' max.in.segment <- rasterZonalStats(seg.selected.maxi, chmchablais3)
#' 
#' # segmentation modification
#' seg.modif1 <- segAdjust(seg.selected.maxi, max.in.segment,
#'  chmchablais3, prop=0.5)
#' seg.modif2 <- segAdjust(seg.selected.maxi, max.in.segment,
#'  chmchablais3, prop=0, min.value=5, min.maxvalue=10)
#' 
#' # plot initial segmented image
#' seg.selected.maxi[seg.selected.maxi==0] <- NA
#' raster::plot(seg.selected.maxi %% 8, main="Initial segments", col=rainbow(8))
#' seg.modif1[seg.modif1==0] <- NA
#' raster::plot(seg.modif1 %% 8, main="Modified segments 1", col=rainbow(8))
#' seg.modif2[seg.modif2==0] <- NA
#' raster::plot(seg.modif2 %% 8, main="Modified segments 2", col=rainbow(8))
#' @seealso \code{\link{maximaDetection}}, \code{\link{maximaSelection}}
#' @export
segAdjust <- function(dem.w, dem.wh, dem.nl, prop=0.3, min.value=2, min.maxvalue=5)
{
  # convert rasterLayer to cimg object if necessary
  if (class(dem.w)[1] =="RasterLayer")
  {
    israster <- TRUE
    dem <- dem.w
    dem.w <- raster2Cimg(dem.w)
    dem.wh <- raster2Cimg(dem.wh)
    dem.nl <- raster2Cimg(dem.nl)
  } else {
    israster <- FALSE
  }
  # removal of segments with maximum value lower than min.maxvalue
  dem.w[dem.wh < min.maxvalue] <- 0
  # removal of pixels with values lower than prop * max value in segment
  dem.w[dem.nl < prop * dem.wh] <- 0
  # removal of pixels with values lower than min.value
  dem.w[dem.nl < min.value] <- 0
  # convert to raster if necessary
  if (israster)
  {
    dem.w <- cimg2Raster(dem.w, dem)
  }
  dem.w
}
################################
#' Preprocessing and segmentation of raster image for tree identification
#' 
#' global function for preprocessing (filtering), maxima detection and selection, segmentation and segmentation adjustment of a raster image.
#'
#' @param dem raster object or string indicating location of raster file (typically a canopy height model or a digital surface model; in the latter case the dtm parameter should be provided)
#' @param nlFilter string. specifies the non-linear filter for image pre-processing, should be an option of function \code{\link{demFiltering}}
#' @param nlSize numeric. width of kernel of non-linear filter in pixels
#' @param sigma numeric or matrix. if a single number is provided, sigmap is the standard deviation of Gaussian filter in meters, 0 corresponds to no smoothing. In case of matrix, the first column corresponds to the standard deviation of the filter, and the second to thresholds for image values (e.g. a filter of standard deviation specified in line \code{i} is applied to pixels in image which values are between thresholds indicated in lines \code{i} and \code{i+1}). Threshold values should be ordered in increasing order.
#' @param dmin numeric. treetop minimum distance to next higher pixel in meters
#' @param dprop numeric. number defining the treetop minimum distance as proportion of height to next higher pixel
#' @param hmin numeric. minimum treetop height
#' @param crownProp numeric. minimum height of tree crown as proportion of treetop height
#' @param crownMinH numeric. minimum crown height
#' @param dtm raster object or string indicating location of raster file with the terrain model. If provided, the maxima extraction and watershed segmentation are performed on the dem (this avoids the deformation of crown because of the normalisation with terrain), but maxima selection and segment adjustement are performed on 'dem-dtm' because the selection criteria is the height to terrain.
#' @references Monnet, J.-M. 2011. Using airborne laser scanning for mountain forests mapping: Support vector regression for stand parameters estimation and unsupervised training for treetop detection. Ph.D. thesis. University of Grenoble, France. Section 6.2 \url{https://tel.archives-ouvertes.fr/tel-00652698/document}
#' 
#' Monnet, J.-M., Mermin, E., Chanussot, J., Berger, F. 2010. Tree top detection using local maxima filtering: a parameter sensitivity analysis. Silvilaser 2010, the 10th International Conference on LiDAR Applications for Assessing Forest Ecosystems, September 14-17, Freiburg, Germany, 9 p. \url{https://hal.archives-ouvertes.fr/hal-00523245/document}
#' @examples
#' data(chmchablais3)
#' 
#' # tree segmentation
#' segments <- treeSegmentation(chmchablais3)
#' segments2<- treeSegmentation(chmchablais3, nlFilter="Median", nlSize=3,
#'   sigma=cbind(c(0.2,0.8), c(0,15)),dmin=0,dprop=0,hmin=10,crownProp=0.5,crownMinH=5)
#' 
#' # plot initial image segments
#' raster::plot(chmchablais3, main="Initial image")
#' raster::plot(segments$smoothed.dem, main="Filtered image")
#' raster::plot(segments$local.maxima, main="Local maxima")
#' #
#' # replace segment with id 0 (not a tree) with NA
#' segments$segments.id[segments$segments.id==0] <- NA
#' raster::plot(segments$segments.id %% 8 , main="Segments", col=rainbow(8))
#' #
#' # plot segmentation with other parameters
#' segments2$segments.id[segments2$segments.id==0] <- NA
#' raster::plot(segments2$segments.id %% 8, main="Segments2", col=rainbow(8))
#' 
#' @seealso \code{\link{demFiltering}}, \code{\link{maximaDetection}}, \code{\link{maximaDetection}}, \code{\link{maximaSelection}}, \code{\link{segmentation}}, \code{\link{segAdjust}}, \code{\link{treeExtraction}}
#' @return A RasterStack with 4 layers: selected local maxima (values = distance to higher pixel), segments, non-linear preprocessed dem, smoothed preprocessed dem
#' @export
treeSegmentation <- function(dem, nlFilter="Closing", nlSize=5, sigma=0.3, dmin=0, dprop=0.05, hmin=5, crownProp=0.3, crownMinH=2, dtm=NULL)
{
  #
  if (crownMinH > hmin) {stop("minimum tree height lower than minimum crown height")}
  if (is.character(dem)) {dem <- raster::raster(dem)}
  # convert NAs to 0s
  if (!is.null(dtm))
  {
    if (is.character(dtm)) {dtm <- raster::raster(dtm)}
    dem[is.na(dem)] <- dtm[is.na(dem)]
  } else {
    dtm <- 0
    dem[is.na(dem)] <- dtm
    }
  #
  # conversion of Gaussian filter standard deviation from meters to pixels
  if (length(sigma==1))
  {
    sigma <- sigma/raster::res(dem)[1]
  } else {
    sigma[,2] <- sigma[,2]/raster::res(dem)[1]
  }
  # dem filtering 
  temp <- demFiltering(dem, nlFilter=nlFilter, nlSize=nlSize, sigmap=sigma)
  dem.nl <- temp$non.linear.image
  dem.gs <- temp$smoothed.image
  #
  # maxima detection
  maxi <- maximaDetection(dem.gs,raster::res(dem)[1])
  #
  # maxima selection
  maxi <- maximaSelection(maxi, dem.nl-dtm, 0, dmin, dprop) # no selection based on top height at this stage, otherwise some artefacts occur in the watershed step
  #
  # segmentation
  dem.w <- segmentation(maxi, dem.nl)
  dem.wh <- rasterZonalStats(dem.w, dem.nl-dtm, fun=max)
  #
  # segmentation adjustment
  dem.w <- segAdjust(dem.w,dem.wh,dem.nl-dtm,crownProp,crownMinH, hmin) # removal of trees with top heigth < hmin
  # remove trees from r.maxi
  maxi[dem.w==0] <- 0
  #
  output <- raster::addLayer(maxi, dem.w, dem.nl, dem.gs)
  names(output) <- c("local.maxima", "segments.id", "filled.dem", "smoothed.dem")
  output
}
################################
#' Tree extraction
#' 
#' creates a dataframe with segment id, height and coordinates of maxima, surface and volume, computed from three images: initial, local maxima and segmented, obtained with \code{\link{treeSegmentation}}
#'
#' @param r.dem.nl raster object. raster of canopy height model, preferably filtered to avoid effect of holes on volume and surface computation
#' @param r.maxi raster object. raster positive values at local maxima
#' @param r.dem.w raster object. segmented raster
#' @param r.mask raster object. only segments which maxima are inside the mask are extracted
#' @return A spatial data.frame with tree id, local maximum stats (height, dominance radius), segment stats (surface and volume).
#' @examples
#' data(chmchablais3)
#' 
#' # tree segmentation
#' segments <- treeSegmentation(chmchablais3)
#' 
#' # tree extraction
#' trees <- treeExtraction(segments$filled.dem, segments$local.maxima,
#'                         segments$segments.id)
#' trees
#' 
#' # plot initial image
#' raster::plot(chmchablais3)
#' 
#' # add treetop positions
#' sp::plot(trees, cex=trees$h/20, add=TRUE, pch=1)
#' \donttest{
#' # add segment contours (vectorization is slow)
#' contours <- raster::rasterToPolygons(segments$segments.id, dissolve=TRUE)
#' sp::plot(contours, add=TRUE, border="white")}
#'  
#' @seealso \code{\link{treeSegmentation}}
#' @export
treeExtraction <- function(r.dem.nl, r.maxi, r.dem.w, r.mask=NULL)
{
  # apply mask to remove maxima outside of region of interest
  if (!is.null(r.mask))
  {
    r.maxi <- r.maxi * r.mask
  }
  # segments surface
  s <- data.frame(t(table(raster::values(r.dem.w))))[,-1]
  names(s) <- c("id","s")
  s$s <- s$s * raster::res(r.dem.nl)[1] * raster::res(r.dem.nl)[2]
  # segments volume
  v <- stats::aggregate(raster::values(r.dem.nl), by=list(raster::values(r.dem.w)), FUN=sum)
  names(v) <- c("id","v")
  v$v <- v$v * raster::res(r.dem.nl)[1] * raster::res(r.dem.nl)[2]
  #
  # compute surface and volume inside mask, if mask is not null
  if (!is.null(r.mask))
  {
    # segments surface in mask
    sp <- data.frame(t(table(raster::values(r.dem.w*r.mask))))[,-1]
    names(sp) <- c("id","sp")
    sp$sp <- sp$sp * raster::res(r.dem.nl)[1] * raster::res(r.dem.nl)[2]
    # segments volume in mask
    vp <- stats::aggregate(raster::values(r.dem.nl*r.mask), by=list(raster::values(r.dem.w)), FUN=sum)
    names(vp) <- c("id","vp")
    vp$vp <- vp$vp * raster::res(r.dem.nl)[1] * raster::res(r.dem.nl)[2]
  }
  #
  # extract segment id, height and dom.radius
  cells <- which(raster::values(r.maxi)>0)
  # check if positive values are present
  if (length(cells)==0)
  {
    segms <- NULL
  } else {
  coord <- raster::xyFromCell(r.maxi,cells)
  segms <- data.frame(coord,id=raster::values(r.dem.w)[cells],h=raster::values(r.dem.nl)[cells],dom.radius=raster::values(r.maxi)[cells])
  segms <- merge(segms, s, all.x = TRUE)
  segms <- merge(segms, v, all.x = TRUE)
  if (!is.null(r.mask))
  {
    segms <- merge(segms, sp, all.x = TRUE)
    segms <- merge(segms, vp, all.x = TRUE)
  }
  sp::coordinates(segms) <- ~ x+y
  segms@proj4string <- r.dem.nl@crs
  }
  segms
}
#
################################
#' Cimg to RasterLayer conversion
#' 
#' converts a cimg object to a RasterLayer object
#'
#' @param cimg raster object. raster of canopy height model, preferably filtered to avoid effect of holes on volume and surface computation
#' @param rasterLayer raster object. defines the extent and projection of conversion result
#' @return A RasterLayer
#' @examples
#' data(chmchablais3)
#' 
#' # convert rasterLayer to cimg object
#' chm.cim <- raster2Cimg(chmchablais3)
#' 
#' # apply filtering
#' chm.cim.filt <- demFiltering(chm.cim,
#'                              nlFilter = "Closing",
#'                              nlSize = 3,
#'                              sigmap = 0)$non.linear.image
#' 
#' # convert to RasterLayer
#' chm.filt <- cimg2Raster(chm.cim.filt, chmchablais3)
#' 
#' # plot rasterLayer
#' raster::plot(chmchablais3)
#' 
#' # plot cimg object
#' plot(chm.cim)
#' 
#' # plot filtered cimg object
#' plot(chm.cim.filt)
#' 
#' # plot filtered rasterLayer
#' raster::plot(chm.filt)
#' @seealso \code{\link{raster2Cimg}}
#' @export
cimg2Raster <- function(cimg, rasterLayer=NULL)
{
  # convert to rasterLayer
  dem <- raster::raster(t(as.matrix(cimg)))
  # if reference is provided
  if (!is.null(rasterLayer))
  {
    # specifiy extent
    raster::extent(dem) <- raster::extent(rasterLayer)
    # specify crs
    raster::crs(dem) <- raster::crs(rasterLayer)
  }
  dem
}
#
################################
#' RasterLayer to Cimg conversion
#' 
#' converts a RasterLayer object to Cimg object. NA values in raster are replaced.
#'
#' @param rasterLayer raster object. raster of canopy height model, preferably filtered to avoid effect of holes on volume and surface computation
#' @param NA.replace numeric. value to replace NA values with.
#' @param maxpixels numeric. maximum number of pixels to be converted to cimg (argument passed to as.cimg).
#' @return A cimg object
#' @examples
#' data(chmchablais3)
#' 
#' chm.cim <- raster2Cimg(chmchablais3)
#' chm.cim
#' summary(chm.cim)
#' 
#' # plot rasterLayer
#' raster::plot(chmchablais3)
#' 
#' # plot cimg object
#' plot(chm.cim)
#' @seealso \code{\link{cimg2Raster}}
#' @export
raster2Cimg <- function(rasterLayer, NA.replace=0, maxpixels=1e+10)
{
  # replace NA values
  rasterLayer[is.na(rasterLayer)] <- NA.replace
  # convert to cimg object
  if (rasterLayer@ncols*rasterLayer@nrows>maxpixels) {warning("Too many rasterLayer pixels: conversion to cimg partial; try with higher maxpixels arguments")}
  imager::as.cimg(rasterLayer, maxpixels=maxpixels)
}