# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
########################### FUNCTIONS FOR TREE DETECTION #######################
#
#-------------------------------------------------------------------------------
#' Disk-shaped matrix mask
#'
#' Creates a matrix with TRUE values shaping a centered disk
#'
#' @param width numeric. disk width in pixels, should be an uneven number
#' @return A matrix with 1 for pixels inside the disk, 0 outside
#' @examples
#' create_disk(7)
#' @export
create_disk <- function(width = 5) {
  if (width %% 2 != 1) {
    stop("Mask width should be uneven")
  }
  # matrix of row indices
  row.mat <- matrix(rep(1:width, width), nrow = width, ncol = width)
  # matrix of column indices
  col.mat <- matrix(rep(1:width, width), nrow = width, ncol = width, byrow = TRUE)
  radius <- width %/% 2
  row.mat <- row.mat - radius - 1
  col.mat <- col.mat - radius - 1
  mask <- (col.mat^2 + row.mat^2) <= radius^2
  mask
}

#-------------------------------------------------------------------------------
#' Image pre-processing (non-linear filtering and Gaussian smoothing)
#'
#' applies two filters to an image:
#' \enumerate{
#' \item A non-linear filter: closing (\code{\link[imager]{mclosing}}) with disk 
#' kernel, or median (\code{\link[imager]{medianblur}}) with square kernel
#' \item A 2D Gaussian smoother (The \code{\link[imager]{deriche}} filter is 
#' applied on both dimensions). Value-dependent smoothing is possible
#' }
#'
#' @param dem cimg object (e.g. obtained with \code{\link[imager]{as.cimg}}) or 
#' Raster Layer object (e.g. obtained with \code{\link[raster]{raster}})
#' @param nl_filter string. type of non-linear filter to apply: "None", "Closing" 
#' or "Median"
#' @param nl_size numeric. kernel width in pixel for non-linear filtering
#' @param sigmap numeric or matrix. if a single number is provided, sigmap is 
#' the standard deviation of the Gaussian filter in pixel, 0 corresponds to no 
#' smoothing. In case of matrix, the first column corresponds to the standard 
#' deviation of the filter, and the second to thresholds for image values (e.g. 
#' a filter of standard deviation specified in line \code{i} is applied to pixels 
#' in image which values are between thresholds indicated in lines \code{i} and 
#' \code{i+1}). Threshold values should be ordered in increasing order.
#' @param padding boolean. Whether image should be padded by duplicating edge 
#' values before filtering to avoid border effects
#' @examples
#' data(chm_chablais3)
#'
#' # filtering with median and Gaussian smoothing
#' im <- dem_filtering(chm_chablais3, nl_filter = "Median", nl_size = 3, sigmap = 0.8)
#'
#' # filtering with median filter and value-dependent Gaussian smoothing
#' # (less smoothing for values between 0 and 15)
#' im2 <- dem_filtering(chm_chablais3,
#'   nl_filter = "Median", nl_size = 3,
#'   sigmap = cbind(c(0.2, 0.8), c(0, 15))
#' )
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Initial image")
#'
#' # plot image after median filter
#' raster::plot(im$non_linear_image, main = "Median filter")
#'
#' # plot image after median and Gaussian filters
#' raster::plot(im$smoothed_image, main = "Smoothed image")
#'
#' # plot image after median and value-dependent Gaussian filters
#' raster::plot(im2$smoothed_image, main = "Value-dependent smoothing")
#' @seealso \code{\link{maxima_detection}}, filters of imager package: 
#' \code{\link[imager]{mclosing}}, \code{\link[imager]{medianblur}}, 
#' \code{\link[imager]{deriche}}
#' @return A list of two cimg or a RasterStack objects: image after non-linear 
#' filter and image after both filters
#' @export
dem_filtering <- function(dem, nl_filter = "Closing", nl_size = 5, sigmap = 0.3, 
                          padding = TRUE) {
  # convert rasterLayer to cimg object if necessary
  if (class(dem)[1] == "RasterLayer") {
    dem.c <- raster2Cimg(dem)
  } else {
    dem.c <- dem
  }
  #
  if (padding) {
    # padding number of cells is maximum of half width of non linear filter or 
    # ceiling value of three times sigmap
    if (!is.null(dim(sigmap))) {
      dummy <- max(sigmap[, 1])
    } else {
      dummy <- sigmap
    }
    border.size <- max((nl_size - 1) / 2 + 1, ceiling(dummy * 3))
    # convert to matrix
    dem.c <- as.matrix(dem.c)
    # duplicate columns
    dem.c <- dem.c[, c(rep.int(1, border.size), 1:ncol(dem.c), 
                       rep.int(ncol(dem.c), border.size))]
    # duplicate rows
    dem.c <- dem.c[c(rep.int(1, border.size), 1:nrow(dem.c), 
                     rep.int(nrow(dem.c), border.size)), ]
    dem.c <- imager::as.cimg(dem.c)
  }
  #
  # non-linear filtering
  dem_nl <- dem.c
  if (nl_filter == "Closing") {
    # closing with structuring element (disk)
    strel <- imager::as.cimg(create_disk(nl_size))
    dem_nl <- imager::mclosing(dem_nl, strel)
  }
  if (nl_filter == "Median") {
    # median on square window
    # strel <- imager::as.cimg(create_disk(nl_size))
    dem_nl <- imager::medianblur(dem_nl, nl_size)
  }
  #
  # linear filtering
  # gaussian smoothing
  # if several values of sigma are provided, value-dependent smoothin is performed
  if (length(sigmap) > 1) {
    dem_gs <- dem_nl
    # for each value of sigma
    for (i in 1:nrow(sigmap))
    {
      # perform 2D smoothing
      dummy <- imager::deriche(dem_nl, sigmap[i, 1], axis = "x")
      dummy <- imager::deriche(dummy, sigmap[i, 1], axis = "y")
      # identify pixels which values are superior to the threshold
      temp <- dem_gs >= sigmap[i, 2]
      # update only those values in the output
      dem_gs[temp] <- dummy[temp]
    }
  } else { # if only one value of sigma is provided
    if (sigmap > 0) {
      dem_gs <- imager::deriche(dem_nl, sigmap, axis = "x")
      dem_gs <- imager::deriche(dem_gs, sigmap, axis = "y")
    } else { # if 0, no smoothing is performed
      dem_gs <- dem_nl
    }
  }
  # remove padding
  if (padding) {
    dem_nl <- as.matrix(dem_nl)
    dem_nl <- imager::as.cimg(dem_nl[(border.size + 1):(nrow(dem_nl) - border.size), 
                                     (border.size + 1):(ncol(dem_nl) - border.size)])
    dem_gs <- as.matrix(dem_gs)
    dem_gs <- imager::as.cimg(dem_gs[(border.size + 1):(nrow(dem_gs) - border.size), 
                                     (border.size + 1):(ncol(dem_gs) - border.size)])
  }
  # convert cimg objects to rasterLayer if necessary
  if (class(dem)[1] == "RasterLayer") {
    output <- raster::addLayer(cimg2Raster(dem_nl, dem), cimg2Raster(dem_gs, dem))
    names(output) <- c("non_linear_image", "smoothed_image")
  } else {
    output <- list(non_linear_image = dem_nl, smoothed_image = dem_gs)
  }
  output
}

#-------------------------------------------------------------------------------
#' Local maxima extraction on image
#'
#' Variable window size maxima detection is performed on the image to extract 
#' local maxima position and calculate the window size where they are global 
#' maxima. Gaussian white noise is added to the image to avoid adjacent maxima 
#' due to neighbor pixels with identical value.
#'
#' @param dem cimg object (e.g. as created by \code{\link[imager]{cimg}}) or 
#' RasterLayer object (e.g. obtained with \code{\link[raster]{raster}})
#' @param dem.res numeric. image resolution, in case \code{dem} is a rasterLayer 
#' object, \code{dem.res} is extracted from the object by \code{\link[raster]{res}}
#' @param max.width numeric. maximum kernel width in pixel to check for local 
#' maximum
#' @param jitter boolean. indicates if noise should be added to image values to 
#' avoid the adjacent maxima due to the adjacent pixels with equal values
#' @return A cimg object or RasterLayer object which values are the radius (n) 
#' in meter of the square window (width 2n+1) where the center pixel is global 
#' maximum
#' @examples
#' data(chm_chablais3)
#'
#' # maxima detection
#' maxi <- maxima_detection(chm_chablais3)
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Initial image")
#'
#' # plot maxima image
#' raster::plot(maxi, main = "Local maxima")
#' @seealso \code{\link{dem_filtering}}, \code{\link{maxima_selection}}
#' @export
maxima_detection <- function(dem, dem.res = 1, max.width = 21, jitter = TRUE) {
  # convert rasterLayer to cimg object if necessary
  if (class(dem)[1] == "RasterLayer") {
    dem_gs <- raster2Cimg(dem)
    dem.res <- raster::res(dem)[1]
  } else {
    dem_gs <- dem
  }
  #
  # add absolute of gaussian white noise, mean=0, sd=sd(dem_gs)/100000 to non 0 pixels
  if (jitter) {
    dem_gs <- dem_gs + abs(imager::imnoise(
      dim = dim(dem_gs), 
      sd = stats::sd(dem_gs) / 100000)) * (dem_gs != 0)
  }
  #
  max.radius <- max.width %/% 2
  # extraction of maxima on variable window size (from 0 to max.radius)
  maxi <- NULL
  for (i in 1:max.radius)
  {
    # create square structuring element for dilation
    strel <- imager::imfill(2 * i + 1, 2 * i + 1, val = 1)
    # if initialization check where pixel values are equal in original image and 
    # in image dilated by structuring element
    if (is.null(maxi)) {
      maxi <- imager::as.cimg(dem_gs == imager::dilate(dem_gs, strel)) * i
      # otherwise perform check and update previous result
    } else {
      maxi <- imager::parmax(list(maxi, 
                                  imager::as.cimg(dem_gs == imager::dilate(dem_gs, strel)) * i))
    }
  }
  # convert window size from pixels to meters
  maxi[maxi > 0] <- (maxi[maxi > 0] + 1) * dem.res
  # convert cimg objects to rasterLayer if necessary
  if (class(dem)[1] == "RasterLayer") {
    maxi <- cimg2Raster(maxi, dem)
  }
  maxi
}

#-------------------------------------------------------------------------------
#' Image maxima selection based on values and neighborhood of local maxima
#'
#' In a maxima image (output of \code{\link{maxima_detection}}), sets values to 
#' zero for pixels which
#' \enumerate{
#' \item value in the initial image (from which maxima were detected) are below 
#' a threshold
#' \item values in the maxima image (corresponding to the radius of the 
#' neighborhood where they are global maxima) are below a threshold depending on 
#' the initial image value.
#' }
#'
#' @param maxi cimg object or RasterLayer object. image with local maxima 
#' (typically output from \code{\link{maxima_detection}}, image values correspond 
#' to neighborhood radius on which pixels are global maxima in the initial image)
#' @param dem_nl cimg object. initial image from which maxima were detected
#' @param hmin numeric. minimum value in initial image for a maximum to be selected
#' @param dmin numeric. intercept term for selection of maxima depending on 
#' neighborhood radius: \code{maxi >= dmin + dem_nl * dprop}
#' @param dprop numeric. proportional term for selection of maxima depending on 
#' neighborhood radius: \code{maxi >= dmin + dem_nl * dprop}
#' @return A cimg object or rasterLayer object which values are the radius (n) 
#' in meter of the square window (width 2n+1) where the center pixel is global 
#' maximum and which fulfill the selection criteria
#' @examples
#' data(chm_chablais3)
#'
#' # maxima detection
#' maxi <- maxima_detection(chm_chablais3)
#'
#' # several maxima selection settings
#' selected_maxi_hmin <- maxima_selection(maxi, chm_chablais3, hmin = 15)
#' selected_maxi_dm <- maxima_selection(maxi, chm_chablais3, dm = 2.5)
#' selected_maxi <- maxima_selection(maxi, chm_chablais3, dm = 1, dprop = 0.1)
#'
#' # corresponding count number of remaining maxima
#' table(raster::values(maxi))
#' table(raster::values(selected_maxi_hmin))
#' table(raster::values(selected_maxi_dm))
#' table(raster::values(selected_maxi))
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Initial image")
#'
#' # plot maxima images, original and first case
#' raster::plot(maxi, main = "Local maxima")
#' raster::plot(selected_maxi, main = "Selected maxima")
#' @seealso \code{\link{maxima_detection}}
#' @export
maxima_selection <- function(maxi, dem_nl, hmin = 0, dmin = 0, dprop = 0) {
  # convert rasterLayer to cimg object if necessary
  if (class(maxi)[1] == "RasterLayer") {
    israster <- TRUE
    dem <- maxi
    maxi <- raster2Cimg(maxi)
    dem_nl <- raster2Cimg(dem_nl)
  } else {
    israster <- FALSE
  }
  # height filter
  maxi[dem_nl < hmin] <- 0
  # distance filter
  maxi[maxi < (dmin + dem_nl * dprop)] <- 0
  # convert to rasterLayer if necessary
  if (israster) {
    cimg2Raster(maxi, dem)
  } else {
    maxi
  }
}

#-------------------------------------------------------------------------------
#' Image segmentation by seed-based watershed algorithm
#'
#' performs a seed-based watershed segmentation (wrapper for imager::watershed)
#'
#' @param maxi cimg or rasterLayer object. image with seed points (e.g. from 
#' \code{\link{maxima_detection}} or \code{\link{maxima_selection}})
#' @param dem_nl cimg or rasterLayer object. image for seed propagation 
#' (typically initial image used for maxima detection).
#' @return A cimg object or rasterlayer object with segments id
#' @examples
#' data(chm_chablais3)
#'
#' # median filter
#' chm_chablais3 <- dem_filtering(chm_chablais3,
#'   nl_filter = "Median", nl_size = 3,
#'   sigmap = 0
#' )$non_linear_image
#'
#' # maxima detection
#' maxi <- maxima_detection(chm_chablais3)
#'
#' # maxima selection
#' selected_maxi <- maxima_selection(maxi, chm_chablais3, dm = 1, dprop = 0.1)
#'
#' # segmentation
#' seg_maxi <- segmentation(maxi, chm_chablais3)
#' seg_selected_maxi <- segmentation(selected_maxi, chm_chablais3)
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Median filter")
#'
#' # plot segmented image
#' # replace segment with id 0 (not a tree) with NA
#' seg_maxi[seg_maxi == 0] <- NA
#' raster::plot(seg_maxi %% 8, main = "Segments, no maxima selection", 
#' col = rainbow(8))
#' seg_selected_maxi [seg_selected_maxi == 0] <- NA
#' raster::plot(seg_selected_maxi %% 8, main = "Segments, maxima selection", 
#' col = rainbow(8))
#' @seealso \code{\link{maxima_detection}}, \code{\link{maxima_selection}}, 
#' \code{\link{seg_adjust}}
#' @export
segmentation <- function(maxi, dem_nl) {
  # convert rasterLayer to cimg object if necessary
  if (class(maxi)[1] == "RasterLayer") {
    israster <- TRUE
    dem <- maxi
    maxi <- raster2Cimg(maxi)
    dem_nl <- raster2Cimg(dem_nl)
  } else {
    israster <- FALSE
  }
  # maxima locations
  a <- which(maxi > 0)
  # tree number
  na <- length(a)
  # create seed image
  dummy <- maxi
  # initialise seeds with id
  dummy[a] <- sample(1:na, na)
  # watershed segmentation
  dem_w <- imager::watershed(dummy, dem_nl)
  #
  # convert to raster if necessary
  if (israster) {
    cimg2Raster(dem_w, dem)
  } else {
    dem_w
  }
}

#-------------------------------------------------------------------------------
#' Image statistic in segment
#'
#' compute zonal statistic of an image
#'
#' @param segms cimg or rasterLayer object. image with segments id (e.g. from 
#' \code{\link{segmentation}})
#' @param dem_nl cimg or rasterLayer object. image to compute statistic from
#' @param fun function to compute statistis from values in each segment
#' @return A cimg object or raster object with values of the statistic
#' @examples
#' data(chm_chablais3)
#'
#' # median filter
#' chm_chablais3 <- dem_filtering(chm_chablais3,
#'   nl_filter = "Median", nl_size = 3,
#'   sigmap = 0
#' )$non_linear_image
#'
#' # maxima detection
#' maxi <- maxima_detection(chm_chablais3)
#'
#' # segmentation
#' seg_maxi <- segmentation(maxi, chm_chablais3)
#'
#' # compute image of maximum value in each segment
#' max_in_segment <- raster_zonal_stats(seg_maxi, chm_chablais3)
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Median filter")
#'
#' # plot segments and image of max value inside segments
#' seg_maxi[seg_maxi == 0] <- NA
#' raster::plot(seg_maxi %% 8, main = "Segments", col = rainbow(8))
#' raster::plot(max_in_segment, main = "Max value in segment")
#' @seealso \code{\link{segmentation}}
#' @export
raster_zonal_stats <- function(segms, dem_nl, fun = max) {
  # convert rasterLayer to cimg object if necessary
  if (class(segms)[1] == "RasterLayer") {
    israster <- TRUE
    dem <- segms
    segms <- raster2Cimg(segms)
    dem_nl <- raster2Cimg(dem_nl)
  } else {
    israster <- FALSE
  }
  dummy <- stats::aggregate(as.numeric(dem_nl), by = list(as.numeric(segms)), 
                            FUN = fun)
  dem_wh <- dummy$x[match(as.numeric(segms), dummy$Group.1)]
  #
  dem_wh <- imager::as.cimg(matrix(dem_wh, nrow(segms)))
  #
  # convert to raster if necessary
  if (israster) {
    cimg2Raster(dem_wh, dem)
  } else {
    dem_wh
  }
}

#-------------------------------------------------------------------------------
#' Modification of segments based on values
#'
#' in a segmented image, removes from segments the pixels which values in a 
#' reference image is below a certain percentage of the highest value inside the 
#' segment. Removed pixels are attributed 0 value.
#'
#' @param dem_w cimg or rasterLayer object. image with segments id, without 0 
#' values
#' @param dem_wh cimg or rasterLayer object. image with max value inside segment
#' @param dem_nl cimg or rasterLayer object. image with initial values
#' @param prop numeric. proportional threshold for removal of pixels which initial 
#' values are lower than the max height of the segment (\code{dem_nl < prop x dem_wh})
#' @param min.value numeric. threshold for removel of pixels which initial values 
#' are lower (\code{dem_nl < min.value})
#' @param min.maxvalue numeric. threshold for complete removal of segments which 
#' maximum value height is smaller to the threshold (\code{dem_wh < min.maxvalue})
#' @return A cimg or rasterLayer object: image with modified segments.
#' @examples
#' data(chm_chablais3)
#'
#' # median filter
#' chm_chablais3 <- dem_filtering(chm_chablais3,
#'   nl_filter = "Median", nl_size = 3,
#'   sigmap = 0
#' )$non_linear_image
#'
#' # maxima detection and selection
#' maxi <- maxima_detection(chm_chablais3)
#' selected_maxi <- maxima_selection(maxi, chm_chablais3, dm = 1, dprop = 0.1)
#'
#' # segmentation
#' seg_selected_maxi <- segmentation(selected_maxi, chm_chablais3)
#'
#' # max value in segments
#' max_in_segment <- raster_zonal_stats(seg_selected_maxi , chm_chablais3)
#'
#' # segmentation modification
#' seg_modif1 <- seg_adjust(seg_selected_maxi , max_in_segment,
#'   chm_chablais3,
#'   prop = 0.5
#' )
#' seg_modif2 <- seg_adjust(seg_selected_maxi , max_in_segment,
#'   chm_chablais3,
#'   prop = 0, min.value = 5, min.maxvalue = 10
#' )
#'
#' # plot initial segmented image
#' seg_selected_maxi [seg_selected_maxi == 0] <- NA
#' raster::plot(seg_selected_maxi %% 8, main = "Initial segments", col = rainbow(8))
#' seg_modif1[seg_modif1 == 0] <- NA
#' raster::plot(seg_modif1 %% 8, main = "Modified segments 1", col = rainbow(8))
#' seg_modif2[seg_modif2 == 0] <- NA
#' raster::plot(seg_modif2 %% 8, main = "Modified segments 2", col = rainbow(8))
#' @seealso \code{\link{maxima_detection}}, \code{\link{maxima_selection}}
#' @export
seg_adjust <- function(dem_w, dem_wh, dem_nl, prop = 0.3, min.value = 2, min.maxvalue = 5) {
  # convert rasterLayer to cimg object if necessary
  if (class(dem_w)[1] == "RasterLayer") {
    israster <- TRUE
    dem <- dem_w
    dem_w <- raster2Cimg(dem_w)
    dem_wh <- raster2Cimg(dem_wh)
    dem_nl <- raster2Cimg(dem_nl)
  } else {
    israster <- FALSE
  }
  # removal of segments with maximum value lower than min.maxvalue
  dem_w[dem_wh < min.maxvalue] <- 0
  # removal of pixels with values lower than prop * max value in segment
  dem_w[dem_nl < prop * dem_wh] <- 0
  # removal of pixels with values lower than min.value
  dem_w[dem_nl < min.value] <- 0
  # convert to raster if necessary
  if (israster) {
    dem_w <- cimg2Raster(dem_w, dem)
  }
  dem_w
}

#-------------------------------------------------------------------------------
#' Preprocessing and segmentation of raster image for tree identification
#'
#' global function for preprocessing (filtering), maxima detection and selection, 
#' segmentation and segmentation adjustment of a raster image.
#'
#' @param dem raster object or string indicating location of raster file 
#' (typically a canopy height model or a digital surface model; in the latter 
#' case the dtm parameter should be provided)
#' @param nl_filter string. specifies the non-linear filter for image pre-processing, 
#' should be an option of function \code{\link{dem_filtering}}
#' @param nl_size numeric. width of kernel of non-linear filter in pixels
#' @param sigma numeric or matrix. if a single number is provided, sigmap is the 
#' standard deviation of Gaussian filter in meters, 0 corresponds to no smoothing. 
#' In case of matrix, the first column corresponds to the standard deviation of 
#' the filter, and the second to thresholds for image values (e.g. a filter of 
#' standard deviation specified in line \code{i} is applied to pixels in image 
#' which values are between thresholds indicated in lines \code{i} and 
#' \code{i+1}). Threshold values should be ordered in increasing order.
#' @param dmin numeric. treetop minimum distance to next higher pixel in meters
#' @param dprop numeric. number defining the treetop minimum distance as 
#' proportion of height to next higher pixel
#' @param hmin numeric. minimum treetop height
#' @param crown_prop numeric. minimum height of tree crown as proportion of 
#' treetop height
#' @param crown_hmin numeric. minimum crown height
#' @param dtm raster object or string indicating location of raster file with 
#' the terrain model. If provided, the maxima extraction and watershed segmentation 
#' are performed on the dem (this avoids the deformation of crown because of the 
#' normalisation with terrain), but maxima selection and segment adjustement are 
#' performed on 'dem-dtm' because the selection criteria is the height to terrain.
#' @references Monnet, J.-M. 2011. Using airborne laser scanning for mountain 
#' forests mapping: Support vector regression for stand parameters estimation 
#' and unsupervised training for treetop detection. Ph.D. thesis. University of 
#' Grenoble, France. Section 6.2 
#' \url{https://tel.archives-ouvertes.fr/tel-00652698/document}
#'
#' Monnet, J.-M., Mermin, E., Chanussot, J., Berger, F. 2010. Tree top detection 
#' using local maxima filtering: a parameter sensitivity analysis. Silvilaser 2010, 
#' the 10th International Conference on LiDAR Applications for Assessing Forest 
#' Ecosystems, September 14-17, Freiburg, Germany, 9 p. 
#' \url{https://hal.archives-ouvertes.fr/hal-00523245/document}
#' @examples
#' data(chm_chablais3)
#'
#' # tree segmentation
#' segments <- tree_segmentation(chm_chablais3)
#' segments2 <- tree_segmentation(chm_chablais3,
#'   nl_filter = "Median", nl_size = 3,
#'   sigma = cbind(c(0.2, 0.8), c(0, 15)), dmin = 0, dprop = 0, hmin = 10, 
#'   crown_prop = 0.5, crown_hmin = 5
#' )
#'
#' # plot initial image segments
#' raster::plot(chm_chablais3, main = "Initial image")
#' raster::plot(segments$smoothed_dem, main = "Filtered image")
#' raster::plot(segments$local_maxima, main = "Local maxima")
#' #
#' # replace segment with id 0 (not a tree) with NA
#' segments$segments_id[segments$segments_id == 0] <- NA
#' raster::plot(segments$segments_id %% 8, main = "Segments", col = rainbow(8))
#' #
#' # plot segmentation with other parameters
#' segments2$segments_id[segments2$segments_id == 0] <- NA
#' raster::plot(segments2$segments_id %% 8, main = "Segments2", col = rainbow(8))
#' @seealso \code{\link{dem_filtering}}, \code{\link{maxima_detection}}, 
#' \code{\link{maxima_detection}}, \code{\link{maxima_selection}}, 
#' \code{\link{segmentation}}, \code{\link{seg_adjust}}, \code{\link{tree_extraction}}
#' @return A RasterStack with 4 layers: selected local maxima (values = 
#' distance to higher pixel), segments, non-linear preprocessed dem, smoothed 
#' preprocessed dem
#' @export
tree_segmentation <- function(dem, nl_filter = "Closing", nl_size = 5, sigma = 0.3, 
                              dmin = 0, dprop = 0.05, hmin = 5, crown_prop = 0.3, 
                              crown_hmin = 2, dtm = NULL) {
  #
  if (crown_hmin > hmin) {
    stop("minimum tree height lower than minimum crown height")
  }
  if (is.character(dem)) {
    dem <- raster::raster(dem)
  }
  # convert NAs to 0s
  if (!is.null(dtm)) {
    if (is.character(dtm)) {
      dtm <- raster::raster(dtm)
    }
    dem[is.na(dem)] <- dtm[is.na(dem)]
  } else {
    dtm <- 0
    dem[is.na(dem)] <- dtm
  }
  #
  # conversion of Gaussian filter standard deviation from meters to pixels
  if (length(sigma == 1)) {
    sigma <- sigma / raster::res(dem)[1]
  } else {
    sigma[, 2] <- sigma[, 2] / raster::res(dem)[1]
  }
  # dem filtering
  temp <- dem_filtering(dem, nl_filter = nl_filter, nl_size = nl_size, sigmap = sigma)
  dem_nl <- temp$non_linear_image
  dem_gs <- temp$smoothed_image
  #
  # maxima detection
  maxi <- maxima_detection(dem_gs, raster::res(dem)[1])
  #
  # maxima selection
  maxi <- maxima_selection(maxi, dem_nl - dtm, 0, dmin, dprop) # no selection based on top height at this stage, otherwise some artefacts occur in the watershed step
  #
  # segmentation
  dem_w <- segmentation(maxi, dem_nl)
  dem_wh <- raster_zonal_stats(dem_w, dem_nl - dtm, fun = max)
  #
  # segmentation adjustment
  dem_w <- seg_adjust(dem_w, dem_wh, dem_nl - dtm, crown_prop, crown_hmin, hmin) # removal of trees with top heigth < hmin
  # remove trees from r_maxi
  maxi[dem_w == 0] <- 0
  #
  output <- raster::addLayer(maxi, dem_w, dem_nl, dem_gs)
  names(output) <- c("local_maxima", "segments_id", "filled_dem", "smoothed_dem")
  output
}

#-------------------------------------------------------------------------------
#' Tree extraction
#'
#' creates a dataframe with segment id, height and coordinates of maxima, surface and volume, computed from three images: initial, local maxima and segmented, obtained with \code{\link{tree_segmentation}}
#'
#' @param r_dem_nl raster object. raster of canopy height model, preferably filtered to avoid effect of holes on volume and surface computation
#' @param r_maxi raster object. raster positive values at local maxima
#' @param r_dem_w raster object. segmented raster
#' @param r_mask raster object. only segments which maxima are inside the mask are extracted
#' @return A spatial data.frame with tree id, local maximum stats (height, dominance radius), segment stats (surface and volume).
#' @examples
#' data(chm_chablais3)
#'
#' # tree segmentation
#' segments <- tree_segmentation(chm_chablais3)
#'
#' # tree extraction
#' trees <- tree_extraction(
#'   segments$filled_dem, segments$local_maxima,
#'   segments$segments_id
#' )
#' trees
#'
#' # plot initial image
#' raster::plot(chm_chablais3)
#'
#' # add treetop positions
#' sp::plot(trees, cex = trees$h / 20, add = TRUE, pch = 1)
#' \donttest{
#' # add segment contours (vectorization is slow)
#' contours <- raster::rasterToPolygons(segments$segments_id, dissolve = TRUE)
#' sp::plot(contours, add = TRUE, border = "white")
#' }
#'
#' @seealso \code{\link{tree_segmentation}}
#' @export
tree_extraction <- function(r_dem_nl, r_maxi, r_dem_w, r_mask = NULL) {
  # apply mask to remove maxima outside of region of interest
  if (!is.null(r_mask)) {
    r_maxi <- r_maxi * r_mask
  }
  # segments surface
  s <- data.frame(t(table(raster::values(r_dem_w))))[, -1]
  names(s) <- c("id", "s")
  s$s <- s$s * raster::res(r_dem_nl)[1] * raster::res(r_dem_nl)[2]
  # segments volume
  v <- stats::aggregate(raster::values(r_dem_nl), by = list(raster::values(r_dem_w)), FUN = sum)
  names(v) <- c("id", "v")
  v$v <- v$v * raster::res(r_dem_nl)[1] * raster::res(r_dem_nl)[2]
  #
  # compute surface and volume inside mask, if mask is not null
  if (!is.null(r_mask)) {
    # segments surface in mask
    sp <- data.frame(t(table(raster::values(r_dem_w * r_mask))))[, -1]
    names(sp) <- c("id", "sp")
    sp$sp <- sp$sp * raster::res(r_dem_nl)[1] * raster::res(r_dem_nl)[2]
    # segments volume in mask
    vp <- stats::aggregate(raster::values(r_dem_nl * r_mask), by = list(raster::values(r_dem_w)), FUN = sum)
    names(vp) <- c("id", "vp")
    vp$vp <- vp$vp * raster::res(r_dem_nl)[1] * raster::res(r_dem_nl)[2]
  }
  #
  # extract segment id, height and dom_radius
  cells <- which(raster::values(r_maxi) > 0)
  # check if positive values are present
  if (length(cells) == 0) {
    segms <- NULL
  } else {
    coord <- raster::xyFromCell(r_maxi, cells)
    segms <- data.frame(coord, id = raster::values(r_dem_w)[cells], 
                        h = raster::values(r_dem_nl)[cells], 
                        dom_radius = raster::values(r_maxi)[cells])
    segms <- merge(segms, s, all.x = TRUE)
    segms <- merge(segms, v, all.x = TRUE)
    if (!is.null(r_mask)) {
      segms <- merge(segms, sp, all.x = TRUE)
      segms <- merge(segms, vp, all.x = TRUE)
    }
    sp::coordinates(segms) <- ~ x + y
    segms@proj4string <- r_dem_nl@crs
  }
  segms
}

#-------------------------------------------------------------------------------
#' Cimg to RasterLayer conversion
#'
#' converts a cimg object to a RasterLayer object
#'
#' @param cimg raster object. raster of canopy height model, preferably filtered to avoid effect of holes on volume and surface computation
#' @param rasterLayer raster object. defines the extent and projection of conversion result
#' @return A RasterLayer
#' @examples
#' data(chm_chablais3)
#'
#' # convert rasterLayer to cimg object
#' chm_cim <- raster2Cimg(chm_chablais3)
#'
#' # apply filtering
#' chm_cim_filt <- dem_filtering(chm_cim,
#'   nl_filter = "Closing",
#'   nl_size = 3,
#'   sigmap = 0
#' )$non_linear_image
#'
#' # convert to RasterLayer
#' chm_filt <- cimg2Raster(chm_cim_filt, chm_chablais3)
#'
#' # plot rasterLayer
#' raster::plot(chm_chablais3)
#'
#' # plot cimg object
#' plot(chm_cim)
#'
#' # plot filtered cimg object
#' plot(chm_cim_filt)
#'
#' # plot filtered rasterLayer
#' raster::plot(chm_filt)
#' @seealso \code{\link{raster2Cimg}}
#' @export
cimg2Raster <- function(cimg, rasterLayer = NULL) {
  # convert to rasterLayer
  dem <- raster::raster(t(as.matrix(cimg)))
  # if reference is provided
  if (!is.null(rasterLayer)) {
    # specifiy extent
    raster::extent(dem) <- raster::extent(rasterLayer)
    # specify crs
    raster::crs(dem) <- raster::crs(rasterLayer)
  }
  dem
}

#-------------------------------------------------------------------------------
#' RasterLayer to Cimg conversion
#'
#' converts a RasterLayer object to Cimg object. NA values in raster are replaced.
#'
#' @param rasterLayer raster object. raster of canopy height model, preferably 
#' filtered to avoid effect of holes on volume and surface computation
#' @param NA_replace numeric. value to replace NA values with.
#' @param maxpixels numeric. maximum number of pixels to be converted to cimg 
#' (argument passed to as.cimg).
#' @return A cimg object
#' @examples
#' data(chm_chablais3)
#'
#' chm_cim <- raster2Cimg(chm_chablais3)
#' chm_cim
#' summary(chm_cim)
#'
#' # plot rasterLayer
#' raster::plot(chm_chablais3)
#'
#' # plot cimg object
#' plot(chm_cim)
#' @seealso \code{\link{cimg2Raster}}
#' @export
raster2Cimg <- function(rasterLayer, NA_replace = 0, maxpixels = 1e+10) {
  # replace NA values
  rasterLayer[is.na(rasterLayer)] <- NA_replace
  # convert to cimg object
  if (rasterLayer@ncols * rasterLayer@nrows > maxpixels) {
    warning("Too many rasterLayer pixels: conversion to cimg partial; try with higher maxpixels arguments")
  }
  imager::as.cimg(rasterLayer, maxpixels = maxpixels)
}
