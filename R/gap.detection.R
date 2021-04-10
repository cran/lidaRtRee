# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
########################### FUNCTIONS FOR GAP DETECTION #############################
#
######################################
#' Gap detection in a Canopy Height Model
#' 
#' Performs gaps detection in a canopy height model. Function \code{\link{demFiltering}} is first applied to the canopy height model to remove artefacts. Gaps are then extracted based on several criteria:
#' \enumerate{
#' \item Vegetation height must be smaller than a threshold
#' \item Gap width must be large enough, depending on surrounding canopy height; distance to surrounding vegetation is tested with morphological closings
#' \item Gap must have a minimum surface
#' }
#' 
#'
#' @param chm raster object. canopy height model
#' @param ratio numeric. maximum ratio between surrounding canopy height and gap distance (a pixel belongs to the gap only if for any vegetation pixel around it, the distance to the vegetation pixel is larger than pixel height/ratio). If ratio is set to NULL, this criterion is not taken into account
#' @param gap.max.height numeric. maximum canopy height to be considered as gap
#' @param min.gap.surface numeric. minimum gap surface
#' @param max.gap.surface numeric. maximum gap surface
#' @param closing.height.bin numeric. height bin width for morphological closing of gaps to test ratio between canopy height and gap distance
#' @param nlFilter string. type of non-linear filter to apply to canopy height model to remove artefacts, should be an option of \code{\link{demFiltering}}
#' @param nlsize numeric. kernel width in pixel for non-linear filtering
#' @param gapReconstruct boolean. default behaviour is that areas that do not fulfill the ratio criterion are removed from gaps. If set to TRUE, in case some pixels of a gap fulfill the distance criterion, the connected pixels that fulfill the height criterion are also integrated to it. 
#' @examples 
#' data(chmchablais3)
#' 
#' # fill NA values in canopy height model
#' chmchablais3[is.na(chmchablais3)] <- 0
#' 
#' # gap detection with distance larger than canopy height / 2
#' gaps <- gapDetection(chmchablais3, ratio=2, gap.max.height=1, min.gap.surface=0)
#' 
#' # gap detection with distance larger than canopy height / 2
#' # and reconstruction of border areas
#' gaps1 <- gapDetection(chmchablais3, ratio=2, gap.max.height=1, min.gap.surface=0,
#' gapReconstruct=TRUE)
#' 
#' # gap detection without distance criterion
#' gaps2 <- gapDetection(chmchablais3, ratio=NULL, gap.max.height=1, min.gap.surface=0)
#' 
#' # gap id and corresponding surface for third detection parameters
#' table(raster::values(gaps2$gap.id))*raster::res(gaps2$gap.id)[1]^2
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Initial image")
#' 
#' # plot binary image of gaps
#' raster::plot(gaps$gap.id>0, main="Gaps", legend=FALSE)
#' raster::plot(gaps1$gap.id>0, main="Gaps, with reconstruction", legend=FALSE)
#' raster::plot(gaps2$gap.id>0, main="Gaps, no width criterion", legend=FALSE)
#' 
#' # plot filtered CHM
#' raster::plot(gaps2$filled.chm, main="Filtered CHM")
#' 
#' @seealso \code{\link{demFiltering}}, \code{\link{edgeDetection}}
#' @return A list of three raster objects: raster with gap labels, raster with gap surface, canopy height model after filter.
#' @export
gapDetection <- function(chm, ratio=2, gap.max.height=1, min.gap.surface=25, max.gap.surface=+Inf, closing.height.bin=1, nlFilter="Median", nlsize=3, gapReconstruct=FALSE)
{
  # convert to cimg object
  c.chm <- raster2Cimg(chm)
  # apply non linear filter to chm
  c.chm <- demFiltering(c.chm,nlFilter,nlsize)[[1]]
  # convert to raster for export
  r.nl <- cimg2Raster(c.chm, chm)
  #
  # check gap width
  if (is.null(ratio))
  {
    l.im <- list(imager::as.cimg(c.chm > gap.max.height))
  } else {
    l.im <- list()
    # loop on dilate size -> canopy height (threshold at 50 m in case of outliers)    
    for(i in seq(from=gap.max.height, to=max(gap.max.height,min(max(c.chm),50)), by=closing.height.bin))
    {
      # create binary image to close (areas where chm> i)
      dummy <- imager::as.cimg(c.chm > i)
      # create stucturing element (uneven disk of radius i/2)
      strel <- imager::as.cimg(createDisk(floor(i/ratio/raster::xres(chm)/2)*2+1))
      # perform morphological closing and store in list
      l.im[[as.character(i)]] <- imager::mclosing(dummy,strel)
    }
  }
  #
  # union of closed images -> non gap areas
  final <- imager::parmax(l.im)
  # compute gap areas not closed
  gaps <- abs(final-1)
  if (gapReconstruct)
  {
    # extend non closed gaps into connected pixels where h < gap.max.height
    # map of pixels which comply with the height criterion (gap candidates)
    gaps.candidate <- c.chm < gap.max.height
    # label all gap candidates
    labels <- (imager::label(gaps.candidate) +1 ) * gaps.candidate
    # list of labels of gaps not closed when applying the distance ratio
    not.closed.labels <- setdiff(unique(labels * gaps), 0)
    # remove closed labels
    gaps <- cimg2Raster(labels, chm)
    raster::values(gaps)[!is.element(raster::values(gaps), not.closed.labels)] <- 0
    gaps <- raster2Cimg(gaps>0)
  }
  #
  # label unconnected gaps
  labels <- (imager::label(gaps) +1 ) * gaps
  # extract gap surface
  gap.surface <- table(as.vector(labels)) * (raster::xres(chm))^2
  gap.surface <- as.data.frame(gap.surface)
  # remove non-gap category (0), except if only one gap
  if (nrow(gap.surface)>1)
  {
    gap.surface <- gap.surface[-1,]
  }
  #
  # convert gap label to raster object
  r.labels <- cimg2Raster(labels, chm)
  # set label of non gaps to NA
  r.labels[r.labels==0] <- NA
  # create map where pixel value is gap size
  r.surface <- r.labels
  raster::values(r.surface) <- NA
  # PROBABLY A BETTER WAY TO DO IT THAN IN A LOOP
  for (i in 1:nrow(gap.surface))
  {
    r.surface[r.labels == as.numeric(as.character(gap.surface$Var1[i]))] <- gap.surface$Freq[i]
  }
  # set surface of non gaps to NA
  r.surface[is.na(r.labels)] <- NA
  # removal of labels with small surface
  dummy <- (r.surface < min.gap.surface) | (r.surface > max.gap.surface)
  r.labels[dummy] <-r.surface[dummy] <- NA
  list(gap.id=r.labels, gap.surface=r.surface, filled.chm=r.nl) 
}
######################################
#' Edge detection in gap image
#' 
#' Performs edge detection on a gap image (e.g. output from function \code{\link{gapDetection}}). The gap image is compared to a gap image which has undergone a dilation or erosion to identify edges of gaps.
#'
#' @param gaps raster object. gaps image where 1 represents gaps and 0 non-gaps areas
#' @param inside boolean. defines where the edge is extracted: either inside the gaps (an erosion is applied to the gaps image) or outside (a dilation is applied)
#' @examples 
#' data(chmchablais3)
#' 
#' # fill NA values in canopy height model
#' chmchablais3[is.na(chmchablais3)] <- 0
#' 
#' # gap detection with distance larger than canopy height / 2
#' gaps <- gapDetection(chmchablais3, ratio=2, gap.max.height=1, min.gap.surface=10,
#' gapReconstruct=TRUE)
#' 
#' # edge detection
#' edges.inside <- edgeDetection(!is.na(gaps$gap.id))
#' edges.outside <- edgeDetection(!is.na(gaps$gap.id), inside=FALSE)
#' 
#' # edge propotion
#' sum(raster::values(edges.inside))/(nrow(edges.inside)*ncol(edges.inside))
#' sum(raster::values(edges.outside))/(nrow(edges.outside)*ncol(edges.outside))
#' 
#' # plot original image
#' raster::plot(chmchablais3, main="Initial image")
#' 
#' # plot binary image of gaps
#' raster::plot(gaps$gap.id>0, main="Gaps", legend=FALSE)
#' 
#' # plot edges
#' raster::plot(edges.inside, main="Edges (inside)", legend=FALSE)
#' raster::plot(edges.outside, main="Edges (outside)", legend=FALSE)
#' 
#' @seealso \code{\link{gapDetection}}
#' @return A raster object where edges are labelled as 1.
#' @export
edgeDetection <- function(gaps, inside=TRUE)
{
  # convert to cimg object
  c.gap <- raster2Cimg(gaps)
  # create structuring element before morphological operation
  mask <- createDisk(width=3)
  # convert to cimg object
  c.mask <- imager::as.cimg(mask)
  # apply morphological operation
  if (inside)
  {
    c.morpho <- imager::erode(c.gap, c.mask)
  } else {
    c.morpho <- imager::dilate(c.gap, c.mask)
  }
  # output edges
  c.edges <- c.morpho != c.gap
  # convert to raster
  cimg2Raster(c.edges, gaps)
}