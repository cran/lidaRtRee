# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
#-------------------------------------------------------------------------------
#' Gap detection in a Canopy Height Model
#'
#' Performs gaps detection in a canopy height model. Function
#' \code{\link{dem_filtering}} is first applied to the canopy height model to
#' remove artefacts. Gaps are then extracted based on several criteria:
#' \enumerate{
#' \item Vegetation height must be smaller than a threshold
#' \item Gap width must be large enough, depending on surrounding canopy height;
#' distance to surrounding vegetation is tested with morphological closings
#' \item Gap must have a minimum surface
#' }
#'
#'
#' @param chm raster object. canopy height model
#' @param ratio numeric. maximum ratio between surrounding canopy height and gap
#' distance (a pixel belongs to the gap only if for any vegetation pixel around
#' it, the distance to the vegetation pixel is larger than pixel height/ratio).
#' If ratio is set to NULL, this criterion is not taken into account
#' @param gap_max_height numeric. maximum canopy height to be considered as gap
#' @param min_gap_surface numeric. minimum gap surface
#' @param max_gap_surface numeric. maximum gap surface
#' @param closing_height_bin numeric. height bin width for morphological closing
#' of gaps to test ratio between canopy height and gap distance
#' @param nl_filter string. type of non-linear filter to apply to canopy height
#' model to remove artefacts, should be an option of \code{\link{dem_filtering}}
#' @param nl_size numeric. kernel width in pixel for non-linear filtering
#' @param gap_reconstruct boolean. default behaviour is that areas that do not
#' fulfill the ratio criterion are removed from gaps. If set to TRUE, in case
#' some pixels of a gap fulfill the distance criterion, the connected pixels that
#' fulfill the height criterion are also integrated to it.
#' @examples
#' data(chm_chablais3)
#'
#' # fill NA values in canopy height model
#' chm_chablais3[is.na(chm_chablais3)] <- 0
#'
#' # gap detection with distance larger than canopy height / 2
#' gaps <- gap_detection(chm_chablais3, ratio = 2, gap_max_height = 1, 
#' min_gap_surface = 0)
#'
#' # gap detection with distance larger than canopy height / 2
#' # and reconstruction of border areas
#' gaps1 <- gap_detection(chm_chablais3,
#'   ratio = 2, gap_max_height = 1, min_gap_surface = 0,
#'   gap_reconstruct = TRUE
#' )
#'
#' # gap detection without distance criterion
#' gaps2 <- gap_detection(chm_chablais3, ratio = NULL, gap_max_height = 1, 
#' min_gap_surface = 0)
#'
#' # gap id and corresponding surface for third detection parameters
#' table(raster::values(gaps2$gap_id)) * raster::res(gaps2$gap_id)[1]^2
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Initial image")
#'
#' # plot binary image of gaps
#' raster::plot(gaps$gap_id > 0, main = "Gaps", legend = FALSE)
#' raster::plot(gaps1$gap_id > 0, main = "Gaps, with reconstruction", legend = FALSE)
#' raster::plot(gaps2$gap_id > 0, main = "Gaps, no width criterion", legend = FALSE)
#'
#' # plot filtered CHM
#' raster::plot(gaps2$filled_chm, main = "Filtered CHM")
#' @seealso \code{\link{dem_filtering}}, \code{\link{edge_detection}}
#' @return A list of three raster objects: raster with gap labels, raster with gap surface, canopy height model after filter.
#' @export
gap_detection <- function(chm, ratio = 2, gap_max_height = 1, min_gap_surface = 25, 
                          max_gap_surface = +Inf, closing_height_bin = 1, 
                          nl_filter = "Median", nl_size = 3, gap_reconstruct = FALSE) {
  # convert to cimg object
  c_chm <- raster2Cimg(chm)
  # apply non linear filter to chm
  c_chm <- dem_filtering(c_chm, nl_filter, nl_size)[[1]]
  # convert to raster for export
  r.nl <- cimg2Raster(c_chm, chm)
  #
  # check gap width
  if (is.null(ratio)) {
    l_im <- list(imager::as.cimg(c_chm > gap_max_height))
  } else {
    l_im <- list()
    # loop on dilate size -> canopy height (threshold at 50 m in case of outliers)
    for (i in seq(from = gap_max_height, 
                  to = max(gap_max_height, min(max(c_chm), 50)), 
                  by = closing_height_bin))
    {
      # create binary image to close (areas where chm> i)
      dummy <- imager::as.cimg(c_chm > i)
      # create stucturing element (uneven disk of radius i/2)
      strel <- imager::as.cimg(create_disk(floor(i / ratio / raster::xres(chm) / 2) * 2 + 1))
      # perform morphological closing and store in list
      l_im[[as.character(i)]] <- imager::mclosing(dummy, strel)
    }
  }
  #
  # union of closed images -> non gap areas
  final <- imager::parmax(l_im)
  # compute gap areas not closed
  gaps <- abs(final - 1)
  if (gap_reconstruct) {
    # extend non closed gaps into connected pixels where h < gap_max_height
    # map of pixels which comply with the height criterion (gap candidates)
    gaps_candidate <- c_chm < gap_max_height
    # label all gap candidates
    labels <- (imager::label(gaps_candidate) + 1) * gaps_candidate
    # list of labels of gaps not closed when applying the distance ratio
    not_closed_labels <- setdiff(unique(labels * gaps), 0)
    # remove closed labels
    gaps <- cimg2Raster(labels, chm)
    raster::values(gaps)[!is.element(raster::values(gaps), not_closed_labels)] <- 0
    gaps <- raster2Cimg(gaps > 0)
  }
  #
  # label unconnected gaps
  labels <- (imager::label(gaps) + 1) * gaps
  # extract gap surface
  gap_surface <- table(as.vector(labels)) * (raster::xres(chm))^2
  gap_surface <- as.data.frame(gap_surface)
  # remove non-gap category (0), except if only one gap
  if (nrow(gap_surface) > 1) {
    gap_surface <- gap_surface[-1, ]
  }
  #
  # convert gap label to raster object
  r_labels <- cimg2Raster(labels, chm)
  # set label of non gaps to NA
  r_labels[r_labels == 0] <- NA
  # create map where pixel value is gap size
  r_surface <- r_labels
  raster::values(r_surface) <- NA
  # PROBABLY A BETTER WAY TO DO IT THAN IN A LOOP
  for (i in 1:nrow(gap_surface))
  {
    r_surface[r_labels == as.numeric(as.character(gap_surface$Var1[i]))] <- gap_surface$Freq[i]
  }
  # set surface of non gaps to NA
  r_surface[is.na(r_labels)] <- NA
  # removal of labels with small surface
  dummy <- (r_surface < min_gap_surface) | (r_surface > max_gap_surface)
  r_labels[dummy] <- r_surface[dummy] <- NA
  list(gap_id = r_labels, gap_surface = r_surface, filled_chm = r.nl)
}

#-------------------------------------------------------------------------------
#' Edge detection in gap image
#'
#' Performs edge detection on a gap image (e.g. output from function 
#' \code{\link{gap_detection}}). The gap image is compared to a gap image which 
#' has undergone a dilation or erosion to identify edges of gaps.
#'
#' @param gaps raster object. gaps image where 1 represents gaps and 0 non-gaps 
#' areas
#' @param inside boolean. defines where the edge is extracted: either inside the 
#' gaps (an erosion is applied to the gaps image) or outside (a dilation is applied)
#' @examples
#' data(chm_chablais3)
#'
#' # fill NA values in canopy height model
#' chm_chablais3[is.na(chm_chablais3)] <- 0
#'
#' # gap detection with distance larger than canopy height / 2
#' gaps <- gap_detection(chm_chablais3,
#'   ratio = 2, gap_max_height = 1, min_gap_surface = 10,
#'   gap_reconstruct = TRUE
#' )
#'
#' # edge detection
#' edges_inside <- edge_detection(!is.na(gaps$gap_id))
#' edges_outside <- edge_detection(!is.na(gaps$gap_id), inside = FALSE)
#'
#' # edge propotion
#' sum(raster::values(edges_inside)) / (nrow(edges_inside) * ncol(edges_inside))
#' sum(raster::values(edges_outside)) / (nrow(edges_outside) * ncol(edges_outside))
#'
#' # plot original image
#' raster::plot(chm_chablais3, main = "Initial image")
#'
#' # plot binary image of gaps
#' raster::plot(gaps$gap_id > 0, main = "Gaps", legend = FALSE)
#'
#' # plot edges
#' raster::plot(edges_inside, main = "Edges (inside)", legend = FALSE)
#' raster::plot(edges_outside, main = "Edges (outside)", legend = FALSE)
#' @seealso \code{\link{gap_detection}}
#' @return A raster object where edges are labelled as 1.
#' @export
edge_detection <- function(gaps, inside = TRUE) {
  # convert to cimg object
  c_gap <- raster2Cimg(gaps)
  # create structuring element before morphological operation
  mask <- create_disk(width = 3)
  # convert to cimg object
  c_mask <- imager::as.cimg(mask)
  # apply morphological operation
  if (inside) {
    c_morpho <- imager::erode(c_gap, c_mask)
  } else {
    c_morpho <- imager::dilate(c_gap, c_mask)
  }
  # output edges
  c_edges <- c_morpho != c_gap
  # convert to raster
  cimg2Raster(c_edges, gaps)
}
