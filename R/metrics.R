# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
################################
#' Computes metrics on list of point clouds
#' 
#' Computes metrics for a list of \code{\link[lidR]{LAS}} objects (should be normalized point clouds). Calls the function \code{\link[lidR]{cloud_metrics}} on each element and then arranges the results in a data.frame.
#'
#' @param llasn list of \code{\link[lidR]{LAS}} objects
#' @param func function. function applied on each element to compute metrics, default function is \code{\link[lidR]{stdmetrics}} from package \code{lidR}
#' @return A data frame with metrics in columns corresponding to LAS objects of the list (lines)
#' @seealso \code{\link[lidR]{cloud_metrics}}, \code{\link[lidR]{stdmetrics}}, \code{\link{ABAmodelMetrics}}
#' @examples
#' data(laschablais3)
#' 
#' # extract four point clouds from LAS object
#' llas <- list()
#' llas[["A"]] <- lidR::clip_circle(laschablais3, 974350, 6581680, 10)
#' llas[["B"]] <- lidR::clip_circle(laschablais3, 974390, 6581680, 10)
#' llas[["C"]] <- lidR::clip_circle(laschablais3, 974350, 6581640, 10)
#' # normalize point clouds
#' llas <- lapply(llas, function(x) {lidR::normalize_height(x, lidR::tin())})
#' 
#' # compute metrics
#' cloudMetrics(llas)
#' 
#' # compute metrics with user-defined function
#' # mean and standard deviation of first return points above 10 m
#' user.func <- function(z, rn, hmin = 10)
#' {
#' # first return above hmin subset
#' dummy <- which(z >= hmin & rn == 1)
#' return(list(
#'  mean.z = mean(z[dummy]),
#'  sd.z = stats::sd(z[z>hmin])
#'  ))
#'  }
#'  cloudMetrics(llas, func=~user.func(Z, ReturnNumber, 10))
#' @export
#'
cloudMetrics <- function(llasn,
                         func =  ~ lidR::stdmetrics(X, Y, Z, Intensity, ReturnNumber, Classification, dz = 1))
  {
    # apply lidR::lasmetrics to compute metrics
    metrics <-
      lapply(llasn, function(x) {
        lidR::cloud_metrics(x, func)
      })
    # identify elements with no metrics computed
    dummy <- which(unlist(lapply(metrics, function(x) {
      is.null(x)
    })))
    # combine results in a data.frame
    metrics <- lapply(metrics, function(x) {
      rbind(unlist(x))
    })
    metrics <- do.call(rbind, metrics)
    metrics <- data.frame(metrics)
    # add row names corresponding to list names, while removing elements with no metrics computed
    if (length(dummy) > 0)
    {
      row.names(metrics) <- names(llasn)[-dummy]
      warning("Some elements with no computed metrics")
    } else {
      row.names(metrics) <- names(llasn)
    }
    return(metrics)
  }
# ##############################################################
#' Function for area-based metrics computation
#' 
#' Predefined function usable in \code{\link[lidR]{cloud_metrics}} or \code{\link{cloudMetrics}}. Applies a minimum height threshold to the point cloud and computes the following metrics:
#' \enumerate{
#' \item for all points: total number \code{ntot}, percentage of points above minimum height \code{p.hmin}, percentage of points in height bins \code{H.propZ1_Z2},
#' \item for first return points: percentage above minimum height \code{p.1st.hmin}, 
#' \item for all points above minimum height: height metrics returned by \code{\link[lidR]{stdmetrics_z}} and intensity metrics returned by \code{\link[lidR]{stdmetrics_i}}
#' \item for first returns above minimum height: \code{mCH} and \code{sdCH} as proposed by Bouvier et al.
#' }
#'
#' @param z,i,rn,c Height, Intensity, ReturnNumber and Classification
#' @param hmin numeric. height threshold for low points removal before metrics computation
#' @param breaksH vector. breaks for height histogram proportion computation
#' @references Bouvier et al. 2015. Generalizing predictive models of forest inventory attributes using an area-based approach with airborne LiDAR data. Remote Sensing of Environment 156, pp. 322-334. \doi{10.1016/j.rse.2014.10.004}
#' @seealso \code{\link[lidR]{cloud_metrics}}, \code{\link[lidR]{stdmetrics}}, \code{\link{cloudMetrics}}
#' @examples
#' data(laschablais3)
#' 
#' # extract four point clouds from LAS object
#' llas <- list()
#' llas[["A"]] <- lidR::clip_circle(laschablais3, 974350, 6581680, 10)
#' llas[["B"]] <- lidR::clip_circle(laschablais3, 974390, 6581680, 10)
#' llas[["C"]] <- lidR::clip_circle(laschablais3, 974350, 6581640, 10)
#' # normalize point clouds
#' llas <- lapply(llas, function(x) {lidR::normalize_height(x, lidR::tin())})
#' 
#' # compute metrics
#' cloudMetrics(llas, ~ABAmodelMetrics(
#' Z, Intensity, ReturnNumber, Classification, 2, c(-Inf, 0, 2, 10, 20, +Inf)))
#' @rdname ABAmodelMetrics
#' @export
ABAmodelMetrics <- function(z, i, rn, c, hmin = 2, breaksH=NULL) 
{
  # first return above hmin subset
  dummy <- which(z >= hmin & rn == 1)
  # above hmin subset
  dummy2 <- which(z >= hmin)
  # non-null output only if points are present above hmin
  if (length(dummy2)>0)
  {
    metrics <- list(
      mCH = mean(z[dummy]), # mu_CH of Bouvier et al.
      sdCH = stats::sd(z[dummy]), # sigma_CH of Bouvier et al
      ntot = length(z), # total point number
      p.1st.hmin = length(dummy) / sum(rn==1), # percentage of 1st point above min height
      p.hmin = length(dummy2) / length(z) # percentage of points above min height
    )
    #
    if (!is.null(breaksH)) # strata relative point count
    {
      strata <- as.list(graphics::hist(z, breaks=breaksH, right=F, plot=F)$counts/length(z))
      names(strata) <- gsub("-","",paste("H.prop", breaksH[c(-length(breaksH))],"_", breaksH[c(-1)],sep=""))
    } else {strata <- NULL}
    return(c(lidR::stdmetrics_z(z[dummy2]), lidR::stdmetrics_i(i[dummy2], z=z[dummy2], class=NULL, rn=rn[dummy2]), metrics, strata))
  } else {return(NULL)}
}

#' @rdname ABAmodelMetrics
#' @export
.ABAmodelMetrics <- ~ ABAmodelMetrics(Z, Intensity, ReturnNumber, Classification, hmin=2, breaksH=NULL)

################################
#' Computation of tree metrics
#' 
#' This function computes summary statistics from a data.frame containing tree-level information as returned by \code{\link{treeExtraction}}.
#'
#' @param x data.frame containing the following columns for each line (segmented tree): \code{h} (height), \code{s} (crown surface), \code{v} (crown volume), typically returned by \code{\link{treeExtraction}}. \code{sp} (crown surface inside region of interest) and \code{vp} (crown volume in region of interest) are not used in this function.
#' @param area.ha numeric. area of region of interest in ha
#' @return a data.frame with one line containing the following tree metrics:
#' \enumerate{
#' \item \code{Tree.meanH}: mean height of detected tree apices (m)
#' \item \code{Tree.sdH}: standard deviation of heights of detected tree apices (m)
#' \item \code{Tree.giniH}: Gini index of heights of detected tree apices
#' \item \code{Tree.density}: density of detected tree apices (/ha)
#' \item \code{TreeInf10.density}: density of detected trees apices with h<=10 (/ha)
#' \item \code{TreeSup10.density}: density of detected trees apices with h>10 (/ha)
#' \item \code{TreeSup20.density}: density of detected trees apices with h>20 (/ha)
#' \item \code{TreeSup30.density}: density of detected trees apices with h>30 (/ha)
#' \item \code{Tree.meanCrownSurface}: mean crown surface of detected trees
#' \item \code{Tree.meanCrownVolume}: mean volume of detected trees
#' \item \code{TreeCanopy.meanH}: mean height of union of crowns of detected trees
#' }
#' @seealso \code{\link{treeExtraction}}
#' @examples
#' # sample 50 height values
#' h <- runif(50, 5, 40)
#' # simulate tree data.frame
#' trees <- data.frame(h=h, s=h, sp=h*0.95, v=h*h*0.6, vp=h*h*0.55)
#' stdTreeMetrics(trees, area.ha=0.1)
#' @export
#'
stdTreeMetrics <- function(x, area.ha = NA)
{
  data.frame(
    Tree.meanH = mean(x$h),
    Tree.sdH = stats::sd(x$h),
    Tree.giniH = ifelse(is.null(x), NA, reldist::gini(x$h)),
    Tree.density = length(x$h) / area.ha,
    TreeInf10.density = sum(x$h <= 10) / area.ha,
    TreeSup10.density = sum(x$h > 10) / area.ha,
    TreeSup20.density = sum(x$h > 20) / area.ha,
    TreeSup30.density = sum(x$h > 30) / area.ha,
    Tree.meanCrownSurface = mean(x$s),
    Tree.meanCrownVolume = mean(x$v),
    TreeCanopy.meanH = sum(x$v) / sum(x$s)
  )
}
#
###############################
#' Computation of terrain metrics
#' 
#' This function computes topographic variables from a point cloud
#' \itemize{
#' \item{exposition}
#' \item{altitude}
#' \item{slope}
#' }
#' values are computed after fitting a plane to the points. It supposes a homogeneous sampling of the plot by points. Points can be cropped on disk if center and radius are provided. In case a centre is provided, the altitude is computed by bilinear interpolation at the center location (\code{\link[lidR]{grid_metrics}} with \code{\link[lidR]{tin}} algorithm), otherwise it is the mean of the points altitude range.
#'
#' @param p matrix, data.frame or \code{\link[lidR]{LAS}} object with ground point coordinates (X, Y, Z). In case of an object which is not LAS, the object is first converted, which issues a warning
#' @param centre vector. x y coordinates of center to extract points inside a disc
#' @param r numeric. radius of disc
#' @return a data.frame with altitude, exposition (gr), slope (gr) and adjR2 of plane fitting
#' @examples
#' # sample points
#' XYZ <- data.frame(x=runif(200, -10, 10), y=runif(200, -10, 10))
#' XYZ$z <- 350 + 0.3 * XYZ$x + 0.1 * XYZ$y + rnorm(200, mean=0, sd=0.5)
#' # compute terrain statistics
#' points2terrainStats(XYZ)
#' points2terrainStats(XYZ, centre=c(5,5), r=5)
#' # with a LAS object
#' data(laschablais3)
#' terrain.points <- lidR::filter_ground(laschablais3)
#' points2terrainStats(terrain.points)
#' points2terrainStats(terrain.points, centre=c(974360, 6581650), r=10)
#' @export
#'
points2terrainStats <- function(p, centre=NULL, r=NULL)
{
  # if LAS object, extract coordinates
  if (!inherits(p,"LAS"))
  {
    # create LAS file from points
    p <- lidR::LAS(data.frame(X = p[, 1], Y=p[, 2], Z=p[, 3], Classification=2L),
                       check = FALSE)
  }
  # extract points inside disc if information is provided
  if (!is.null(r) & !is.null(centre))
  {
    # data.frame
    # p <- p[(p$X-as.numeric(centre[1]))^2+(p$Y-as.numeric(centre[2]))^2<=r^2,]
    # point cloud if present
    p <- lidR::clip_circle(p,
                           as.numeric(centre[1]),
                           as.numeric(centre[2]),
                           r)
    }
  # compute statistics if enough points are present
  if (nrow(p@data)<=1) return(NULL)
  # compute plane equation
  modlin <- stats::lm(Z~X+Y,data=p@data)
  # model z=a+bx+cy
  # normal vector: b c -1 (under plane)
  a <- modlin$coefficients[1]
  b <- modlin$coefficients[2]
  c <- modlin$coefficients[3]
  # compute slope
  slope <- atan(sqrt(b^2+c^2))*400/(2*pi)
  # compute azimuth if slope is not zero
  if (abs(slope)>0)
  {
    azimut <- ((pi/2-atan2(c,b))*400/(2*pi)+200)%%400
  } else {
    azimut <- NA
  }
  # extract altitude at centre if provided
  # else output mean of range
  if (!is.null(centre))
  {
    # use grid_terrain function to perform bilinear interpolation on one point
    # create raster with one cell at location of interest
    dummyRaster <- raster::raster(raster::extent(as.numeric(centre[1])-0.5, as.numeric(centre[1])+0.5,
                                                 as.numeric(centre[2])-0.5, as.numeric(centre[2])+0.5))
    raster::res(dummyRaster) <- 1
    altitude <- raster::values(lidR::grid_terrain(p, dummyRaster, lidR::tin()))
  } else {altitude <- NA}
  if (is.na(altitude)) {altitude <- mean(range(p@data$Z))}
  # output
  round(data.frame(altitude=altitude, azimut.gr=azimut, slope.gr=slope,
                   adjR2.plane=summary(modlin)$adj.r.squared * 100),1)
}
################################
#' Computes metrics on trees detected in list of point clouds.  
#' 
#' Extracts summary statistics on trees for each LAS object in a list:
#' 
#' \itemize{
#' \item{calls \code{\link{treeSegmentation}} to segment trees and then \code{\link{treeExtraction}} to extract their features}
#' \item{computes `TreeCanopy.coverInPlot` (proportion of surface of disk of interest which is covered by segmented trees), `TreeCanopy.meanHeightInPlot` (mean canopy height inside intersection of tree segments and disk of interest)}
#' \item{removes detected trees located outside of the disk of interest defined by their centers and radius}
#' \item{computes summary statistics of extracted tree features based on a user-defined function (default is \code{\link{stdTreeMetrics}})}
#' }
#'
#' @param llasn list of \code{\link[lidR]{LAS}} objects
#' @param XY a dataframe or matrix with XY coordinates of plot centers
#' @param plot.radius numeric. plot radius in meters
#' @param res numeric. resolution of canopy height model computed with \code{\link{points2DSM}} before tree segmentation
#' @param func a function to be applied to the attributes of extracted trees (return from internal call to \code{\link{treeExtraction}} function) to compute plot level metrics
#' @param ... other parameters to be passed to \code{\link{treeSegmentation}}
#' @return a dataframe with tree metrics in columns corresponding to LAS objects of the list (lines)
#' @seealso \code{\link{treeSegmentation}}, \code{\link{treeExtraction}}, \code{\link{stdTreeMetrics}}
#' @examples
#' data(laschablais3)
#' 
#' # extract three point clouds of 10 m radius from LAS object
#' llas <- list()
#' llas[[1]] <- lidR::clip_circle(laschablais3, 974350, 6581680, 10)
#' llas[[2]] <- lidR::clip_circle(laschablais3, 974390, 6581680, 10)
#' llas[[3]] <- lidR::clip_circle(laschablais3, 974350, 6581640, 10)
#' # normalize point clouds
#' llas <- lapply(llas, function(x) {lidR::normalize_height(x, lidR::tin())})
#' 
#' # compute tree metrics restricted to disks of radius 8 m.
#' cloudTreeMetrics(llas, 
#'                  cbind(c(974350, 974390, 974350), c(6581680, 6581680, 6581640)),
#'                  8, res=0.5)
#' 
#' # compute metrics with user-defined function
#' # number of detected trees between 20 and 30 meters and their mean height
#' user.func <- function(x)
#' {
#'   dummy <- x$h[which(x$h>20 & x$h<30)]
#'   data.frame(Tree.between.20.30=length(dummy), Tree.meanH=mean(dummy))
#' }
#'  cloudTreeMetrics(llas,
#'                   cbind(c(974350, 974390, 974350), c(6581680, 6581680, 6581640)),
#'                   8, res=0.5, func=user.func)
#' @export
#'
cloudTreeMetrics <- function(llasn, XY, plot.radius, res=0.5, func, ...)
{
  plot.area.ha <- pi*plot.radius^2/10000
  if (missing(func))
  {
    func <- function(x){stdTreeMetrics(x, area.ha=plot.area.ha)}
  }
  xy.list <- split(XY, seq(nrow(XY)))
  # LOOP on list <- replace by lapply
  ltrees <- list()
  lcanopy <- list()
  # add row names
  if (is.null(names(llasn))) names(llasn) <- as.character(1:length(llasn))
  #
  for (i in 1:length(llasn))
  {
    x <- llasn[[i]]
    coord <- xy.list[[i]]
    # compute dsm
    dummy <- points2DSM(x,res=res)
    # replace NA, low and high values
    dummy[is.na(dummy) | dummy<0] <- 0
    # tree detection
    dummy <-  treeSegmentation(dummy, ...)
    # compute mask of area of interest
    mask <- rasterXYMask(coord, plot.radius, dummy$local.maxima, binary = TRUE)
    # tree extraction
    ltrees[[names(llasn)[i]]] <- treeExtraction(dummy$filled.dem, dummy$local.maxima, dummy$segments.id, mask)
    # compute tree canopy cover fraction and mean height inside area of interest
    TreeCanopy.coverInPlot <- sum(raster::values((dummy$segments.id>0)*mask), na.rm=TRUE)/sum(raster::values(mask))
    mask[mask==0 | dummy$segments.id == 0] <- NA
    TreeCanopy.meanHeightInPlot <- mean(raster::values(dummy$filled.dem *mask), na.rm=TRUE)
    lcanopy[[names(llasn)[i]]] <- data.frame(TreeCanopy.coverInPlot, TreeCanopy.meanHeightInPlot)
  }
  #
  if (length(ltrees)>0)
  {
    tree.metrics  <- lapply(ltrees, FUN=func)
    tree.metrics <- do.call(rbind,tree.metrics)
    treeCanopy.metrics <- do.call(rbind,lcanopy)
    # merge datasets
    # create merge column
    treeCanopy.metrics$merge <- row.names(treeCanopy.metrics)
    tree.metrics$merge <- row.names(tree.metrics)
    # merge data.frames
    tree.metrics <- merge(tree.metrics, treeCanopy.metrics, all = TRUE)
    # row.names
    row.names(tree.metrics) <- tree.metrics$merge
    # remove merge column
    tree.metrics$merge <- NULL
  } else {tree.metrics <- NULL}
  return(tree.metrics)
}
