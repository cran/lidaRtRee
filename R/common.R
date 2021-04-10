# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
################################
#' Digital Surface Model
#' 
#' Creates a Digital Surface Model from LAS object. Raster extent is specified by the coordinates of lower left and upper right corners. Default extent covers the full range of points, and aligns on multiple values of the resolution. Cell value is the maximum height of points contained in the cell.
#'
#' @param .las \code{\link[lidR]{LAS}} object or XYZ matrix/data.frame
#' @param res numeric. raster resolution
#' @param xmin numeric. lower left corner easting coordinate for output raster.
#' @param xmax numeric. upper right corner easting coordinate for output raster.
#' @param ymin numeric. lower left corner northing coordinate for output raster.
#' @param ymax numeric. upper right corner northing coordinate for output raster.
#' @return A raster object.
#' @seealso \code{\link{points2DTM}} for Digital Terrain Model computation.
#' @examples
#' data(laschablais3)
#' 
#' # create a digital surface model with first-return points, resolution 0.5 m
#' dsm <- points2DSM(lidR::filter_first(laschablais3), res=0.5)
#' 
#' # display raster
#' raster::plot(dsm,asp=1)
#' @export
points2DSM <- function(.las, res=1, xmin, xmax, ymin, ymax)
{
  # in case input is not a las object
  if (class(.las)[1] != "LAS")
  {
    .las <- lidR::LAS(data=data.frame(X=.las[,1], Y=.las[,2], Z=.las[,3]))
  }
  #
  if (missing(xmin) | missing(xmax) | missing(ymin) | missing(ymax)) # if no extent info
  {
    xmin <- floor(min(.las@data$X)/res)*res
    xmax <- ceiling(max(.las@data$X)/res)*res
    ymin <- floor(min(.las@data$Y)/res)*res
    ymax <- ceiling(max(.las@data$Y)/res)*res
  }
  # create empty raster
  r <- raster::raster()
  raster::extent(r) <- c(xmin,xmax,ymin,ymax)
  raster::res(r) <- res
  raster::crs(r) <- .las@proj4string
  # convert LAS coordinates to spatial data
  points <- as.data.frame(.las@data[,1:2])
  val <- .las@data[,3]
  sp::coordinates(points)=c(1,2)
  points@proj4string <- .las@proj4string
  # rasterize with max function
  raster::rasterize(points,r,val,fun=max)
}
###################################
#' Digital Terrain Model
#' 
#' Creates a Digital Terrain Model from LAS object or XYZ data. Raster extent is specified by the coordinates of lower left and upper right corners. Default extent covers the full range of points, and aligns on multiple values of the resolution. Cell value is compute as the bilinear interpolation at the cell center form an Delaunay triangulation. Relies on \code{\link[lidR]{grid_terrain}} with algorithm \code{\link[lidR]{tin}}. In case a LAS object is provided, only points classified as ground or water (2 or 9) will be used; a warning is issued when the raster output inherits the projection info from LAS input.
#'
#' @param .las \code{\link[lidR]{LAS}} object or XYZ matrix/data.frame containing only ground points
#' @param res numeric. raster resolution
#' @param xmin numeric. lower left corner easting coordinate for output raster.
#' @param xmax numeric. upper right corner easting coordinate for output raster.
#' @param ymin numeric. lower left corner northing coordinate for output raster.
#' @param ymax numeric. upper right corner northing coordinate for output raster.
#' @seealso \code{\link{points2DSM}} for Digital Surface Model computation.
#' @return A raster object
#' @examples
#' data(laschablais3)
#' 
#' # create digital terrain model with points classified as ground
#' dtm <- points2DTM(laschablais3)
#' 
#' # display raster
#' raster::plot(dtm,asp=1)
#' @export
points2DTM = function(.las, res = 1, xmin, xmax, ymin, ymax)
{
  # in case input is not a las object
  if (class(.las)[1] != "LAS")
  {
    .las <- lidR::LAS(data=data.frame(X=.las[,1], Y=.las[,2], Z=.las[,3], Classification=2L))
  }
  #
  if (missing(xmin) | missing(xmax) | missing(ymin) | missing(ymax)) # if no extent info
  {
    xmin <- floor(min(.las@data$X)/res)*res
    xmax <- ceiling(max(.las@data$X)/res)*res
    ymin <- floor(min(.las@data$Y)/res)*res
    ymax <- ceiling(max(.las@data$Y)/res)*res
  }

    #points <- as.matrix(.las@data[,c("X","Y","Z")])
    # interpolation: value estimated at pixel center by bilinear interpolation
    #mnt <- akima::interp(points[,1],points[,2],points[,3], xo=seq(xmin+res/2, xmax-res/2,by=res), yo=seq(ymax-res/2,ymin+res/2,by=-res),linear=TRUE,extrap=FALSE, duplicate="user", dupfun=function(x){min(x)})
    #
  # dtm = raster::raster(t(mnt[[3]]), xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, crs=.las@proj4string)
  # create target raster
  dtm <- raster::raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, resolution = res)
  dtm <- lidR::grid_terrain(.las, dtm, lidR::tin())
  raster::crs(dtm) <- sp::proj4string(.las)
  return(dtm)
}
# ################################
# #' creates Digital Elevation Model from LAS object and extent info (parallelized version)
# #' tiles extent and calls points2DTM or points2DSM
# #'
# #' @param .las LAS object (from \code{lidR} package)
# #' @param type model to return ("dsm" or "dtm")
# #' @param res the raster resolution
# #' @param xmin the raster lower left corner easting coordinate
# #' @param xmax the raster upper right corner easting coordinate
# #' @param ymin the raster lower left corner northing coordinate
# #' @param ymax the raster upper right corner northing coordinate
# #' @param method a string specifying the interpolation method
# #' @param buffer the size of buffer area for tiling
# #' @param tile.size the size of tiles for parallelisation
# #' @param n.cores the number of cores to use
# #' @return a raster object
# #' @export
# points2DEM = function(.las, type="dtm", res = 1, xmin, xmax, ymin, ymax, method="delaunay", buffer=5, tile.size=250, n.cores=4)
# {
#   doParallel::registerDoParallel(cores=n.cores)
#   if (missing(xmin) | missing(xmax) | missing(ymin) | missing(ymax)) # if no extent info
#   {
#     xmin <- floor(min(.las@data$X)/res)*res
#     xmax <- ceiling(max(.las@data$X)/res)*res
#     ymin <- floor(min(.las@data$Y)/res)*res
#     ymax <- ceiling(max(.las@data$Y)/res)*res
#   }
#   # low left corner
#   xlow <- floor(xmin/tile.size)*tile.size
#   ylow <- floor(ymin/tile.size)*tile.size
#   # tiles indices which intersect the extent
#   i <- 0
#   j <- 0
#   n <- 1
#   l <- list()
#   while(xlow+i*tile.size<xmax)
#   {
#     while(ylow + j * tile.size < ymax)
#     {
#       l[[n]] <- c(i,j)
#       n <- n+1
#       j <- j+1
#     }
#     i <- i+1
#     j <- 0
#   }
#   l <- foreach::foreach(i=1:length(l), .errorhandling="remove") %dopar%
#   {
#     coord <- c(max(xmin, xlow+l[[i]][1]*tile.size), min(xmax, xlow+(l[[i]][1]+1)*tile.size), max(ymin, ylow+l[[i]][2]*tile.size), min(ymax, ylow+(l[[i]][2]+1)*tile.size))
#     dummy <- .las %>% lidR::lasfilter(X >= coord[1]-buffer & X <= coord[2]+buffer & Y >= coord[3]-buffer & Y <= coord[4]+buffer)
#     if (type=="dtm") {return(points2DTM(dummy, res = res, coord[1], coord[2], coord[3], coord[4], method=method))}
#     else {return(points2DSM(dummy, res = res, coord[1], coord[2], coord[3], coord[4]))}
#   }
#   l$fun <- mean
#   l$na.rm <- TRUE
#   do.call(raster::mosaic, l)
# }
################################
#' Polar to cartesian coordinates conversion
#' 
#' Computes projected coordinates (Easting, Northing, Altitude) from polar coordinates (Azimuth, Slope, Distance) and center position (Easting, Northing, Altitude). Magnetic declination and meridian convergence are optional parameters. In case distance is measured to the border of objects (e.g. trees), the diameter can be added to compute the coordinates of object center.
#'
#' @param x vector. easting coordinates of centers in meter
#' @param y vector. northing coordinates of centers in meter
#' @param z vector. altitudes of centers in meters
#' @param declination vector. magnetic declination values in radian
#' @param convergence vector. meridian convergence values in radian
#' @param azimuth vector. azimuth values from centers in radian
#' @param slope vector. slope values from centers in radian
#' @param dist vector. distances between centers and objects in meter
#' @param diameter vector. diameters in meter (e.g. in case a radius should be added to the distance)
#' @seealso \code{\link{plotTreeInventory}} for tree inventory display
#' @return A data.frame with easting, northing and altitude coordinates, and horizontal distance from centers to objects centers
#' @examples
#' # create data.frame of trees with polar coordinates and diameters
#' trees <- data.frame (x=rep(c(0,10), each=2),
#'  y = rep(c(0,10),each=2),
#'  z=rep(c(0,2),each=2),
#'  azimuth=rep(c(0,pi/3)),
#'  dist=rep(c(2,4)),
#'  slope=rep(c(0,pi/6)),
#'  diameter.cm=c(15,20,25,30))
#' trees
#'  
#'  # compute projected coordinates
#' polar2Projected(trees$x, trees$y, trees$z, trees$azimuth, trees$dist,
#'  trees$slope, declination=0.03, convergence=0.02, trees$diameter.cm/100)
#' @export
polar2Projected <- function(x, y, z=0, azimuth, dist, slope=0, declination=0, convergence=0, diameter=0)
{
  # compute horizontal distance
  d <- dist*cos(slope)+diameter/2
  data.frame(x = x + d*sin(azimuth+convergence+declination),
             y = y + d*cos(azimuth+convergence+declination),
             z = z + dist*sin(slope),
             d)
}
###############################
#' Table of species names, abreviations and display colors
#' 
#' table for species names, abreviations and type (coniferous/broadleaf), and display color
#' @return A data frame with species name, color, coniferous (C) / broadleaf (B) type, and name abreviation GESP of GEnus and SPecies 
#' @seealso \code{\link{plotTreeInventory}} for tree inventory display
#' @examples
#' # load table
#' tab.species <- speciesColor()
#' head(tab.species)
#' summary(tab.species)
#' @export
speciesColor <- function()
{
  d <- rbind(c("Abies alba","purple", "C"),
             c("Acer","orange", "B"),
             c("Acer campestre","orange1", "B"),
             c("Acer opalus","orange2", "B"),
             c("Acer platanoides","orange3", "B"),
             c("Acer pseudoplatanus","orange4", "B"),
             c("Acer sp.","orange", "B"),
             c("Alnus sp.","violet", "B"),
             c("Alnus incana","violet", "B"),
             c("Alnus viridis","violet", "B"),
             c("Betula pendula","darkgreen", "B"),
             c("Betula pubescens","darkgreen", "B"),
             c("Betula sp.","darkgreen", "B"),
             c("Buxus sempervirens","black", "B"),
             c("Carpinus betulus","seagreen", "B"),
             c("Castanea sativa","pink", "B"),
             c("Corylus avellana","cadetblue2", "B"),
             c("Cornus mas","black", "B"),
             c("Cotinus sp.","black", "B"),
             c("Crataegus sp.","black", "B"),
             c("Crataegus monogyna","black", "B"),
             c("Euonymus latifolius","black", "B"),
             c("Fagus sylvatica","green", "B"),
             c("feuillus","black", "B"),
             c("Fraxinus excelsior","yellow", "B"),
             c("Ilex aquifolium","cyan", "B"),
             c("inconnu","black", NA),
             c("Juglans regia", "black", "B"),
             c("Juniperus communis","black", "C"),
             c("Juniperus sp.","black", "C"),
             c("Laburnum anagyroides","black", "B"),
             c("Larix decidua","pink", "C"),
             c("Larix kaempferi","pink", "C"),
             c("Malus sylvestris","plum", "B"),
             c("Picea abies","blue", "C"),
             c("Pinus cembro","salmon3", "C"),
             c("Pinus mugo","salmon2", "C"),
             c("Pinus nigra","salmon1", "C"),
             c("Pinus sp.","salmon1", "C"),
             c("Pinus uncinata","salmon4", "C"),
             c("Pinus sylvestris","salmon", "C"),
             c("Populus alba","slateblue", "B"),
             c("Populus nigra","slateblue", "B"),
             c("Populus tremula","slateblue", "B"),
             c("Populus sp.","slateblue", "B"),
             c("Pseudotsuga menziesii","darkblue", "C"),
             c("Prunus avium","grey", "B"),
             c("Prunus sp.","grey", "B"),
             c("Pyrus communis","grey", "B"),
             c("Quercus sp.","turquoise", "B"),
             c("Quercus petraea","turquoise", "B"),
             c("Quercus pubescens","turquoise", "B"),
             c("Quercus robur","turquoise", "B"),
             c("Robinier pseudoacacia","pink", "B"),
             c("Salix caprea","darkgoldenrod2", "B"),
             c("Salix sp.","darkgoldenrod2", "B"),
             c("Salix nigra","darkgoldenrod2", "B"),
             c("Sorbus aria","red3", "B"),
             c("Sorbus aucuparia","red", "B"),
             c("Taxus baccata","burlywood3", "C"),
             c("Tilia cordata","chocolate4", "B"),
             c("Tilia platyphyllos","chocolate4", "B"),
             c("Tilia sp.","chocolate4", "B"),
             c("Ulmus glabra","brown", "B"),
             c("Ulmus sp.","brown", "B"),
             c("Viburnum lantana","black", "B"))
  d <- as.data.frame(d,stringsAsFactors = F)
  names(d) <- c("name","col","broad.conif")
  d$broad.conif <- factor(d$broad.conif)
  # abreviation GEsp (GEnus species)
  dummy <- strsplit(d$name," ")
  dummy2 <- lapply(dummy, function(x){
    ifelse(length(x)==2 && x[2]=="sp.", 
           paste(toupper(substr(x[1],1,2)),"sp",sep=""),
           paste(toupper(substr(x,1,2)),collapse=""))})
  d$abvr <- row.names(d) <- unlist(dummy2)
  d
}
#######################################
#' Displays a map of tree inventory data
#' 
#' displays tree inventory data
#'
#' @param xy data.frame. contains two columns with the X, Y  coordinates of tree centers
#' @param height vector.  tree heights in meters 
#' @param diam vector. tree diameters in centimeters
#' @param species vector. species abbreviation as in \code{\link{speciesColor}} for display with corresponding color
#' @param ... Arguments to be passed to methods, as in \code{\link[graphics]{plot}}
#' @seealso \code{\link{speciesColor}} for a table of species and associated colors
#' @examples
#' # load tree inventory data from plot Chablais 3
#' data("treeinventorychablais3")
#' 
#' # display tree inventory
#' plotTreeInventory(treeinventorychablais3[,c("x","y")],
#' diam = treeinventorychablais3$d, col = "red",
#' pch = treeinventorychablais3$e,
#' xlab = "X", ylab = "Y")
#' 
#' # display tree inventory with CHM background
#' data("chmchablais3")
#' raster::plot(chmchablais3, col=gray(seq(0,1,1/255)))
#' plotTreeInventory(treeinventorychablais3[,c("x","y")],
#' height = treeinventorychablais3$h,
#' species = treeinventorychablais3$s,
#' add = TRUE)
#' 
#' @return no return
#' @export
#'
plotTreeInventory <-function(xy, height=NULL, diam=NULL, species=NULL, ...)
{
  # retrieve dots
  dots <- list(...)
  if (!methods::hasArg("add")) dots$add <- FALSE
  # set size of symbol
  # proportional to height if present
  if(!is.null(height))
  {
    size <- height/20
    # otherwise proportional to diameter
  } else {
    if (!is.null(diam/40))
    {
      size <- diam/20
    } else { size <- 1}
  }
  # is user-specified cex is present
  if (!methods::hasArg("cex")) dots$cex = size
  # apply species-specific colors
  if (!is.null(species))
  {
    # load palette
    color <- speciesColor()
    # extract corresponding colors
    col1 <- color[as.character(species), c("col","abvr")]
    # add /overwrite in dots arguments
    dots$col <- col1$col
  } else { 
    if (!methods::hasArg("col")) # is user-specified colors are not present
    {
      dots$col <- "black"
    }
  }
  #
  if (dots$add==FALSE)
  {
    args <- list(x = xy[,1], y = xy[,2])
    # call plot with those arguments and those in dots except "add"
    dots$add <- NULL
    # if user-specified asp are not present
    if (!methods::hasArg("asp")) dots$asp <- 1
    # if user-specified xlab is not present
    if (!methods::hasArg("xlab")) dots$xlab = "Easting (m)"
    # if user-specified ylab is not present
    if (!methods::hasArg("ylab")) dots$ylab = "Northing (m)"
    #
    do.call(plot, c(args, dots))
  } else { 
    args <- list(x = xy[,1], y = xy[,2])
    # call plot with those arguments and those in dots except "add"
    dots$add <- NULL
    do.call(graphics::points, c(args, dots))
  }
  # add color legend
  if (!is.null(species))
  {
    texte <- sort(unique(col1$abvr))
    graphics::legend("topleft", texte, fill=color[texte,"col"], y.intersp = 1)
  }
}
##################################
#' Raster mask by union of buffers around xy positions
#' 
#' creates a raster mask by union of circular buffers around xy positions
#'
#' @param xy 2 columns matrix or data.frame. xy positions
#' @param buff vector.  buffers to apply to the xy positions
#' @param r raster object. target raster
#' @param binary boolean. should the output mask be boolean (TRUE) or greyscale (FALSE)
#' @return a raster object
#' @seealso \code{\link{rasterChullMask}}
#' @examples
#' # create raster
#' r <- raster::raster()
#' raster::extent(r) <- c(0,40,0,40)
#' raster::res(r) <- 1
#' 
#' # xy positions
#' xy <- data.frame(x=c(10,20,31.25,15),
#'   y=c(10,20,31.25,25))
#' # compute mask
#' mask1 <- rasterXYMask(xy, c(5, 8, 5, 5), r)
#' mask2 <- rasterXYMask(xy, c(5, 8, 5, 5), r, binary=FALSE)
#' 
#' # display binary raster
#' raster::plot(mask1)
#' graphics::points(xy)
#' 
#' # display distance raster
#' raster::plot(mask2)
#' graphics::points(xy)
#' @export
rasterXYMask <- function(xy, buff, r, binary=TRUE)
{
  # convert vector to data.frame in case only one pair of coodinates is provided
  if (is.null(dim(xy)) & length(xy)==2) {xy <- data.frame(matrix(xy, 1, 2))}
  # compute squared buffers
  buff2 <- buff^2
  # compute XY coordinates of cell centers
  dummyXY <- raster::xyFromCell(r,1:length(r))
  # create matrix of cell values
  val <- matrix(NA,nrow=length(r), ncol=nrow(xy))
  # compute cell value for each position
  for (i in 1:nrow(xy))
  {
    val[,i] <- sqrt(pmax(0,buff2[i] - ((dummyXY[,1]-xy[i,1])^2 + (dummyXY[,2]-xy[i,2])^2)))
  }
  # take maximum for each cell
  val <- apply(val,1,max)
  # convert to binary if required
  if (binary) {val <- val >0}
  raster::values(r) <- val
  r
}
####################################
#' Raster mask of convex hull
#' 
#' creates raster mask corresponding to the convex hull of xy positions
#'
#' @param xy 2 columns matrix or data.frame. xy positions
#' @param r raster object. target raster
#' @return a raster with 0 or 1
#' @examples
#' # create raster
#' r <- raster::raster()
#' raster::extent(r) <- c(0,40,0,40)
#' raster::res(r) <- 1
#' 
#' # xy positions
#' xy <- data.frame(c(10,20,31.25,15),
#'   c(10,20,31.25,25))
#' # compute mask
#' mask1 <- rasterChullMask(xy, r)
#' 
#' # display binary raster
#' raster::plot(mask1)
#' graphics::points(xy)
#' @seealso \code{\link{rasterXYMask}}
#' @export
rasterChullMask <- function(xy, r)
{
  # points on convexHull
  idchull <- grDevices::chull(xy)
  #
  # create chull polygon
  Po <- sp::Polygons(list(sp::Polygon(xy[c(idchull,idchull[1]),])),1)
  sPo <- sp::SpatialPolygons(list(Po),proj4string=r@crs)
  # rasterize polygon
  dummy <- raster::rasterize(sPo,r)
  dummy[is.na(dummy)] <- 0
  dummy
}
####################################
#' Create elliptical polygons from centres and extensions in four directions
#' 
#' creates polygons from the union of four quarters of ellipses, specified by the ellipse center, and maximum extension in two directions
#'
#' @param x,y vectors of numerics. Coordinates of ellipses centers
#' @param n,s,e,w vectors of numerics. Coordinates of ellipses extention in the north, south, east and west directions
#' @param id vector of strings. id of each polygon
#' @param step numeric. Angular step for the modelling of ellipses
#' @param angle.offset numeric. Angle offset to tilt ellipses, positive values rotates clockwise
#' @return a list of data.frame containing the coordinates of polygons
#' @examples
#' # compute coordinates of ellipses
#' ellipses1 <- ellipses4Crown(c(0,10), c(0,10), c(2,2), c(3,4), c(2.5,3), c(2,3),
#'                             id=c("A", "B"))
#' ellipses1[["A"]]
#' # tilted ellipse
#' ellipses2 <- ellipses4Crown(c(0,10), c(0,10), c(2,2), c(3,4), c(2.5,3), c(2,3),
#'                             angle.offset=pi/6)
#' ellipses2[[2]]
#' 
#' # draw ellipses in black, tilted ellipses in red
#' plot(ellipses1[[1]],type="l", asp=1, xlim=c(-5,15), ylim=c(-5,15))
#' lines(ellipses1[[2]])
#' lines(ellipses2[[1]], col="red")
#' lines(ellipses2[[2]], col="red")
#' @seealso \code{\link{pointList2SPDF}}
#' @export
ellipses4Crown <- function(x, y, n, s, e, w, id=NULL, step=pi/12, angle.offset=0)
{
  # create polygons by union of 4 quarter of ellipses
  # NE quadrant
  a <- seq(from=0,to=pi/2-step,by=step)
  X1 <- matrix(e, nrow=length(e)) %*% matrix(cos(a), ncol=length(a)) 
  Y1 <- matrix(n, nrow=length(n)) %*% matrix(sin(a), ncol=length(a)) 
  # NW quadrant
  a <- seq(from=pi/2,to=pi-step,by=step)
  X2 <- matrix(w, nrow=length(w)) %*% matrix(cos(a), ncol=length(a))
  Y2 <- matrix(n, nrow=length(n)) %*% matrix(sin(a), ncol=length(a))
  # SW quadrant
  a <- seq(from=pi,to=3*pi/2-step,by=step)
  X3 <- matrix(w, nrow=length(w)) %*% matrix(cos(a), ncol=length(a))
  Y3 <- matrix(s, nrow=length(s)) %*% matrix(sin(a), ncol=length(a))
  # SE quadrant
  a <- seq(from=3*pi/2,to=2*pi,by=step)
  X4 <- matrix(e, nrow=length(e)) %*% matrix(cos(a), ncol=length(a))
  Y4 <- matrix(s, nrow=length(s)) %*% matrix(sin(a), ncol=length(a))
  #
  X <- cbind(X1, X2, X3, X4)
  Y <- cbind(Y1, Y2, Y3, Y4)
  # create list of cordinates for each crown
  l <- list()
  # translate and rotation polygons
  for (i in 1:length(x))
  {
    l[[i]] <- cbind(X[i,], Y[i,]) %*% matrix(c(cos(angle.offset), -sin(angle.offset), sin(angle.offset), cos(angle.offset)), nrow=2, byrow=TRUE) + cbind(rep(x[i], times=ncol(X)), rep(y[i], times=ncol(Y)))
  }
  # add names
  if (!is.null(id)) {names(l) <- id}
  # remove polygons with NA values
  to.remove <- NULL
  for (i in 1:length(l))
  {
    if (!all(!is.na(l[[i]]))) {to.remove <- c(to.remove, i)}
  }
  if (is.null(to.remove))
  {
    l
  } else {
    l[-to.remove]
  }
}
##########################################################
#' Convert list of points into Spatial Polygons DataFrame object
#'
#' Converts a list of points specifying polygons into a Spatial Polygons DataFrame object
#'
#' @param points.list list of dataframes of xy coordinates. The first and last coordinates in each dataframe must be the same
#' @param df data.frame. Optional data.frame to be associated to Spatial Polygons
#' @param ... arguments to be passed to \code{\link[sp]{SpatialPolygons}}, e.g. proj4string info
#' @return an object of class \link[sp]{SpatialPolygons-class}, or class \link[sp]{SpatialPolygonsDataFrame-class} if input data.frame is specified.
#' @examples
#' # Compute coordinates of polygons
#' ellipses <- ellipses4Crown(c(0,10), c(0,10), c(2,2), c(3,4), c(2.5,3), c(2,3),
#'                            id=c("A", "B"))
#' # Convert to Spatial object
#' ellipses1 <- pointList2SPDF(ellipses)
#' ellipses1
#' # Convert to Spatial object with data.frame
#' ellipses2 <- pointList2SPDF(ellipses, df=data.frame(info=1:2),
#'                             proj4string = sp::CRS(SRS_string = "EPSG:2154"))
#' ellipses2
#' 
#' # draw ellipses
#' sp::plot(ellipses1)
#' sp::plot(ellipses2)
#' @seealso \code{\link{ellipses4Crown}}
#' @export
pointList2SPDF <- function(points.list, df=NULL, ...)
{
  # convert each element of list of coodinates to Polygon
  H <- lapply(points.list, FUN=function(x) {sp::Polygon(x, hole=F)})
  # convert to Polygons
  Hs <- lapply(H, function(x) {sp::Polygons(list(x), "temp")})
  # change ID
  if (!is.null(names(points.list))) 
  {
    for (i in 1:length(Hs)) {Hs[[i]]@ID <- names(points.list)[i]}
  }
  # convert to spatial polygons
  sp.H <- sp::SpatialPolygons(Hs, ...)
  if (!is.null(df))
  {
    sp::SpatialPolygonsDataFrame(sp.H, df, match.ID=FALSE)
  } else {
    sp.H
  }
}
##########################################################
# # plot tree inventory with ggplot2
# # # charger donnees chablais3
# data(treeinventorychablais3, package = "lidaRtRee")
# ggplot(treeinventorychablais3, aes(x = x, y = y, group = s)) +
#   geom_point(aes(color = s, size = d / 20)) +
#   coord_sf(datum = 2154)
# #
# # avec couleur personnalisée
# # charger la table de reference
# table.especes <- lidaRtRee::speciesColor()
# # extraire de la table les couleurs
# # correspondant aux especes presentes sur la placette
# table.placette <-
#   table.especes[levels(treeinventorychablais3$s), "col"]
# 
# ggplot(treeinventorychablais3, aes(x = x, y = y, group = s)) +
#   geom_point(aes(color = s, size = d)) +
#   coord_sf(datum = 2154) +
#   scale_color_manual(values = table.placette) +
#   scale_radius(name="Diamètre") + # taille de symbole proportionnelle au diamètre
#   geom_text(aes(label=n, size=20), hjust=0, vjust=1) +
#   labs(color = "Essence") # titre de la légende
# # table.especes[levels(treeinventorychablais3$s), "col"] 
# 
# # comparer le résultat avec ceci (éviter les surprises de permutation des lignes)
# lidaRtRee::plotTreeInventory(treeinventorychablais3[,1:2], diam = treeinventorychablais3$d, species=treeinventorychablais3$s)