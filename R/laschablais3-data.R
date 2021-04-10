#' las data in France (Chablais 3 plot)
#'
#' @description Airborne laser scanning data over the Chablais 3 plot, acquired in 2009 by Sintegra, copyright INRAE
#' 
#' @details Additional information about the data
#' \itemize{
#' \item{Sensor: RIEGL LMS-Q560)}
#' \item{EPSG code of coordinates system: 2154}
#'}
#'
#' @docType data
#'
#' @usage data(laschablais3)
#'
#' @format An object of class \code{\link[lidR]{LAS}}
#'
#' @keywords datasets
#'
#' @references Monnet, J.-M. 2011. Using airborne laser scanning for mountain forests mapping: Support vector regression for stand parameters estimation and unsupervised training for treetop detection. Ph.D. thesis. University of Grenoble, France. pp. 21-22. \url{https://tel.archives-ouvertes.fr/tel-00652698/document}
#'
#' @source Monnet J.-M. INRAE
#'
#' @examples
#' data(laschablais3)
#' laschablais3
#' # display point cloud
#' lidR::plot(laschablais3)
#' @name laschablais3
NULL
"laschablais3"
