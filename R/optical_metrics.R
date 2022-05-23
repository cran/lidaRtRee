# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
#-------------------------------------------------------------------------------
#' Add vegetation indices on a IRC image
#'
#' Computes vegetations indices from the Red, Green and Infra-Red bands of an IRC
#' image and adds them as additionnal bands or columns. Indices are listed on
#' \url{https://www.l3harrisgeospatial.com/docs/broadbandgreenness.html}
#'
#' @param r raster or data.frame. Should contain bands or columns with
#' names nir, r, g
#' @param all boolean. indicates whether all indices should be computed;
#' default:FALSE, only grvi, sr and ndvi are calculated
#' @return a raster or data.frame with added bands or columns
#' @examples
#' df <- data.frame(nir = c(110, 150, 20), r = c(25, 50, 30), g = c(10, 60, 10))
#' add_vegetation_indices(df, all = TRUE)
#' @export
add_vegetation_indices <- function(r, all = FALSE) {
  # Green Ratio Vegetation Index (GRVI)
  r$grvi <- r$nir / r$g
  # Simple ratio (SR)
  r$sr <- r$nir / r$r
  # Normalized Difference Vegetation Index (NDVI)
  r$ndvi <- (r$nir - r$r) / (r$nir + r$r)
  #
  # additional indices
  if (all) {
    # Difference Vegetation Index (DVI)
    r$dvi <- r$nir - r$r
    # Green Difference Vegetation Index (GDVI)
    r$gdvi <- r$nir - r$g
    # Infrared Percentage Vegetation Index (IPVI)
    r$ipvi <- r$nir / (r$nir + r$r)
    # modified non-linear index (MNLI)
    r$mnli <- (r$nir^2 - r$r) * (1 + 0.5) / (r$nir^2 + r$r + 0.5)
    # Modified Simple Ratio (MSR)
    r$msr <- ((r$nir / r$r) - 1) / (sqrt(r$nir / r$r) - 1)
    # Non-Linear Index (NLI)
    r$nli <- (r$nir^2 - r$r) / (r$nir^2 + r$r)
    # Optimized Soil Adjusted Vegetation Index (OSAVI)
    r$osavi <- (1.5 * (r$nir - r$r)) / (r$nir + r$r + 0.16)
    # Renormalized Difference Vegetation Index (RDVI)
    r$rdvi <- (r$nir - r$r) / sqrt(r$nir + r$r)
    # Soil Adjusted Vegetation Index (SAVI)
    r$savi <- (1.5 * (r$nir - r$r)) / (r$nir + r$r + 0.5)
    # Transformed Difference Vegetation Index (TDVI)
    r$ndvi <- sqrt(0.5 + (r$nir - r$r) / (r$nir + r$r))
  }
  r
}
