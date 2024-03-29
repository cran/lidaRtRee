% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aba.R
\name{aba_combine_strata}
\alias{aba_combine_strata}
\title{Combines a list of ABA models into a single ABA model object}
\usage{
aba_combine_strata(model.list, plotsId = NULL)
}
\arguments{
\item{model.list}{list. stratum-specific models returned by
\code{\link{aba_build_model}}}

\item{plotsId}{vector. "plotsId" for ordering row names in the "values" element
of the output list}
}
\value{
a list with three elements
\itemize{
\item \code{model}: a list of regression models corresponding to each stratum
 (output from \code{\link[stats]{lm}}),
\item \code{stats}: model statistics of each stratum-specific model (as in
\code{\link{aba_build_model}}) plus one line corresponding to statistics for all
strata (COMBINED)
\item \code{values}: data.frame with observed and values predicted in
cross-validation, and information on which stratum it belongs to.
}
}
\description{
Combines a list of models (obtained with \code{\link{aba_build_model}}) into a
single object. Typically used to merge stratum-specific models into one object.
Validation statistics are computed for the combined strata, making it easier
to compare prediction performance with an unstratified model.
}
\examples{
# load Quatre Montagnes dataset
data(quatre_montagnes)
# initialize list of models
model_aba_stratified <- list()
# calibrate basal area prediction model for each stratum
for (i in levels(quatre_montagnes$stratum))
{
  subsample <- which(quatre_montagnes$stratum == i)
  model_aba_stratified[[i]] <-
    aba_build_model(quatre_montagnes[subsample, "G_m2_ha"],
      quatre_montagnes[subsample, 9:76],
      transform = "boxcox", nmax = 4,
      xy = quatre_montagnes[subsample, c("X", "Y")]
    )
}
# combine models in single object
model_aba_stratified <- aba_combine_strata(
  model_aba_stratified,
  quatre_montagnes$plotId
)
# display content of output list
model_aba_stratified$model
model_aba_stratified$stats
summary(model_aba_stratified$values)

# plot field values VS predictions in cross-validation
aba_plot(model_aba_stratified)
}
\seealso{
\code{\link{aba_build_model}} for calibrated ABA model,
\code{\link{aba_plot}} for plotting model cross-validation results.
}
