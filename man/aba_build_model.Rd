% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aba.R
\name{aba_build_model}
\alias{aba_build_model}
\title{Calibrates and validates area-based models}
\usage{
aba_build_model(
  variable,
  predictors,
  transform = "none",
  nmax = 3,
  test = c("partial_p", "vif", "gvlma"),
  xy = NULL,
  threshold = NULL
)
}
\arguments{
\item{variable}{vector. dependent variable values}

\item{predictors}{data.frame. independent variables (columns: metrics,
lines: observations). Row names are used for the output predicted values}

\item{transform}{string. transformation to be applied to data (\code{"none"},
\code{"boxcox"}: Box-Cox transformation applied only to the dependent variable,
\code{"log"}: log transformation applied to both dependent and independent
variables)}

\item{nmax}{numeric. maximum number of independent variables in the model}

\item{test}{vector. which tests should be satisfied by the models, one to
three in \code{"partial_p"}, \code{"vif"}, \code{"gvlma"}}

\item{xy}{data.frame or matrix of easting and northing coordinates of
observations: not used in the function but exported in the result for use in
 further inference functions}

\item{threshold}{vector of length two. minimum and maximum values of threshold
to apply to predicted values}
}
\value{
a list with three elements
\itemize{
\item \code{model}: list with one regression model (output from
 \code{\link[stats]{lm}}),
\item \code{stats}: model statistics (root mean square error estimated in
leave-one-out cross validation, coefficient of variation of rmse, p-value of
wilcoxon test of observed and predicted values, p-value of t-test of observed
and predicted values, p-value of anova of observed and predicted values,
correlation of observed and predicted values, R2 of observed and predicted
values, variance of regression residuals)
\item \code{values}: data.frame with observed and values predicted in
cross-validation.
}
}
\description{
The function can first apply a Box-Cox transformation to the dependent variable,
in order to normalize its distribution, or a log transformation to the whole
dataset. Then it uses \code{\link[leaps]{regsubsets}} to find the 20 linear
regressions with the best adjusted-R2 among combinations of at most \code{nmax}
independent variables. Each model can then be tested regarding the following
linear model assumptions are checked:
\itemize{
\item tests performed by \code{\link[gvlma]{gvlma}}
\item the variance inflation factor is below 5 (models with two or more
independent variables)
\item no partial p.value of variables in the model is below 0.05
}
The model with the highest adjusted-R2 among those fulfilling the required
conditions is selected. A leave-one-out cross validation (LOO CV) is performed
by fitting the model coefficients using all observations except one and applying
the resulting model to predict the value for the remaining observation. In
case a transformation was performed beforehand, a bias correction is applied.
LOO CV statistics are then computed.
}
\examples{
data(quatre_montagnes)
# build ABA model for basal area, with all metrics as predictors
model_aba <- aba_build_model(quatre_montagnes$G_m2_ha, quatre_montagnes[, 9:76],
  transform = "boxcox", nmax = 3
)
# summary of regression model
summary(model_aba$model)
# validation statistics
model_aba$stats
# observed and predicted values
summary(model_aba$values)

# plot field values VS predictions in cross-validation
aba_plot(model_aba, main = "Basal area")
}
\seealso{
\code{\link{aba_combine_strata}} for combining models calibrated
on different strata, \code{\link{aba_plot}} for plotting model
cross-validation results, \code{\link[leaps]{regsubsets}} for variable selection,
 \code{\link{lma_check}} for linear model assumptions check,
 \code{\link{boxcox_itr_bias_cor}} for reverse Box-Cox transformation with bias
  correction.
}
