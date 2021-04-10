# package lidaRtRee
# Copyright INRAE
# Author(s): Jean-Matthieu Monnet
# Licence: GPL-3
################################ 
#' Calibrates and validates area-based models
#' 
#' The function can first apply a Box-Cox transformation to the dependent variable, in order to normalize its distribution, or a log transformation to the whole dataset. Then it uses \code{\link[leaps]{regsubsets}} to find the 20 linear regressions with the best adjusted-R2 among combinations of at most \code{nmax} independent variables. Each model can then be tested regarding the following linear model assumptions are checked:
#' \itemize{
#' \item tests performed by \code{\link[gvlma]{gvlma}}
#' \item the variance inflation factor is below 5 (models with two or more independent variables)
#' \item no partial p.value of variables in the model is below 0.05
#' }
#' The model with the highest adjusted-R2 among those fulfilling the required conditions is selected. A leave-one-out cross validation (LOO CV) is performed by fitting the model coefficients using all observations except one and applying the resulting model to predict the value for the remaining observation. In case a transformation was performed beforehand, a bias correction is applied. LOO CV statistics are then computed.
#' 
#' @param variable vector. dependent variable values
#' @param predictors data.frame. independent variables (columns: metrics, lines: observations). Row names are used for the output predicted values
#' @param transform string. transformation to be applied to data (\code{"none"}, \code{"boxcox"}: Box-Cox transformation applied only to the dependent variable, \code{"log"}: log transformation applied to both dependent and independent variables)
#' @param nmax numeric. maximum number of independent variables in the model
#' @param test vector. which tests should be satisfied by the models, one to three in \code{"partial_p"}, \code{"vif"}, \code{"gvlma"}
#' @param xy data.frame or matrix of easting and northing coordinates of observations: not used in the function but exported in the result for use in further inference functions
#' @param threshold vector of length two. minimum and maximum values of threshold to apply to predicted values
#' @seealso \code{\link{ABAmodelCombineStrata}} for combining models calibrated on different strata, \code{\link{ABAmodelPlot}} for plotting model cross-validation results, \code{\link[leaps]{regsubsets}} for variable selection, \code{\link{lmaCheck}} for linear model assumptions check, \code{\link{iBoxcoxTrBiasCor}} for reverse Box-Cox transformation with bias correction.
#' @examples
#' # load Quatre Montagnes dataset
#' data(quatremontagnes)
#' # build ABA model for basal area, with all metrics as predictors
#' model.aba <- ABAmodel(quatremontagnes$G.m2.ha, quatremontagnes[, 9:76],
#' transform = "boxcox", nmax = 3)
#' # summary of regression model
#' summary(model.aba$model)
#' # validation statistics
#' model.aba$stats
#' # observed and predicted values
#' summary(model.aba$values)
#' 
#' # plot field values VS predictions in cross-validation
#' ABAmodelPlot(model.aba, main = "Basal area")
#' @return a list with three elements
#' \itemize{
#' \item \code{model}: list with one regression model (output from \code{\link[stats]{lm}}),
#' \item \code{stats}: model statistics (root mean square error estimated in leave-one-out cross validation, coefficient of variation of rmse, p-value of wilcoxon test of observed and predicted values, p-value of t-test of observed and predicted values, p-value of anova of observed and predicted values, correlation of observed and predicted values, R2 of observed and predicted values, variance of regression residuals)
#' \item \code{values}: data.frame with observed and values predicted in cross-validation.
#' }
#' @export
#'
ABAmodel <- function(variable, predictors, transform="none", nmax=3, test=c("partial_p", "vif", "gvlma"), xy=NULL, threshold=NULL)
{
  ##################
  # build data.frame with dependent variable and predictors
  df <- data.frame(dep.var=variable, predictors)
  row.names(df) <- row.names(predictors)
  # remove lines where dependent variable is negative or missing
  dummy <- which(!is.na(df$dep.var) & df$dep.var>0)
  if (length(dummy)<nrow(df)) { print("Negative or NA observations of dependent variable have been removed") }
  df <- df[dummy,]
  #
  # column to remove from data.frame to leave only independent variables
  var.out.indices <- 1
  #
  ##################
  # variable transformation
  # no transformation
  if (transform=="none")
  {
    df.transform <- df
    lambda <- NA
  }
  # Box-Cox transformation of the dependent variable
  if (transform=="boxcox")
  {
    # selecting the best Box-Cox parameter for normalization of the dependent variable distribution
    lambda <- car::powerTransform(df[ , 1])$lambda
    # create a copy of data
    df.transform <- df
    # apply Box-Cox transformation to dependent variable
    df.transform[,1] <- BoxcoxTr(df.transform[,1], lambda)
  }
  # log transformation of all data
  if (transform=="log")
  {
    df.transform <- log(df)
    lambda <- NA
    # check that no non-finite values have been created
    test.finite <- apply(df.transform, 2, function(x){all(is.finite(x))})
    # remove variable with non-finite values if present
    if (!all(test.finite))
    {
      warning(paste("Removed variables because of non-finite values after log transform: ", paste(names(df.transform[which(!test.finite)]))))
    }
    # add variables with non-finite values to vector of columns to be removed before model building
    var.out.indices <- c(var.out.indices, which(!test.finite))
  }
  #
  ##################
  # selection of 20 combinations of variables with best adj-R2
  var.combi <-
    leaps::regsubsets(
      y = df.transform[, 1],
      x = df.transform[, -var.out.indices],
      nbest = 20,
      nvmax = nmax - 1,
      method = "exhaustive",
      really.big = TRUE
    )
  ##################
  # check linear model assumptions
  # matrix of booleans indicating formulas
  fomula.bool <- summary(var.combi)$which
  # check all formulas for modeling hypotheses
  test_lin <- apply(fomula.bool, 1, function(x)
  {
    formule <-
      stats::as.formula(paste("dep.var~", paste(colnames(fomula.bool)[which(x == TRUE)[-1]], collapse ="+")))
    lmaCheck(formule,
             df.transform,
             max.pvalue = 0.05,
             max.vif = 5)
  })
  # bind results
  test_lin <- do.call(rbind, test_lin)
  # removal of models which do not satisfy hypotheses
  if (length(test)>0)
  {
    test_lin <- test_lin[which(apply(test_lin[,test], 1, all)), ]
  }
  ##################
  # output the best remaining model
  if (dim(test_lin)[1]>0)
  {
    # keep only best adjR2
    ir <- which(test_lin$adjR2==max(test_lin$adjR2))[1]
    nbobs <- nrow(df)
    formule <- as.character(test_lin$formula[ir])
    # create result data.frame with model information
    modeles <-
      data.frame(
        n = nbobs,
        formula = formule,
        adjR2 = test_lin$adjR2[ir],
        transform = transform,
        lambda = lambda
      )
    # build formula
    formule <- stats::as.formula(paste0("dep.var ~ ", formule))
    # save global regression model
    regr <- eval(parse(text=paste0("stats::lm(", formule[2], "~", formule[3], ",data = df.transform)")))
    ##################
    # leave one out cross validation
    # create vectors of outputs (predictions, residuals variance)
    prediction <- var.res <- rep(NA, nbobs)
    # Leave-One-Out Cross Validation
    for (i in 1:nbobs)
    {
      # fit lm with formula and n-1 values
      dummy <- stats::lm(formula = formule, data=df.transform[-i,])
      # predict on remaining value
      prediction[i] <- stats::predict.lm(dummy,df.transform[i,])
      # record residuals variance for model of prediction i
      var.res[i] <- sum(dummy$residuals^2)/dummy$df.residual
    }
    ##################
    # apply back-transformation if required
    if (transform=="boxcox")
    {
      # bias correction using the variance of residuals [Tiefelsdorf (January 2013): A Note on the Reverse Box-Cox Transformation]
      prediction <- iBoxcoxTrBiasCor(prediction, lambda, var.res)
    }
    if (transform=="log")
    {
      # bias correction as in [Bouvier et al.: Generalizing...]
      prediction <- exp(prediction) * exp(var.res/2)
    }
    #
    # threshold for predicted values
    if (!is.null(threshold))
    {
      prediction <- pmin(prediction, threshold[2])
      prediction <- pmax(prediction, threshold[1])
    }
    #
    ##################
    # validation statistics
    # vector of dependent variable
    dummy <- df[,1]
    # compute root mean square error
    modeles$rmse <- (sum((prediction-df[,1])^2)/nbobs)^0.5
    # compute coefficient of variation of rmse
    modeles$cvrmse <- modeles$rmse/mean(df[,1])
    # compute p-value of wilcoxon test of observed and predicted values
    modeles$pwil <- stats::wilcox.test(prediction,df[,1],paired=TRUE)$p.value
    # compute p-value of t-test of observed and predicted values
    modeles$pttest <- stats::t.test(prediction,df[,1],paired=TRUE)$p.value
    # compute p-value of anova of observed and predicted values
    modeles$paov <- stats::anova(stats::lm(c(prediction,df[,1])~as.factor(c(rep(0,nbobs),rep(1,nbobs)))))[5][1,1]
    # compute correlation of observed and predicted values
    modeles$cor <- stats::cor(df[,1],prediction,use="pairwise.complete.obs")
    # compute R2 of observed and predicted values
    modeles$looR2 <- 1-sum((prediction-df[,1])^2)/sum((df[,1]-mean(df[,1]))^2)
    # compute variance of regression residuals
    modeles$var.res <- sum(regr$residuals^2)/regr$df.residual
    # row.names
    # rename first column in case data.frame
    row.names(modeles)[1] <- "dep.var"
    #
    # prepare output list
    output <- list("model"=regr, "stats"=modeles, "values"=data.frame(field=df[,1], predicted=prediction, residual=df[,1]-prediction, row.names=row.names(df)))
    # add observation coordinates if provided
    if (!is.null(xy)) 
    {
      output$values$x <- xy[,1]
      output$values$y <- xy[,2]
    }
    return(output)
  } else {return(NULL)}
}

#########################
#' Checks linear model assumptions of a multiple regression model
#' 
#' The performed tests are:
#' \itemize{
#' \item partial p.values calculated by \code{\link[stats]{lm}} are all below a given value
#' \item tests implemented by \code{\link[gvlma]{gvlma}}
#' \item variance inflation factors calculated by \code{\link[car]{vif}} are all below a given value
#' }
#'
#'@param formule formula. model to be evaluated
#'@param df data.frame. data to evaluate the model
#'@param max.pvalue numeric. maximum p-value of variables included in the model
#'@param max.vif numeric. maximum variance inflation factor of variables included in the model
#'@return a one line data.frame with 5 columns.
#'\itemize{
#' \item a string: evaluated formula
#' \item a numeric: the adjusted R squared of the model
#' \item a boolean: do all variables in the model have a partial p-value < \code{max.pvalue}
#' \item a boolean: are all tests implemented by \code{\link[gvlma]{gvlma}} false
#' \item a boolean: is the variance inflation factor computed with \code{\link[car]{vif}} of all variables < \code{max.vif} 
#' }
#' @examples
#' # load Quatre Montagnes dataset
#' data(quatremontagnes)
#' # fit lm model
#' model <- lm(G.m2.ha ~ zmax + zq95, data = quatremontagnes)
#' lmaCheck(eval(model$call[[2]]), quatremontagnes)
#' # trying with Box-Cox transformation of dependent variable  
#' # and other independent variables
#' model <- lm(BoxcoxTr(G.m2.ha, -0.14) ~ Tree.meanH + Tree.density + zpcum7 , data = quatremontagnes)
#' lmaCheck(eval(model$call[[2]]), quatremontagnes)
#'@export
#'
lmaCheck <- function(formule,
                     df,
                     max.pvalue = 0.05,
                     max.vif = 5)
{
  # build multiple regression model
  reg <- stats::lm(formule, data = df)
  # check linear model assumptions with gvlma
  dummy <- gvlma::gvlma(reg, df, alphalevel = 0.1)
  # select tests results
  testlma <- c(
    dummy$GlobalTest$GlobalStat4$Decision,
    dummy$GlobalTest$DirectionalStat1$Decision,
    dummy$GlobalTest$DirectionalStat2$Decision,
    dummy$GlobalTest$DirectionalStat3$Decision,
    dummy$GlobalTest$DirectionalStat4$Decision
  )
  # invert test results
  testlma <- abs(testlma - 1)
  # compute variance inflation factor if more than one variable in the model
  if (length(strsplit(as.character(formule)[3], "\\+")[[1]]) > 1)
  {
    testvif <- max(car::vif(reg))
  } else {
    testvif <- 0
  }
  # prepare regression model summary
  reg <- summary(reg)
  # output result data.frame with test results
  data.frame(
    "formula" = as.character(formule)[[3]],
    "adjR2" = reg$adj.r.squared,
    "partial_p" = (max(reg$coefficients[-1, 4]) < max.pvalue),
    "gvlma" = all(as.logical(testlma)),
    "vif" = testvif < max.vif
  )
}
#########################
#' Box-Cox Transformation
#'
#'@param x vector or RasterLayer. values to be transformed
#'@param lambda numeric. parameter of Box-Cox transformation
#'@return a vector or RasterLayer of transformed values
#'@seealso \code{\link{iBoxcoxTr}} inverse Box-Cox transformation, \code{\link{iBoxcoxTrBiasCor}} inverse Box-Cox transformation with bias correction.
#' @examples
#' x <- 1:10
#' BoxcoxTr(x, -2)
#' BoxcoxTr(x, 0)
#' BoxcoxTr(x, 0.5)
#' BoxcoxTr(x, 2)
#' 
#' # plot functions
#' curve(BoxcoxTr(x, 1.5), 1, 5, main = "Box Cox transform", xlab = "x",
#' ylab = "Boxcox(x, lambda)", col = "red")
#' curve(BoxcoxTr(x, -2), 1, 5, col = "green", add = TRUE)
#' curve(BoxcoxTr(x, 0), 1, 5, col = "blue", add = TRUE)
#' curve(BoxcoxTr(x, 0.5), 1, 5,  col = "black", add = TRUE)
#' curve(BoxcoxTr(x, 1), 1, 5, col = "pink", add = TRUE)
#' legend("topleft", legend = rev(c(-2, 0, 0.5, 1, 1.5, "lambda")),
#' col = rev(c("green", "blue", "black", "pink", "red", NA)), lty = 1)
#'@export
BoxcoxTr <- function(x,lambda)
{
  count.neg <- ifelse(class(x)[1]=="RasterLayer", sum(raster::values(x<0), na.rm=TRUE), sum(x<0, na.rm=TRUE))
  if (count.neg >0)
  {
    x[x<0] <- NA
    warning(paste0(count.neg, " negative value(s) set to NA"))
  }
  if (lambda !=0)
  {
    (x^lambda-1)/lambda
  } else {
    log(x)
  }
}
########################
#' Inverse Box-Cox transformation
#'
#'@param x vector or RasterLayer. values to be transformed
#'@param lambda numeric. parameter of Box-Cox transformation
#'@return a vector or RasterLayer of transformed values
#'@seealso \code{\link{BoxcoxTr}} Box-Cox transformation, \code{\link{iBoxcoxTrBiasCor}} inverse Box-Cox transformation with bias correction.
#' @examples
#' x <- 1:10
#' iBoxcoxTr(x, 0)
#' iBoxcoxTr(x, 0.5)
#' iBoxcoxTr(x, 2)
#' iBoxcoxTr(BoxcoxTr(x,2), 2)
#' 
#' # plot functions
#' curve(iBoxcoxTr(x, 0), 0, 3, col = "blue", main = "inverse Box Cox transf.",
#' xlab = "x", ylab = "inverse Boxcox(x, lambda)")
#' curve(iBoxcoxTr(x, 1.5), 0, 3, col = "red", add = TRUE)
#' curve(iBoxcoxTr(x, 0.5), 0, 3,  col = "black", add = TRUE)
#' curve(iBoxcoxTr(x, 1), 0, 3, col = "pink", add = TRUE)
#' legend("topleft", legend = c("lambda", 0, 0.5, 1, 1.5),
#' col = c(NA, "blue", "black", "pink", "red"), lty = 1)
#'@export
iBoxcoxTr <- function(x,lambda)
{
  count.neg <- ifelse(class(x)[1]=="RasterLayer", sum(raster::values(x<0), na.rm=TRUE), sum(x<0, na.rm=TRUE))
  if (count.neg >0)
  {
    x[x<0] <- NA
    warning(paste0(count.neg, " negative value(s) set to NA"))
  }
  if (lambda !=0)
  {
    (lambda * x +1)^(1/lambda)
  } else {
    exp(x)
  }
}
#########################
#' Inverse Box-Cox transformation with bias correction
#' 
#' Inverse Box-Cox transform with bias correction as suggested by Pu & Tiefelsdorf (2015). Here `varmod` is not the local prediction variance as suggested in the paper but the model residuals variance. For variance computation, use `n-p` instead of `n-1`, with `p` the number of variables in the model.
#' 
#'@param x vector or RasterLayer. values to be transformed
#'@param lambda numeric. parameter of Box-Cox transformation
#'@param varmod numeric. model residuals variance
#'@references Xiaojun Pu and Michael Tiefelsdorf, 2015. A variance-stabilizing transformation to mitigate biased variogram estimation in heterogeneous surfaces with clustered samples. \doi{10.1007/978-3-319-22786-3_24}
#'@return a vector or RasterLayer
#'@seealso \code{\link{BoxcoxTr}} Box-Cox transformation, \code{\link{iBoxcoxTr}} inverse Box-Cox transformation.
#'@examples
#' x <- 1:10
#' iBoxcoxTr(x, 0.3)
#' iBoxcoxTrBiasCor(x, 0.3, 0)
#' iBoxcoxTrBiasCor(x, 0.3, 2)
#' 
#' # plot functions
#' curve(iBoxcoxTr(x, 0.3), 0, 3, col = "blue",
#' main = "inverse Box Cox transf., lambda = 0.3",
#' xlab = "x", ylab = "inverse Boxcox(x, lambda = 0.3)")
#' curve(iBoxcoxTrBiasCor(x, 0.3, 1), 0, 3, col = "red", add = TRUE)
#' curve(iBoxcoxTrBiasCor(x, 0.3, 2), 0, 3,  col = "black", add = TRUE)
#' legend("topleft", legend = c("residuals variance  = 2",
#' "residuals variance  = 1", "residuals variance not accounted for"),
#' col = c("black", "red", "blue"), lty = 1)
#'@export
iBoxcoxTrBiasCor <- function(x, lambda, varmod)
{
  if (lambda !=0)
  {
    iBoxcoxTr(x,lambda) * (1+(varmod*(1-lambda)/(2*(lambda*x+1)^2)))
  } else {
    exp(x) * exp(varmod/2)
  }
}
####################################
#' Combines a list of ABA models into a single ABA model object
#' 
#' Combines a list of models (obtained with \code{\link{ABAmodel}}) into a single object. Typically used to merge stratum-specific models into one object. Validation statistics are computed for the combined strata, making it easier to compare prediction performance with an unstratified model.
#'
#' @param model.list list. stratum-specific models returned by \code{\link{ABAmodel}}
#' @param plotsId vector. "plotsId" for ordering row names in the "values" element of the output list
#' @return a list with three elements
#' \itemize{
#' \item \code{model}: a list of regression models corresponding to each stratum (output from \code{\link[stats]{lm}}),
#' \item \code{stats}:model statistics of each stratum-specific model (as in \code{\link{ABAmodel}}) plus one line corresponding to statistics for all strata (COMBINED)
#' \item \code{values}: data.frame with observed and values predicted in cross-validation, and information on which stratum it belongs to.
#' }
#' @seealso \code{\link{ABAmodel}} for calibrated ABA model, \code{\link{ABAmodelPlot}} for plotting model cross-validation results.
#' @examples
#' # load Quatre Montagnes dataset
#' data(quatremontagnes)
#' # initialize list of models
#' model.ABA.stratified <- list()
#' # calibrate basal area prediction model for each stratum
#' for (i in levels(quatremontagnes$stratum))
#' {
#'   subsample <- which(quatremontagnes$stratum==i)
#'   model.ABA.stratified[[i]] <-
#'   ABAmodel(quatremontagnes[subsample, "G.m2.ha"],
#'   quatremontagnes[subsample, 9:76], transform="boxcox", nmax=4,
#'   xy = quatremontagnes[subsample,c("X", "Y")])
#' }
#' # combine models in single object
#' model.ABA.stratified <- ABAmodelCombineStrata(model.ABA.stratified,
#' quatremontagnes$plotId)
#' # display content of output list
#' model.ABA.stratified$model
#' model.ABA.stratified$stats
#' summary(model.ABA.stratified$values)
#' 
#' # plot field values VS predictions in cross-validation
#' ABAmodelPlot(model.ABA.stratified)
#' @export
#'
ABAmodelCombineStrata <- function(model.list, plotsId=NULL)
{
  model.combined <- list()
  # extract models into separate list
  model.combined$model <- lapply(model.list, function(x) {x$model})
  # bind statistics into single data.frame
  model.combined$stats <- do.call(rbind, lapply(model.list, function(x) {x$stats}))
  # bind predicted values and add stratum
  model.combined$values <- do.call(rbind, lapply(as.list(names(model.list)), function(x){data.frame(model.list[[x]]$values, stratum=x)}))
  # 
  # compute rmse of combined models
  rmse <- (sum((model.combined$values$residual)^2)/nrow(model.combined$values))^0.5
  # add statistics
  model.combined$stats[nrow(model.combined$stats)+1,] <- c(nrow(model.combined$values),
                                           NA,
                                           NA,
                                           NA,
                                           NA,
                                           rmse,
                                           rmse/mean(model.combined$values$field),
                                           stats::wilcox.test(model.combined$values$residual)$p.value,
                                           stats::t.test(model.combined$values$residual)$p.value,
                                           NA,
                                           stats::cor(model.combined$values$field,model.combined$values$predicted, use="pairwise.complete.obs"), 
                                           1-sum((model.combined$values$residual)^2)/sum((model.combined$values$field-mean(model.combined$values$field))^2),
                                           NA)
  # update names of models 
  row.names(model.combined$stats) <- c(names(model.list),"COMBINED")
  names(model.combined$model) <- names(model.list)
  # order values same as original data
  if (!is.null(plotsId))
  { model.combined$values <- model.combined$values[plotsId,]}
  #
  # convert stratum field to factor
  model.combined$values$stratum <- as.factor(model.combined$values$stratum)
  model.combined
}
####################################
#' Plots observed VS values predicted in leave one out cross validation of an \code{\link{ABAmodel}}
#'
#' @param modell list. as returned by \code{\link{ABAmodel}}
#' @param disp.text boolean. indicates if points should be labeled with id
#' @param col color to be passed to \code{\link[graphics]{plot}}, default is black for single models, depends on stratum in stratified models
#' @param ... other parameters to be passed to \code{\link[graphics]{plot}}, \code{xlab} and \code{ylab} are automatically setup
#' @examples
#' # load Quatre Montagnes dataset
#' data(quatremontagnes)
#' # build ABA model for basal area, with all metrics as predictors
#' model.aba <- ABAmodel(quatremontagnes$G.m2.ha, quatremontagnes[, 9:76],
#' transform = "boxcox", nmax = 3)
#' 
#' # plot field values VS predictions in cross-validation
#' ABAmodelPlot(model.aba, main = "Basal area")
#' @return nothing
#' @export
#'
ABAmodelPlot <- function(modell, disp.text=F, col = NULL, ...)
{
  # color
  if (is.null(col))
  {
    if (nrow(modell$stats)>1) # if stratified model, color values by stratum
    {
      col <- as.numeric(modell$values$stratum)
    } else { col <- "black"}
  }
  #
  # display points, symbol color is white in case text is displayed
  graphics::plot(modell[["values"]]$field, modell[["values"]]$predicted, asp=1,
                 xlab="Field", ylab="Predicted in LOOCV", col=ifelse(rep(disp.text, length(col)), "white", col), ...)
  # add 1:1 line
  graphics::abline(c(0,1))
  # add legend if stratified model
  if (nrow(modell$stats)>1)
  {
    graphics::legend("topleft", paste0(levels(modell$values$stratum), " model"), fill=1:length(levels(modell$values$stratum)))
  }
  # display text if required
  if (disp.text)
  {
    graphics::text(modell[["values"]]$field, modell[["values"]]$predicted, labels=row.names(modell[["values"]]), cex=0.8, col = col)
  }
}
####################################
#' Mapping of ABA prediction models
#' 
#' Applies calibrated area-based prediction models output of \code{\link{ABAmodel}} to a RasterStack of metrics to obtain a raster of predictions
#'
#' @param aba.model model returned by \code{\link{ABAmodel}} or \code{\link{ABAmodelCombineStrata}}
#' @param metrics.map RasterStack. metrics returned e.g by \code{\link[lidR]{grid_metrics}}
#' @param stratum string. indicates which layer of metrics.map contains the \code{stratum} in \code{modell}
#' @param addError boolean. indicates whether errors sampled from a normal distribution N(0, sigma(residuals)) should be added to fitted values; implemented only for \code{log} transformation case
#' @examples 
#' # load data
#' data(quatremontagnes)
#' # build model
#' ABA.model <- ABAmodel(quatremontagnes$G.m2.ha, quatremontagnes[,9:76],
#' transform = "boxcox")
#' # build example raster to apply model
#' quatremontagnes$X <- rep(1:8, 12)
#' quatremontagnes$Y <- rep(1:12, each = 8)
#' metrics.map <- raster::rasterFromXYZ(quatremontagnes[, c(2,3,9:76)])
#' predict.map <- lidaRtRee::ABApredict(ABA.model, metrics.map)
#' 
#' # plot map
#' raster::plot(predict.map, main = "predictions")
#' @seealso \link{ABAmodel} for model fitting and \link{ABAmodelCombineStrata} for combining stratified models, \link{cleanRaster} for applying spatial mask and value thresholds to a raster.
#' @return a raster of predictions obtained by applying the model \code{aba.model} to the observations in \code{metrics.map}
#' @export
#'
ABApredict <- function(aba.model, metrics.map, stratum=NULL, addError=FALSE)
{
  # create factor of stratum if not existing
  if (is.null(stratum))
  {
    metrics.map$stratum <- factor("all")
    metrics.map$stratum@data@attributes[[1]] <- data.frame(ID=1, stratum="all")
    row.names(aba.model$stats)[1] <- "all"
    aba.model$model <- list("all"=aba.model$model)
  }
  #
  dummy <- list()
  # loop on strata
  for (i in 1:nrow(metrics.map[["stratum"]]@data@attributes[[1]]))
  {
    # extract stratum
    stratum <- as.character(metrics.map[["stratum"]]@data@attributes[[1]][i, "stratum"])
    stratum.ID <- metrics.map[["stratum"]]@data@attributes[[1]][i, "ID"]
    # predict on all cells
    if (aba.model$stats[stratum,"transform"] == "boxcox") # case of Box-Cox transform
    {
      # apply linear model
      variables <- names(aba.model$model[[stratum]]$coefficients)
      variables <- variables[variables !="(Intercept)"]
      dummy[[stratum]] <- raster::predict(metrics.map[[variables]], aba.model$model[[stratum]])
      # back-transform
      dummy[[stratum]] <- iBoxcoxTrBiasCor(dummy[[stratum]], aba.model$stats[stratum,"lambda"], aba.model$stats[stratum,"var.res"])
      if (addError==TRUE) {warning("Error sampling not implemented in the boxcox transformation case")}
    } 
    if (aba.model$stats[stratum,"transform"] == "log") # case of case of log-log transform
    {
      variables <- names(aba.model$model[[stratum]]$coefficients)
      variables <- variables[variables !="(Intercept)"]
      newdata <- log(metrics.map[[variables]])
      names(newdata) <- variables
      dummy[[stratum]] <- raster::predict(newdata, aba.model$model[[stratum]])
      if (addError==TRUE) # sampling of errors
      {
        rastResidual <- dummy[[stratum]]
        raster::values(rastResidual) <- stats::rnorm(length(rastResidual), 0, sqrt(aba.model$stats[stratum,"var.res"]))
        dummy[[stratum]]<- exp(dummy[[stratum]] + rastResidual)
      } else { # bias correction in the case of log transformation
        dummy[[stratum]]<- exp(dummy[[stratum]]) * exp(aba.model$stats[stratum,"var.res"]/2)
      }
    }
    if (aba.model$stats[stratum,"transform"] == "none") # case of case of no transform
    {
      variables <- names(aba.model$model[[stratum]]$coefficients)
      variables <- variables[variables !="(Intercept)"]
      dummy[[stratum]] <- raster::predict(metrics.map[[variables]], aba.model$model[[stratum]])
      if (addError==TRUE) {warning("Error sampling not implemented in the none transformation case")}
    }
    # set prediction outside of strata to NA
    dummy[[stratum]][metrics.map$stratum != stratum.ID] <- NA
  }
  # merge strata results
  return(Reduce(function(...)raster::merge(..., tolerance=1), dummy))
}
####################################
#' Applies thresholds and mask to a rasterLayer object
#' 
#' Applies a lower and upper thresholds to the values of the input raster. If the mask input is provided, first all NA values in the raster are set to 0, then the raster in multiplied by the mask. Cells to be masked should therefore have a NA value in the mask raster object.
#' 
#' @param rast raster object. 
#' @param minmax vector of two numeric values. minimum and maximum thresholds to apply to `rast` values
#' @param mask raster object. mask to be applied (multiplication with input raster `rast`)
#' @examples
#' # load data
#' data(quatremontagnes)
#' # build model
#' ABA.model <- lidaRtRee::ABAmodel(quatremontagnes$G.m2.ha,
#' quatremontagnes[, 9:76], transform = "boxcox")
#' # build example raster to apply model
#' quatremontagnes$X <- rep(1:8, 12)
#' quatremontagnes$Y <- rep(1:12, each = 8)
#' metrics.map <- raster::rasterFromXYZ(quatremontagnes[, c(2,3,9:76)])
#' predict.map <- lidaRtRee::ABApredict(ABA.model, metrics.map)
#' # create raster mask
#' mask <- predict.map
#' # set values to 1 or NA
#' raster::values(mask) <- rep(c(1, 1, NA), each = 32)
#' # apply thresholds and mask
#' predict.map.clean <- cleanRaster(predict.map, c(40, 70), mask)
#' 
#' # plot maps
#' raster::plot(predict.map, main = "Predictions")
#' raster::plot(mask, main = "Mask", legend = FALSE)
#' raster::plot(predict.map.clean, main = "Cleaned predictions")
#' @return a raster object
#' @export
#'
cleanRaster <- function(rast, minmax = c(-Inf, +Inf), mask = NULL)
{
  # if mask is present
  if (! is.null(mask))
  {
    # fill NA values in rast with 0
    rast[is.na(rast)] <- 0
    # then apply mask
    rast <- rast * mask
  }
  # 
  rast[rast<=minmax[1]] <- minmax[1]
  rast[rast>=minmax[2]] <- minmax[2]
  #
  return(rast)
}
####################################
#' computes inference from area-based model and predicted values
#' @param modell a model returned by by \code{\link{ABAmodel}} or \code{\link{ABAmodelCombineStrata}}
#' @param r.predictions raster of predicted values 
#' @param type string vector specifying which estimators should be computed (one or several in "SRS", "ED", "D", "STR", "SYNT")
#' @param r.mask raster to mask region of interest (NA values), may contain post-stratification categories (should be integer, positive values)
#' @return a dataframe with estimation of parameter value and standard deviation of estimation for all required estimators.
#' @export
#'
ABAinference <- function(modell, r.predictions, type=c("SRS", "ED", "D", "STR", "SYNT"), r.mask=NULL)
{
  inference <- list()
  # extract observations
  observations <- modell$values
  # extract predicted values
  # apply mask if present
  if (!is.null(r.mask))
  {
    pixels <- r.predictions * (r.mask>=0)
    # apply mask also to observations
    coord <- modell$values[,c("x", "y")]
    sp::coordinates(coord) <- ~ x+y
    observations$mask <- raster::extract(r.mask, coord)
    observations <-  observations[!is.na(observations$mask),]
  } else {
    pixels <- r.predictions
  }
  # raster pixels may contain na values
  #
  # number of observations
  n <- nrow(observations)
  # number of predictions (no NA values in raster)
  N <- sum(!is.na(raster::values(pixels)))
  # number of parameters in the model
  # TAKE MAX NUMBER OF PARAMETERS IN CASE OF STRATIFIED MODELS
  n.para <- length(strsplit(as.character(modell$stats$formula[1]), "+", fixed=TRUE)[[1]])+1
  #
  # Simple random sampling (SRS)
  if (is.element("SRS", type))
  {
    inference[["SRS"]] <- data.frame(mean=mean(observations$field), var=stats::var(observations$field))
  }
  # 
  # Generalized difference estimator
  if (is.element("ED", type))
  {
    dummy.prediction.bias <- mean(observations$predicted-observations$field)
    dummy.mean <- mean(raster::values(pixels), na.rm=TRUE)-dummy.prediction.bias
    dummy.var <- sum(((observations$predicted-observations$field)-dummy.prediction.bias)^2)/(n-n.para)#/n
    inference[["ED"]] <- data.frame(mean=dummy.mean, var=dummy.var)
    # use n number of observations in ROI or in calibration dataset
  }
  #
  # Model-assisted difference estimator (D)
  # caution : no true correspondance between the spatial location of observations (calibration step) and predictions (mapping step)
  # here pixels containing an observation are replaced by the observation, but bias is still computed by substracting the observation and the prediction at the exact location
  if (is.element("D", type))
  {
    observed.pixels <- raster::cellFromXY(pixels, observations[,c("x","y")])
    dummy.prediction.bias <- mean(observations$predicted-observations$field)
    dummy.var <- sum(((observations$predicted-observations$field)-dummy.prediction.bias)^2)/(n-n.para)#/n (erreur)
    inference[["D"]] <- data.frame(
      mean=(sum(observations$field)+sum(raster::values(pixels)[-observed.pixels], na.rm=TRUE))/N - mean(observations$predicted-observations$field),
      var=dummy.var)
  }
  #
  # Stratified estimator (STR)
  if (is.element("STR", type) & !is.null(r.mask))
  {
    # check that there are pixels and observations in all categories
    if (setequal(stats::na.omit(unique(raster::values(r.mask))), stats::na.omit(unique(observations$mask))))
    {
      dummy <- list()
      # compute weights of each category
      surface <- table(raster::values(r.mask))
      surface <- surface/sum(surface)
      # for each category
      for (i in (names(surface)))
      {
        # selec observations
        temp.selec <- which(observations$mask==as.numeric(i))
        # compute mean and variance in category
        dummy[[i]] <- data.frame(W=surface[i], m=mean(observations$field[temp.selec]), n=length(temp.selec), var=stats::var(observations$field[temp.selec]))
      }
      dummy <- do.call(rbind,dummy)
      inference[["STR"]] <- data.frame(mean=sum(dummy$W*dummy$m), var=sum((dummy$var*dummy$W)^2/dummy$n))
    } else {warning("Impossible to compute STR inference: categories not present in both observations and pixels")}
  }
  # synthetic regression estimator (SYNT)
  if (is.element("SYNT", type))
  {
    inference[["SYNT"]] <- data.frame(mean=mean(raster::values(pixels), na.rm=TRUE), var=NA)
  }
  # convert list to data.frame
  inference <- do.call(rbind, inference)
  # compute sd
  inference$sd <- sqrt(inference$var)
  return(inference)
}
