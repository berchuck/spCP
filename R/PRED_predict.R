###Function to get model fit diagnostics given a STBDwDM object
#'
#' predict.spCP
#'
#' Predicts future observations from the \code{\link{spCP}} model.
#'
#' @param object a \code{\link{spCP}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.spCP} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged observed outcomes along with the
#'  observational level parameters (\code{mu}, \code{tau}, and \code{alpha}).
#'
#' @return \code{predict.spCP} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{MuTauAlpha}}{A \code{list} containing three matrices, \code{mu},
#'   \code{tau} and \code{alpha}. Each matrix is dimension \code{NKeep x s}, where
#'   \code{s} is the number of new time points. Each matrix contains posterior
#'   samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{s} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x s}, where \code{s}
#'   is the number of new time points. Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spCP function
predict.spCP <- function(object, NewTimes, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.spCP(object)) stop('"object" must be of class spCP')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M

  ###Update DatObj
  DatObj$NewTimes <- NewTimes
  DatObj$NNewTimes <- length(NewTimes)

  ###Set mcmc object
  NKeep <- dim(object$delta)[1]

  ###Create parameter object
  Para <- list()
  Para$Beta0 <- object$beta0
  Para$Beta1 <- object$beta1
  Para$Lambda0 <- object$lambda0
  Para$Lambda1 <- object$lambda1
  Para$Eta <- object$eta

  ###Obtain samples of mu, tau and alpha using Bayesian krigging
  YFuture <- PredictFuture(DatObj, Para, NKeep)

  ###Format prediction samples for output
  Out <- list()
  for (i in 1:DatObj$NNewTimes) {
    Out[[i]] <- YFuture[ , , i]
    colnames(Out[[i]]) <- 1:M
    rownames(Out[[i]]) <- 1:NKeep
  }
  names(Out) <- NewTimes

  ###Return formated samples
  return(Out)

}



###Function to get model fit diagnostics given a CP object
#'
#' predict.CP
#'
#' Predicts future observations from the \code{\link{spCP}} model.
#'
#' @param object a \code{\link{spCP}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.spCP} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged observed outcomes along with the
#'  observational level parameters (\code{mu}, \code{tau}, and \code{alpha}).
#'
#' @return \code{predict.spCP} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{MuTauAlpha}}{A \code{list} containing three matrices, \code{mu},
#'   \code{tau} and \code{alpha}. Each matrix is dimension \code{NKeep x s}, where
#'   \code{s} is the number of new time points. Each matrix contains posterior
#'   samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{s} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x s}, where \code{s}
#'   is the number of new time points. Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spBDwDM function
predict.CP <- function(object, NewTimes, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.CP(object)) stop('"object" must be of class CP')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M

  ###Update DatObj
  DatObj$NewTimes <- NewTimes
  DatObj$NNewTimes <- length(NewTimes)

  ###Set mcmc object
  NKeep <- dim(object$delta)[1]

  ###Create parameter object
  Para <- list()
  Para$Beta0 <- object$beta0
  Para$Beta1 <- object$beta1
  Para$Lambda0 <- object$lambda0
  Para$Lambda1 <- object$lambda1
  Para$Eta <- object$eta

  ###Obtain samples of mu, tau and alpha using Bayesian krigging
  YFuture <- PredictFuture(DatObj, Para, NKeep)

  ###Format prediction samples for output
  Out <- list()
  for (i in 1:DatObj$NNewTimes) {
    Out[[i]] <- YFuture[ , , i]
    colnames(Out[[i]]) <- 1:M
    rownames(Out[[i]]) <- 1:NKeep
  }
  names(Out) <- NewTimes

  ###Return formated samples
  return(Out)

}



###Function to get model fit diagnostics given a spCP_novar object
#'
#' predict.spCP_novar
#'
#' Predicts future observations from the \code{\link{spCP}} model.
#'
#' @param object a \code{\link{spCP}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.STBDwDM} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged observed outcomes along with the
#'  observational level parameters (\code{mu}, \code{tau}, and \code{alpha}).
#'
#' @return \code{predict.STBDwDM} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{MuTauAlpha}}{A \code{list} containing three matrices, \code{mu},
#'   \code{tau} and \code{alpha}. Each matrix is dimension \code{NKeep x s}, where
#'   \code{s} is the number of new time points. Each matrix contains posterior
#'   samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{s} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x s}, where \code{s}
#'   is the number of new time points. Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spBDwDM function
predict.spCP_novar <- function(object, NewTimes, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.spCP_novar(object)) stop('"object" must be of class spCP_novar')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M

  ###Update DatObj
  DatObj$NewTimes <- NewTimes
  DatObj$NNewTimes <- length(NewTimes)

  ###Set mcmc object
  NKeep <- dim(object$delta)[1]

  ###Create parameter object
  Para <- list()
  Para$Beta0 <- object$beta0
  Para$Beta1 <- object$beta1
  Para$Lambda0 <- object$lambda
  Para$Lambda1 <- object$lambda
  Para$Eta <- object$eta

  ###Obtain samples of mu, tau and alpha using Bayesian krigging
  YFuture <- PredictFuture_novar(DatObj, Para, NKeep)

  ###Format prediction samples for output
  Out <- list()
  for (i in 1:DatObj$NNewTimes) {
    Out[[i]] <- YFuture[ , , i]
    colnames(Out[[i]]) <- 1:M
    rownames(Out[[i]]) <- 1:NKeep
  }
  names(Out) <- NewTimes

  ###Return formated samples
  return(Out)

}



###Function to get model fit diagnostics given a spCP_lmc object
#'
#' predict.spCP_lmc
#'
#' Predicts future observations from the \code{\link{spCP}} model.
#'
#' @param object a \code{\link{spCP}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.STBDwDM} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged observed outcomes along with the
#'  observational level parameters (\code{mu}, \code{tau}, and \code{alpha}).
#'
#' @return \code{predict.STBDwDM} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{MuTauAlpha}}{A \code{list} containing three matrices, \code{mu},
#'   \code{tau} and \code{alpha}. Each matrix is dimension \code{NKeep x s}, where
#'   \code{s} is the number of new time points. Each matrix contains posterior
#'   samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{s} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x s}, where \code{s}
#'   is the number of new time points. Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spBDwDM function
predict.spCP_lmc <- function(object, NewTimes, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.spCP_lmc(object)) stop('"object" must be of class spCP_lmc')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M

  ###Update DatObj
  DatObj$NewTimes <- NewTimes
  DatObj$NNewTimes <- length(NewTimes)

  ###Set mcmc object
  NKeep <- dim(object$delta)[1]

  ###Create parameter object
  Para <- list()
  Para$Beta0 <- object$beta0
  Para$Beta1 <- object$beta1
  Para$Lambda0 <- object$lambda0
  Para$Lambda1 <- object$lambda1
  Para$Eta <- object$eta

  ###Obtain samples of mu, tau and alpha using Bayesian krigging
  YFuture <- PredictFuture(DatObj, Para, NKeep)

  ###Format prediction samples for output
  Out <- list()
  for (i in 1:DatObj$NNewTimes) {
    Out[[i]] <- YFuture[ , , i]
    colnames(Out[[i]]) <- 1:M
    rownames(Out[[i]]) <- 1:NKeep
  }
  names(Out) <- NewTimes

  ###Return formated samples
  return(Out)

}


