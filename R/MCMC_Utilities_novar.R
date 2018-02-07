###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples_novar <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  PhiIndeces <- DatObj$PhiIndeces
  tNu <- DatObj$tNu

  ###Format raw samples
  RawSamples <- t(RawSamples)
  Alpha <- RawSamples[, 1, drop = FALSE]
  Delta <- RawSamples[, 2:5]
  Sigma <- RawSamples[, 6:15]
  Phi <- RawSamples[, 16:((4 * M) + 15)]
  Beta0 <- Phi[, PhiIndeces[1, ] + 1]
  Beta1 <- Phi[, PhiIndeces[2, ] + 1]
  Lambda <- Phi[, PhiIndeces[3, ] + 1]
  Eta <- Phi[, PhiIndeces[4, ] + 1]
  Theta <- apply(Eta, 2, function(x) pmin(tNu, exp(x)))
  colnames(Alpha) <- "Alpha"
  colnames(Delta) <- paste0("Delta", 1:4)
  colnames(Sigma) <- c(paste0("Sigma", 1:4, "1"), paste0("Sigma", 2:4, "2"), paste0("Sigma", 3:4, "3"), "Sigma44")
  colnames(Beta0) <- paste0("Beta0(", 1:M, ")")
  colnames(Beta1) <- paste0("Beta1(", 1:M, ")")
  colnames(Lambda) <- paste0("Lambda(", 1:M, ")")
  colnames(Eta) <- paste0("Eta(", 1:M, ")")
  colnames(Theta) <- paste0("Theta(", 1:M, ")")
  Out <- list(Alpha = Alpha, Delta = Delta, Sigma = Sigma, Beta0 = Beta0, Beta1 = Beta1, Lambda = Lambda, Eta = Eta, Theta = Theta)
  return(Out)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis_novar <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  M <- DatObj$M

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropLambdaVec <- MetropRcpp$MetropLambdaVec
  AcceptanceLambdaVec <- MetropRcpp$AcceptanceLambdaVec
  MetropEtaVec <- MetropRcpp$MetropEtaVec
  AcceptanceEtaVec <- MetropRcpp$AcceptanceEtaVec
  MetropAlpha <- MetropRcpp$MetropAlpha
  AcceptanceAlpha <- MetropRcpp$AcceptanceAlpha
  OriginalTuners <- MetrObj$OriginalTuners

  ###Summarize and output
  TuningParameters <- c(MetropLambdaVec, MetropEtaVec, MetropAlpha)
  AcceptanceCount <- c(AcceptanceLambdaVec, AcceptanceEtaVec, AcceptanceAlpha)
  AcceptancePcts <- AcceptanceCount / NSims
  MetrSummary <- cbind(AcceptancePcts, TuningParameters, OriginalTuners)
  rownames(MetrSummary) <- c(paste0("Lambda(", 1:M, ")"), paste0("Eta(", 1:M, ")"), "Alpha")
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  return(MetrSummary)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.spCP_novar
#'
#' \code{is.spCP_novar} is a general test of an object being interpretable as a
#' \code{\link{spCP_novar}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{spCP_novar}} class is defined as the regression object that
#'  results from the \code{\link{spCP_novar}} regression function.
#'
#' @export
is.spCP_novar <- function(x) {
  identical(attributes(x)$class, "spCP_novar")
}


