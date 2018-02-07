###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples_lms <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  PhiIndeces <- DatObj$PhiIndeces
  tNu <- DatObj$tNu

  ###Format raw samples
  RawSamples <- t(RawSamples)
  Alpha <- RawSamples[, 1:5]
  Delta <- RawSamples[, 6:10]
  Sigma <- RawSamples[, 11:25]
  Phi <- RawSamples[, 26:((5 * M) + 25)]
  Beta0 <- Phi[, PhiIndeces[, 1]]
  Beta1 <- Phi[, PhiIndeces[, 2]]
  Lambda0 <- Phi[, PhiIndeces[, 3]]
  Lambda1 <- Phi[, PhiIndeces[, 4]]
  Eta <- Phi[, PhiIndeces[, 5]]
  Theta <- apply(Eta, 2, function(x) pmin(tNu, exp(x)))
  colnames(Alpha) <- paste0("Alpha", 1:5)
  colnames(Delta) <- paste0("Delta", 1:5)
  colnames(Sigma) <- c(paste0("Sigma", 1:5, "1"), paste0("Sigma", 2:5, "2"), paste0("Sigma", 3:5, "3"), paste0("Sigma", 4:5, "4"), "Sigma55")
  colnames(Beta0) <- paste0("Beta0(", 1:M, ")")
  colnames(Beta1) <- paste0("Beta1(", 1:M, ")")
  colnames(Lambda0) <- paste0("Lambda0(", 1:M, ")")
  colnames(Lambda1) <- paste0("Lambda1(", 1:M, ")")
  colnames(Eta) <- paste0("Eta(", 1:M, ")")
  colnames(Theta) <- paste0("Theta(", 1:M, ")")
  Out <- list(Alpha = Alpha, Delta = Delta, Sigma = Sigma, Beta0 = Beta0, Beta1 = Beta1, Lambda0 = Lambda0, Lambda1 = Lambda1, Eta = Eta, Theta = Theta)
  return(Out)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis_lms <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  M <- DatObj$M

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropBeta0Vec <- MetropRcpp$MetropBeta0Vec
  AcceptanceBeta0Vec <- MetropRcpp$AcceptanceBeta0Vec
  MetropBeta1Vec <- MetropRcpp$MetropBeta1Vec
  AcceptanceBeta1Vec <- MetropRcpp$AcceptanceBeta1Vec
  MetropLambda0Vec <- MetropRcpp$MetropLambda0Vec
  AcceptanceLambda0Vec <- MetropRcpp$AcceptanceLambda0Vec
  MetropLambda1Vec <- MetropRcpp$MetropLambda1Vec
  AcceptanceLambda1Vec <- MetropRcpp$AcceptanceLambda1Vec
  MetropEtaVec <- MetropRcpp$MetropEtaVec
  AcceptanceEtaVec <- MetropRcpp$AcceptanceEtaVec
  MetropSigma <- MetropRcpp$MetropSigma
  AcceptanceSigma <- MetropRcpp$AcceptanceSigma
  MetropAlpha <- MetropRcpp$MetropAlpha
  AcceptanceAlpha <- MetropRcpp$AcceptanceAlpha
  OriginalTuners <- MetrObj$OriginalTuners

  ###Summarize and output
  TuningParameters <- c(MetropBeta0Vec, MetropBeta1Vec, MetropLambda0Vec, MetropLambda1Vec, MetropEtaVec, MetropSigma, MetropAlpha)
  AcceptanceCount <- c(AcceptanceBeta0Vec, AcceptanceBeta1Vec, AcceptanceLambda0Vec, AcceptanceLambda1Vec, AcceptanceEtaVec, AcceptanceSigma, AcceptanceAlpha)
  AcceptancePcts <- AcceptanceCount / NSims
  MetrSummary <- cbind(AcceptancePcts, TuningParameters, OriginalTuners)
  rownames(MetrSummary) <- c(paste0("Beta0(", 1:M, ")"), paste0("Beta1(", 1:M, ")"), paste0("Lambda0(", 1:M, ")"), paste0("Lambda1(", 1:M, ")"), paste0("Eta(", 1:M, ")"), paste0("Sigma", 1:5, "1"), paste0("Sigma", 2:5, "2"), paste0("Sigma", 3:5, "3"), paste0("Sigma", 4:5, "4"), "Sigma55", paste0("Alpha", 1:5))
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  return(MetrSummary)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.spCP_lmc
#'
#' \code{is.spCP_lmc} is a general test of an object being interpretable as a
#' \code{\link{spCP_lmc}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{spCP_lmc}} class is defined as the regression object that
#'  results from the \code{\link{spCP_lmc}} regression function.
#'
#' @export
is.spCP_lmc <- function(x) {
  identical(attributes(x)$class, "spCP_lmc")
}



