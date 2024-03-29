% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_CP.R
\name{PlotCP}
\alias{PlotCP}
\title{PlotCP}
\usage{
PlotCP(
  object,
  data,
  location = "Location",
  time = "Time",
  dls = "DLS",
  line = TRUE,
  ci = TRUE,
  lwd = 1,
  lty = 1,
  col = 2,
  ci.lwd = 1,
  ci.lty = 2,
  ci.col = 2,
  cp.line = FALSE,
  cp.ci = FALSE,
  cp.lwd = 1,
  cp.lty = 1,
  cp.col = 4,
  cp.ci.lwd = 1,
  cp.ci.lty = 2,
  cp.ci.col = 4,
  main = "Estimated visual field sensitivity\\nusing change points",
  xlab = "Years from baseline visit",
  ylab = "Sensitivity (dB)"
)
}
\arguments{
\item{object}{a spCP regression object.}

\item{data}{a dataframe containing the raw sensitivities.}

\item{location}{a character string indicating the column name in data of the variable of locations on the visual field.}

\item{time}{a character string indicating the column name in data of the variable of visual field testing times.}

\item{dls}{a character string indicating the column name in data of the variable of the raw visual field sensitivities.}

\item{line}{logical, determines if there are regression lines printed (default = TRUE).}

\item{ci}{logical, determines if there are confidence intervals printed (default = TRUE).}

\item{lwd}{integer, specifies the width of the regression line (default = 1).}

\item{lty}{integer, specifies the type of regression line (default = 1).}

\item{col}{color for the regression line, either character string corresponding
to a color or a integer (default = "red").}

\item{ci.lwd}{integer, specifies the width of the confidence intervals (default = 1).}

\item{ci.lty}{integer, specifies the type of confidence intervals (default = 2).}

\item{ci.col}{color for the confidence intervals for the regression line, either character string corresponding
to a color or a integer (default = "red").}

\item{cp.line}{logical, determines if there the change point is printed (default = FALSE).}

\item{cp.ci}{logical, determines if there the change point confidence intervals are printed (default = FALSE).}

\item{cp.lwd}{integer, specifies the width of the change point line (default = 1).}

\item{cp.lty}{integer, specifies the type of change point line (default = 1).}

\item{cp.col}{color for the optional change point, either character string corresponding
to a color or a integer (default = "blue").}

\item{cp.ci.lwd}{integer, specifies the width of the change point confidence intervals (default = 1).}

\item{cp.ci.lty}{integer, specifies the type of change point confidence intervals (default = 2).}

\item{cp.ci.col}{color for the confidence intervals of the optional change point, either character string corresponding
to a color or a integer (default = "blue").}

\item{main}{an overall title for the plot.}

\item{xlab}{a title for the x axis.}

\item{ylab}{a title for the y axis.}
}
\description{
Plots estimated visual field sensitivities using change point model.
}
\details{
\code{PlotCP} is used in the application of glaucoma progression.
 The function is capable of plotting the observed DLS values across the visual field,
 along with the estimated mean process (with 95 percent credible intervals) and the
 estimated mean posterior change point location (with 95 percent credible intervals).
}
\author{
Samuel I. Berchuck
}
