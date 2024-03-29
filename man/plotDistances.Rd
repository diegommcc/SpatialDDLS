% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSpatialCoor.R
\name{plotDistances}
\alias{plotDistances}
\title{Plot distances between intrinsic and extrinsic profiles}
\usage{
plotDistances(
  object,
  index.st,
  mid.scale = "mean",
  size.point = 1,
  title = NULL,
  theme = NULL
)
}
\arguments{
\item{object}{A \code{\linkS4class{SpatialDDLS}} object.}

\item{index.st}{Index of the spatial transcriptomics data to be plotted. It
can be either a position or a name if a named list was provided.}

\item{mid.scale}{The midpoint of the diverging scale. it may be \code{'mean'}
or \code{'median'} (the former by default).}

\item{size.point}{Size of points (0.1 by default).}

\item{title}{Title of plot.}

\item{theme}{\pkg{ggplot2} theme.}
}
\value{
A ggplot object.
}
\description{
Color spots on the spatial coordinates according to distances between
intrinsic and extrinsic transcriptional profiles.
}
\seealso{
\code{\link{deconvSpatialDDLS}} \code{\link{trainDeconvModel}}
}
