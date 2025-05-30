% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{bundle_up}
\alias{bundle_up}
\title{Assign numeric values to bundles summing up to a set value}
\usage{
bundle_up(
  x,
  maxval,
  group.high = c("own", "na", "error"),
  group.na = c("own", "na", "same", "error")
)
}
\arguments{
\item{x}{A numeric vector to bundle up.}

\item{maxval}{The maximum value a bundle can have.}

\item{group.high}{What should be done if an individual value exceeds \code{maxval}?
\code{"own"} assigns it to its own bundle, \code{"na"} assigns it to \code{NA}, and 
\code{"error"} gives an error.}

\item{group.na}{What should be done if a value is \code{NA}? 
\code{"own"} assigns it to its own bundle, \code{"na"} assigns it to \code{NA}, 
\code{"same"} assigns it to the same bundle as the previous value, and 
\code{"error"} gives an error.}
}
\value{
An integer vector of the same length as \code{x}, where each value indicates
to which bundle the value in \code{x} has been assigned.
}
\description{
Assign numeric values to bundles summing up to a set value
}
\examples{
testvec <- sample(c(NA,NA,1:14))
t(cbind(testvec,bundle_up(testvec,maxval=12,group.high="own",group.na="own")))

bundle_up(testvec,maxval=10,group.high="na",group.na="na")

}
\author{
Sercan Kahveci
}
