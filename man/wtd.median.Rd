% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{wtd.median}
\alias{wtd.median}
\title{Weighted Median}
\usage{
wtd.median(x, wts, na.rm = T)
}
\arguments{
\item{x}{an input vector}

\item{wts}{a vector of weights}

\item{na.rm}{Logical indicating whether NA values in the input and weight vectors should be stripped.}
}
\value{
A weighted median of the input values and weights.
}
\description{
Weighted Median
}
\examples{
wtd.median(1:5,c(.5,4,1,2,1))

}
