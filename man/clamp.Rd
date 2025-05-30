% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{clamp}
\alias{clamp}
\alias{clamp0}
\alias{scale2range}
\title{Clamp}
\usage{
clamp(x, minval = -Inf, maxval = Inf)

clamp0(x, minval = 0, maxval = 1)

scale2range(x, minval = 0, maxval = 1, na.rm = TRUE)
}
\arguments{
\item{x}{The vector/matrix to clamp.}

\item{minval}{Minimum value; all lower values are clamped to this value.}

\item{maxval}{Maximum value; all higher values are clamped to this value.}

\item{na.rm}{for \code{scale2range}, whether \code{NA}s should be removed before
applying the function (default \code{TRUE}. 
If \code{FALSE} and \code{x} contains \code{NA}s, then the
complete vector will be \code{NA}.}
}
\value{
Clamped vector.
}
\description{
Clamp a numeric vector between a minimum and maximum value.
}
\section{Functions}{
\itemize{
\item \code{clamp()}: Set all values exceeding the minimum and maximum to 
the minimum and maximum value, respectively.

\item \code{clamp0()}: Same as \code{clamp())} but with a default range of 0 to 1.

\item \code{scale2range()}: Rescales the vector such that it fits neatly between 
the given minimum and maximum values.

}}
\examples{
clamp(0:10,2,8)
clamp0(rnorm(10))
scale2range(rnorm(10))

}
\author{
Sercan Kahveci
}
