% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{vec2columns}
\alias{vec2columns}
\title{Split a character column into multiple values}
\usage{
vec2columns(x, sep = ";")
}
\arguments{
\item{x}{a character vector to split into columns.}

\item{sep}{a character separating the different values.}
}
\value{
a \code{data.frame} of boolean values, with each row representing 
a value of x and each column representing a unique value 
in \code{x} following splitting. A column is marked TRUE in a specific row if
the value representing that column was present in that row.
}
\description{
Split a character column into multiple values
}
\examples{
unsplit<-c("flour;salt;baking soda;steak;sugar;water;sauce;vinegar",
"flour;sauce;mustard;salt;pepper;vinegar;baking soda;water;tomatoes;onion;steak")
vec2columns(unsplit)

}
\author{
Sercan Kahveci
}
