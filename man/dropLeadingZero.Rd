% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatting.R
\name{dropLeadingZero}
\alias{dropLeadingZero}
\title{Drop leading zeros
Remove leading zeroes and return as a character object.}
\usage{
dropLeadingZero(x)
}
\arguments{
\item{x}{Numeric vector, matrix, or list to remove leading zeros from.}
}
\value{
A character vector, matrix, or list, with numbers with leading zeros removed.
}
\description{
Drop leading zeros
Remove leading zeroes and return as a character object.
}
\examples{
dropLeadingZero(c(-1,0,1,0.5,-0.5,1.5,-1.5))

}
