% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{compcorr}
\alias{compcorr}
\title{Test if two correlation coefficients significantly differ}
\usage{
compcorr(cor1, cor2, n1, n2)
}
\arguments{
\item{cor1, cor2}{Correlation values being compared}

\item{n1, n2}{Sample sizes of the correlation coefficients}
}
\value{
List containing the z-score and p-value
}
\description{
Uses Fisher's r to z transformation, then performs a z-test on the resulting z-scores
}
\examples{
compcorr(.1,.6,50,100)
}
\references{
http://vassarstats.net/rdiff.html
}
