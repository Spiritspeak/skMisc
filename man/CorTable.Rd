% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{CorTable}
\alias{CorTable}
\title{Create a Correlation Table}
\usage{
CorTable(df, rowids, columnids, rowdf, columndf)
}
\arguments{
\item{df}{A data.frame.}

\item{rowids, columnids}{character vectors containing column names from \code{df} that need to be correlated.}

\item{rowdf, columndf}{data.frames whose columns need to be correlated. 
Either \code{df, rowids, & columnids} or \code{rowdf & columndf} are required.}
}
\value{
A formatted markdown table containing correlation coefficients, p-values, and 
the number and percentage of cases that need to be removed to flip the sign of each correlation coefficient.
}
\description{
Create a Correlation Table
}
\examples{
CorTable(mtcars,rowids=c("mpg","disp","hp"),columnids=c("drat","wt","qsec"))

CorTable(rowdf=mtcars[,c(1,3,4)],columndf=mtcars[,5:7])
}
