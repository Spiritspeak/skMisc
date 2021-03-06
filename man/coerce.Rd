% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{coerce}
\alias{coerce}
\title{coerce a vector to contain only TRUE and FALSE}
\usage{
coerce(x, default = FALSE)
}
\arguments{
\item{x}{Numeric/logical vector/matrix to coerce into TRUE/FALSE}

\item{default}{default returned value if NULL or NA is encountered}
}
\value{
logical vector or matrix with only T and F
}
\description{
coerce a vector to contain only TRUE and FALSE
}
\examples{
coerce(NULL)
# FALSE

coerce(c(T,F,NA,NA,T))
# T F F F T

coerce(matrix(c(T,T,F,F,NA,NA),nrow=2))
#     [,1]  [,2]  [,3]
#[1,] TRUE FALSE FALSE
#[2,] TRUE FALSE FALSE
}
