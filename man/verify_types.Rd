% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{verify_types}
\alias{verify_types}
\title{Verify variable types in bulk}
\usage{
verify_types(...)
}
\arguments{
\item{...}{Named arguments, where the argument is the object to be checked and 
the name of the argument is the mode (numeric, list, character, etc).}
}
\value{
Returns true on success, causes error if not.
}
\description{
Verify variable types in bulk
}
\examples{
try(verify_types(character="test",numeric=0000,character=12345))

}
\author{
Sercan Kahveci
}
