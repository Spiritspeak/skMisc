% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skUtils.R
\name{tokens_compound_stepwise}
\alias{tokens_compound_stepwise}
\title{Compound tokens without overflowing memory and crashing R}
\usage{
tokens_compound_stepwise(
  x,
  pattern,
  stepsize = 100,
  concatenator = "_",
  valuetype = c("glob", "regex", "fixed"),
  case_insensitive = TRUE,
  join = TRUE
)
}
\description{
A wrapper around \link[quanteda]{tokens_compound} that processes your tokens in chunks, 
set by argument \code{stepsize}. See \link[quanteda]{tokens_compound} for more info.
}
\examples{
toks<-tokens(data_corpus_inaugural)
compounded<-tokens_compound_stepwise(x=toks,pattern="I am",stepsize=10)

#note: does not work?

}
