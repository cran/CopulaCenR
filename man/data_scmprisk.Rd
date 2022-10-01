\name{data_scmprisk}
\alias{data_scmprisk}
\docType{data}
\title{
  data_scmprisk
}
\description{
  A simulated real dataset of interval-censored semi-competing risks data with 500 subjects and 10 columns.
}
\usage{data("data_scmprisk")}
\format{
  A data frame with 500 subjects on the following 10 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{Left}}{left bound of observed interval of non-terminal event}
    \item{\code{Right}}{right bound of observed interval of non-terminal event}
    \item{\code{status}}{censoring indicator of non-terminal event; 1=interval-censor, 0=right censor.}
    \item{\code{timeD}}{observed time of non-terminal event}
    \item{\code{statusD}}{censoring indicator of terminal event; 1=interval-censor, 0=right censor.}
    \item{\code{A}}{study entry time (i.e., left truncation time).}
    \item{\code{x1}}{first covariate}
    \item{\code{x2}}{second covariate}
    \item{\code{x3}}{third covariate}
  }
}

\keyword{datasets}
