\name{data_sim_ic}
\alias{data_sim_ic}
\docType{data}
\title{
  data_sim_ic
}
\description{
  A simulated real dataset of bivariate interval-censored data with 500 subjects and 7 columns.
}
\usage{data("data_sim_ic")}
\format{
  A data frame with 500 subjects on the following 7 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{Left.L}}{left bound of the first event}
    \item{\code{Right.L}}{right bound of the first event}
    \item{\code{status.L}}{censoring indicator of the first event; 1=interval censor, 0=right censor.}
    \item{\code{Left.R}}{left bound of the second event}
    \item{\code{Right.R}}{right bound of the second event}
    \item{\code{status.R}}{censoring indicator of the second event; 1=interval censor, 0=right censor.}
  }
}

\keyword{datasets}
