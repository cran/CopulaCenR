\name{data_sim_RC}
\alias{data_sim_RC}
\docType{data}
\title{
  data_sim_RC
}
\description{
  A simulated real dataset of bivariate right-censored data with 500 subjects and 5 columns.
}
\usage{data("data_sim_RC")}
\format{
  A data frame with 500 subjects on the following 5 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{event_time.L}}{observed time of the first event}
    \item{\code{status.L}}{censoring indicator of the first event; 1=exactly observed, 0=right censor.}
    \item{\code{event_time.R}}{observed time of the second event}
    \item{\code{status.R}}{censoring indicator of the second event; 1=exactly observed, 0=right censor.}
  }
}

\keyword{datasets}
