\name{data_sim_rec}
\alias{data_sim_rec}
\docType{data}
\title{
  data_sim_rec
}
\description{
  A simulated real dataset of bivariate recurrent events data with 500 subjects and 6 columns.
}
\usage{data("data_sim_rec")}
\format{
  A data frame with 500 subjects on the following 6 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{gap1}}{observed time of the first gap time}
    \item{\code{status1}}{censoring indicator of the first event; 1=exactly observed, 0=right censor.}
    \item{\code{gap2}}{observed time of the second gap time}
    \item{\code{status2}}{censoring indicator of the second event; 1=exactly observed, 0=right censor.}
    \item{\code{d}}{cluster size}
  }
}

\keyword{datasets}
