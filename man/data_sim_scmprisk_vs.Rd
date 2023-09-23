\name{data_sim_scmprisk_vs}
\alias{data_sim_scmprisk_vs}
\docType{data}
\title{
  data_sim_scmprisk_vs
}
\description{
  A simulated real dataset of bivariate time-to-event data based on specific copula and marginal distributions with 500 subjects and 10 columns.
}
\usage{data("data_sim_scmprisk_vs")}
\format{
  A data frame with 500 subjects on the following 10 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{ind}}{1,2 for two margins}
    \item{\code{x1}}{first covariate}
    \item{\code{x2}}{second covariate}
    \item{\code{x3}}{third covariate}
    \item{\code{x4}}{fourth covariate}
    \item{\code{x5}}{fifth covariate}
    \item{\code{time}}{time to event}
    \item{\code{obs_time}}{observed time}
    \item{\code{status}}{0 for right-censoring, 1 for event}
  }
}

\keyword{datasets}
