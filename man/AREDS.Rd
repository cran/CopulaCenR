\name{AREDS}
\alias{AREDS}
\docType{data}
\title{
AREDS
}
\description{
A real dataset of bivariate interval-censored data with 629 subjects and 4 non-genetic covariates and 1 genetic covariate. The dataset is selected from the Age-related Eye Disease Study (AREDS) (AREDS Group, 1999). Two events are the progression times (in years) to late-AMD in the left and right eyes.
}
\usage{data("AREDS")}
\format{
  A data frame with 1258 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{ind}}{margin indicator, 1=shedding time, 2=colonization time}
    \item{\code{Left}}{left bound of observed interval}
    \item{\code{Right}}{right bound of observed interval}
    \item{\code{status}}{censoring indicator; 1=interval-censor, 0=right censor.}
    \item{\code{SevScaleBL}}{baseline AMD severity score, margin-specific}
    \item{\code{ENROLLAGE}}{age at baseline}
    \item{\code{rs2284665}}{a SNP covariate highly associated with late-AMD progression, coded as 0,1,2}
  }
}
\source{
AREDS Group (1999). The Age-Related Eye Disease Study (AREDS): design implications. AREDS report no. 1. \emph{Control Clinical Trials} \strong{20}, 573-600.
}
\examples{
data(AREDS)
}
\keyword{datasets}
