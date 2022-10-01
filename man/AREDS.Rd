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
  A data frame with 1258 observations (629 subjects with 2 eyes) on the following 8 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{ind}}{margin indicator, 1=left eye, 2=right eye}
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
copula2_sp <- ic_spTran_copula(data = AREDS,
              copula = "Copula2",
              var_list = c("ENROLLAGE",
                           "rs2284665",
                           "SevScaleBL"),
              l = 0, u = 15, m = 3, r = 3)
summary(copula2_sp)
}
\keyword{datasets}
