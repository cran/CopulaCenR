\name{ACTG181}
\alias{ACTG181}
\docType{data}
\title{
ACTG181
}
\description{
An example real dataset of bivariate interval-censored data with 204 subjects and one covariate. The data come from the AIDS Clinical Trials Group protocol ACTG 181. Two events are the shedding time (in months) of cytomegalovirus (CMV) in urine/blood and the colonization time (in months) of mycobacterium avium complex (MAC) in sputum/stool (Betensky and Finkelstein, 1999).
}
\usage{data("ACTG181")}
\format{
  A data frame with 408 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{ind}}{margin indicator, 1=shedding time, 2=colonization time}
    \item{\code{Left}}{left bound of observed interval}
    \item{\code{Right}}{right bound of observed interval}
    \item{\code{status}}{censoring indicator; 1=interval-censor, 0=right censor. All observations are interval-censored in this dataset}
    \item{\code{x}}{covariate}
  }
}
\source{
Betensky and Finkelstein (1999). A non-parametric maximum likelihood estimator for bivariate interval censored data. \emph{Statistics in Medicine} \strong{18} 3089-3100.
}
\examples{
data(ACTG181)
clayton_sp <- ic_spTran_copula(data = ACTG181,
                           copula = "Clayton",
                           l = 0, u = 25,
                           r = 3, m = 3,
                           var_list = "x")
summary(clayton_sp)
}
\keyword{datasets}
