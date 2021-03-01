\name{Kidney}
\alias{Kidney}
\docType{data}
\title{
Kidney
}
\description{
A real dataset of bivariate right-censored data with 38 subjects and 3 covariates. The data are the recurrence times to infection, at the point of insertion of the catheter, for kidney patients using portable dialysis equipment. Catheters may be removed for reasons other than infection, in which case the observation is censored. Each patient has exactly 2 observations.
}
\usage{data("Kidney")}
\format{
  A data frame with 76 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{ind}}{margin indicator}
    \item{\code{obs_time}}{observed time}
    \item{\code{status}}{event status}
    \item{\code{age}}{in years}
    \item{\code{sex}}{1=male, 2=female}
    \item{\code{disease}}{disease type with 4 levels \code{Other} \code{GN} \code{AN} \code{PKD}}
  }
}
\note{
This data has often been used to illustrate the use of random effects (frailty) in a survival model. However, one of the males (id 21) is a large outlier, with much longer survival than his peers. If this observation is removed no evidence remains for a random subject effect.
}
\source{
https://github.com/therneau/survival
}
\references{
CA McGilchrist, CW Aisbett (1991), Regression with frailty in survival analysis. \emph{Biometrics} \strong{47}, 461-66.
}
\examples{
data(Kidney)
clayton_cox <- rc_spCox_copula(data = Kidney,
               var_list = c("age","sex","disease"),
               copula = "Clayton",
               method = "BFGS",
               B = 2)
summary(clayton_cox)
}
\keyword{datasets}
