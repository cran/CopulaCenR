\name{DRS}
\alias{DRS}
\docType{data}
\title{
DRS
}
\description{
A real dataset of bivariate right-censored data with 197 subjects and 3 covariates. The patients were a 50\% random sample of the patients with "high-risk" diabetic retinopathy as defined by the Diabetic Retinopathy Study (DRS) (Huster, 1989).
}
\usage{data("DRS")}
\format{
  A data frame with 394 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{ind}}{margin indicator, 1 for right eye, 2 for left eye}
    \item{\code{obs_time}}{time of blindness (in months) from treatment}
    \item{\code{status}}{censoring indicator, 1=blindness, 0=right censoring}
    \item{\code{treat}}{laser treatment type, 0=no treatment, 1=xenon, 2=argon}
    \item{\code{age}}{age at diagnosis of diabetes}
    \item{\code{type}}{type of diabetes, 1= juvenile (age at treatment < 20), 2=adult}
  }
}
\details{
Each patient had one eye randomized to laser treatment and the other eye received no treatment. For each eye, the event of interest was the time from initiation of treatment to the time when visual acuity dropped below 5/200 two visits in a row (call it "blindness"). Thus there is a built-in lag time of approximately 6 months (visits were every 3 months). Survival times in this dataset are therefore the actual time to blindness in months, minus the minimum possible time to event (6.5 months). Censoring was caused by death, dropout, or end of the study.
}
\source{
https://www.mayo.edu/research/documents/DRShtml/DOC-10027460
}
\references{
Huster WJ, Brookmeyer R, Self SG (1989). Modelling paired survival data with covariates. \emph{Biometrics} \strong{45}, 145-156.
}
\examples{
data(DRS)
}
\keyword{datasets}
