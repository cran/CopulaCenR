\name{CopulaCenR}
\alias{CopulaCenR}
\title{
Copula-based regression models for bivariate censored data
}
\description{
Bivariate time-to-event data frequently arise in research areas such as clinical trials and epidemiological studies, where the occurrence of two events are correlated. In many cases, the exact event times are unknown due to censoring. Specifically, bivariate right-censored data occur when the study ends prior to the occurrence of one or both events. In another situation, bivariate interval-censored data occur when the status of both events are periodically examined at intermittent assessment times. In this case, the right censoring could also happen if the event still does not occur at the last assessment time. A special case of interval-censored data is the current status data if there is only one assessment time and the event is only known to occur or not by its assessment time. \cr

The copula model is a popular approach for modeling correlated bivariate censored data. One unique feature of copula is that it models the two marginal distributions and the between-margin dependence separately, allowing flexibility in marginal models and straightforward interpretation for covariate effects. Moreover, the challenge from censoring can be naturally handled through the marginal distributions within the copula function. Besides, the joint, marginal and conditional distributions can be obtained based on the copula model. However, there is a lack of R package implementing copula-based regression models for bivariate data under both right- and interval-censoring. \cr

The CopulaCenR package can build copula-based regression models for both bivariate right-censored data and bivariate interval-censored data (including the special case of bivariate current status data). The package is flexible in terms of various built-in Archimedean copulas (such as Clayton, Gumbel, Frank, Joe, AMH), together with a two-parameter copula (Copula2) that incorporates the popular Clayton and Gumbel copulas as special cases. It also allows a broad class of marginal distributions, including parametric (i.e.,Weibull, Gompertz, Loglogistic) and semiparametric (i.e., unspecified functions, such as Cox, transformation) margins. \cr

The main features of CopulaCenR are: \cr
(1) Copula models with parametric margins: \code{rc_par_copula} and \code{ic_par_copula}) for bivariate right-censored data and bivariate interval-censored data, respectively. \cr

(2) Copula models with semiparametric margins: \code{rc_spCox_copula} and \code{ic_spTran_copula}) for bivariate right-censored data and bivariate interval-censored data, respectively. \cr

(3) Wald (performed during model fitting), score (\code{score_copula}) and likelihood-ratio (\code{lrt_copula}) tests for covariate effects. \cr

(4) Calculate Kendall's tau from a fitted model: \code{tau_copula}. \cr

(5) Calculate other model parameters using S3 methods: \code{print}, \code{summary}, \code{coef}, \code{logLik}, \code{AIC}, \code{BIC}. \cr

(6) Extract fitted values (such as linear predictors and survival probabilities) from a fitted model: \code{fitted}. \cr

(7) Predict in new observations (such as linear predictors and survival probabilities) from a fitted model: \code{predict}. \cr

(8) Plot joint, conditional and marginal survival probability profiles for given subjects based on a fitted model: \code{plot}, \code{lines}. \cr

(9) A user-friendly function to simulate bivariate event times: \code{data_sim_copula}. \cr
}
\source{
Tao Sun, Yi Liu, Richard J. Cook, Wei Chen and Ying Ding (2018). Copula-based Score Test for Bivariate Time-to-event Data, with Application to a Genetic Study of AMD Progression. #' \emph{Lifetime Data Analysis} 25(3), 546???568. \cr

Tao Sun and Ying Ding (In Press). Copula-based Semiparametric Regression Method for Bivariate Data under General Interval Censoring. http://arxiv.org/abs/1901.01918. \cr
}
\examples{
### bivariate right-censored data
data(DRS)
# fit a Clayton-Weibull model
clayton_wb <- rc_par_copula(data = DRS,
                            var_list = "treat",
                            copula = "Clayton",
                            m.dist = "Weibull")
summary(clayton_wb)

### bivariate interval-censored data
data(AREDS)
copula2_sp <- ic_spTran_copula(data = AREDS,
              copula = "Copula2", l = 0, u = 15,
              m = 3, r = 3, var_list =
              c("ENROLLAGE", "rs2284665",
              "SevScaleBL"))
summary(copula2_sp)
}

