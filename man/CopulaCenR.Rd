\name{CopulaCenR}
\alias{CopulaCenR}
\title{
Copula-based regression models for bivariate censored data
}
\description{
Copula-based regression models for both bivariate right-censored data and bivariate interval-censored data: \cr
1. Copula models with parametric (e.g. Weibull) margins: estimation/wald test (\code{rc_par_copula}, \code{ic_par_copula}), score test (\code{score_copula}) and likelihood ratio test (\code{lrt_copula}); \cr
2. Copula models with semiparametric sieve margins: estimation/wald test (\code{ic_sp_copula}), score test (\code{score_copula}) and likelihood ratio test (\code{lrt_copula}). \cr
3. Predict and plot joint/conditional/marginal survival probabilities based on a fitted copula model using \code{fitted}, \code{predict}, \code{plot}, \code{lines}. \cr
4. A user-friendly function (\code{data_sim_copula}) to simulate correlated event times based on copula and marginal specifications. \cr

Data must be a data frame. For bivariate right-censored data, it shall include \code{id} (subject id), \code{ind} (1,2 to indicate units of one subject), \code{obs_time}, \code{status} (0 for right-censoring and 1 for event). \cr
For bivariate interval-censored data, it shall include \code{id} (subject id), \code{ind} (1,2 to indicate units of one subject), \code{Left} (great than 0), \code{Right} (can include \code{Inf}), \code{status} (0 for right-censoring and 1 for left- and interval-censoring). \cr

}
\source{
Tao Sun, Yi Liu, Richard J. Cook, Wei Chen and Ying Ding (2018). Copula-based Score Test for Bivariate Time-to-event Data, with Application to a Genetic Study of AMD Progression. \emph{Lifetime Data Analysis} doi:10.1007/s10985-018-09459-5. \cr
Tao Sun and Ying Ding (2019). Copula-based Semiparametric Transformation Model for Bivariate Data Under General Interval Censoring. http://arxiv.org/abs/1901.01918. \cr
}
\examples{
### bivariate right-censored data
data(DRS)
# fit a Clayton-Weibull model
clayton_wb <- rc_par_copula(data = DRS, var_list = "treat",
              copula = "Clayton", m.dist = "Weibull")
summary(clayton_wb)

### bivariate interval-censored data
data(AREDS)
copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2",
              l = 0, u = 15, m = 3, r = 3,
              var_list = c("ENROLLAGE","rs2284665","SevScaleBL"))
summary(copula2_sp)
}

