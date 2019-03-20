#' Simulate bivariate time-to-event times based on specific copula and marginal distributions
#'
#' To generate a sample of subjects with two correlated event times based on specific copula and marginal models
#'
#' @name data_sim_copula
#' @aliases data_sim_copula
#' @param n sample size
#' @param copula types of copula, including \code{"Clayton"}, \code{"Gumbel"},
#' \code{"Frank"}, \code{"AMH"}, \code{"Joe"}
#' @param eta copula parameter \eqn{\eta}
#' @param dist marginal distributions, including \code{"Weibull"},
#' \code{"Gompertz"}, \code{"Loglogistic"}
#' @param baseline marginal distribution parameters.
#' For \code{Weibull} and \code{Loglogistic}, it shall be \eqn{\lambda} (scale) and \eqn{k} (shape);
#' for \code{Gompertz}, it shall be \eqn{a} (shape) and \eqn{b} (rate)
#' @param var_list a vector of covariate names;
#' assume the same covariates for two margins
#' @param COV_beta a vector of regression coefficients corresponding to \code{var_list};
#' assume the same coefficients between two margins
#' @param x1 a data frame of covariates for margin 1; it shall have n rows,
#' with columns corresponding to the \code{var_list}
#' @param x2 a data frame of covariates for margin 2
#' @importFrom copula archmCopula
#' @importFrom copula rCopula
#' @return a data frame of bivariate time-to-event data with covariates
#' @export
#'
#' @details The parametric generator functions of copula functions are list below:
#'
#' The Clayton copula has a generator \deqn{\phi_{\eta}(t) = (1+t)^{-1/\eta},}
#' with \eqn{\eta > 0} and Kendall's \eqn{\tau = \eta/(2+\eta)}.
#'
#' The Gumbel copula has a generator \deqn{\phi_{\eta}(t) = \exp(-t^{1/\eta}),}
#' with \eqn{\eta \geq 1} and Kendall's \eqn{\tau = 1 - 1/\eta}.
#'
#' The Frank copula has a generator \deqn{\phi_{\eta}(t) = -\eta^{-1}\log \{1+e^{-t}(e^{-\eta}-1)\},}
#' with \eqn{\eta \geq 0} and Kendall's \eqn{\tau = 1+4\{D_1(\eta)-1\}/\eta},
#' in which \eqn{D_1(\eta) = \frac{1}{\eta} \int_{0}^{\eta} \frac{t}{e^t-1}dt}.
#'
#' The AMH copula has a generator \deqn{\phi_{\eta}(t) = (1-\eta)/(e^{t}-\eta),}
#' with \eqn{\eta \in [0,1)} and Kendall's \eqn{\tau =  1-2\{(1-\eta)^2 \log (1-\eta) + \eta\}/(3\eta^2)}.
#'
#' The Joe copula has a generator \deqn{\phi_{\eta}(t) = 1-(1-e^{-t})^{1/\eta},}
#' with \eqn{\eta \geq 1} and Kendall's \eqn{\tau = 1 - 4 \sum_{k=1}^{\infty} \frac{1}{k(\eta k+2)\{\eta(k-1)+2\}}}. \cr
#'
#'
#'
#' The marginal survival distributions are listed below:
#'
#' The Weibull (PH) survival distribution is \deqn{\exp \{-(t/\lambda)^k  e^{Z^{\top}\beta}\},}
#' with \eqn{\lambda > 0} as scale and \eqn{k > 0} as shape.
#'
#' The Gompertz (PH) survival distribution is \deqn{\exp \{-\frac{b}{a}(e^{at}-1) e^{Z^{\top}\beta}\},}
#' with \eqn{a > 0} as shape and \eqn{b > 0} as rate
#'
#' The Loglogistic (PO) survival distribution is \deqn{\{1+(t/\lambda)^{k} e^{Z^{\top}\beta} \}^{-1},}
#' with \eqn{\lambda > 0} as scale and \eqn{k > 0} as shape.
#'
#' @examples
#' library(CopulaCenR)
#' set.seed(1)
#' dat <- data_sim_copula(n = 500, copula = "Clayton", eta = 3,
#'                        dist = "Weibull", baseline = c(0.1,2),
#'                        var_list = c("var1", "var2"),
#'                        COV_beta = c(0.1, 0.1),
#'                        x1 = cbind(rnorm(500, 6, 2),
#'                                   rbinom(500, 1, 0.5)),
#'                        x2 =  cbind(rnorm(500, 6, 2),
#'                                    rbinom(500, 1, 0.5)))
#' plot(x = dat$time[dat$ind == 1], y = dat$time[dat$ind == 2],
#'      xlab = expression(t[1]), ylab = expression(t[2]),
#'      cex.axis = 1, cex.lab = 1.3)


data_sim_copula <- function(n, copula, eta, dist, baseline,
                            var_list, COV_beta, x1, x2) {

  # simulate
  cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
  Cop <- rCopula(n, cl)
  u <- Cop[ , 1]
  v <- Cop[ , 2]

  # baseline parameters
  if (dist != "Gompertz"){
    lambda <- baseline[1]
    k <- baseline[2]
  }
  if (dist == "Gompertz"){
    a <- baseline[1]
    b <- baseline[2]
  }

  # regression parameters and design matrix
  COV_beta <- matrix(COV_beta, ncol = 1)
  colnames(x1) = colnames(x2) = var_list
  dat1 <- data.frame(id = seq(1:n),ind = rep(1,n), x1)
  dat2 <- data.frame(id = seq(1:n),ind = rep(2,n), x2)


  # solve for t1 and t2
  if (dist == "Weibull") {
    t1 <- as.vector(((-log(u))/(exp(matrix(unlist(dat1[, var_list]),
                                           ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
    t2 <- as.vector(((-log(v))/(exp(matrix(unlist(dat2[, var_list]),
                                           ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
  }

  if (dist == "Loglogistic") {
    t1 <- as.vector(((-1 + 1/u)/(exp(matrix(unlist(dat1[, var_list]),
                                         ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
    t2 <- as.vector(((-1 + 1/v)/(exp(matrix(unlist(dat2[, var_list]),
                                         ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
  }

  if (dist == "Gompertz") {
    t1 <- as.vector(1/a*log(1 - a*log(u)*exp(-1*matrix(unlist(dat1[, var_list]),
                                                       ncol = length(var_list))%*%COV_beta)/b))
    t2 <- as.vector(1/a*log(1 - a*log(v)*exp(-1*matrix(unlist(dat2[, var_list]),
                                                       ncol = length(var_list))%*%COV_beta)/b))

  }


  dat1 <- cbind(dat1, data.frame(time = t1))
  dat2 <- cbind(dat2, data.frame(time = t2))

  dat <- rbind(dat1, dat2)
  dat <- dat[order(dat$id), ]

  return(dat)
}

