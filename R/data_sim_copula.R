#' Simulate bivariate time-to-event times based on specific copula and marginal distributions
#'
#' This function (data_sim_copula) is to generate a sample of subjects with two correlated event times based on specific copula and marginal models
#'
#' @name data_sim_copula
#' @aliases data_sim_copula
#' @param n sample size
#' @param copula types of copula, including "Clayton", "Gumbel", "Frank", "AMH" and "Joe"
#' @param eta copula parameter
#' @param dist marginal distributions, including "Weibull", "Gompertz" and "Loglogistic"
#' @param baseline marginal distribution parameters. For weibull and loglogistic, it shall be lambda (scale) and k (shape); for gompertz, it shall be a (shape) and b (rate)
#' @param var_list a vector of covariate names; assume same covariate names for two margins
#' @param COV_beta a vector of regression coefficients corresponding to var_list; assume same COV_beta between two margins
#' @param x1 a data frame of covariates for margin 1; it shall have n rows, with columns corresponding to the var_list
#' @param x2 a data frame of covariates for margin 2
#' @importFrom copula archmCopula
#' @importFrom copula rCopula
#' @return a data frame of bivariate time-to-event data with covariate
#' @export
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


data_sim_copula <-function(n, copula, eta, dist, baseline, var_list,
                           COV_beta, x1, x2) {

  # simulate
  cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
  Cop <- rCopula(n, cl)
  u<-Cop[,1]
  v<-Cop[,2]

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
  COV_beta <- matrix(COV_beta,ncol = 1)
  colnames(x1) = colnames(x2) = var_list
  dat1<-data.frame(id=seq(1:n),ind=rep(1,n),x1)
  dat2<-data.frame(id=seq(1:n),ind=rep(2,n),x2)


  # solve for t1 and t2
  if (dist == "Weibull") {
    t1<- as.vector(((-log(u))/(exp(matrix(unlist(dat1[,var_list]),ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
    t2<- as.vector(((-log(v))/(exp(matrix(unlist(dat2[,var_list]),ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
  }

  if (dist == "Loglogistic") {
    t1<- as.vector(((-1+1/u)/(exp(matrix(unlist(dat1[,var_list]),ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
    t2<- as.vector(((-1+1/v)/(exp(matrix(unlist(dat2[,var_list]),ncol = length(var_list))%*%COV_beta)))^(1/k)/lambda)
  }

  if (dist == "Gompertz") {
    t1<- as.vector(1/a*log(1-a*log(u)*exp(-1*matrix(unlist(dat1[,var_list]),ncol = length(var_list))%*%COV_beta)/b))
    t2<- as.vector(1/a*log(1-a*log(v)*exp(-1*matrix(unlist(dat2[,var_list]),ncol = length(var_list))%*%COV_beta)/b))

  }


  dat1<-cbind(dat1,data.frame(time=t1))
  dat2<-cbind(dat2,data.frame(time=t2))

  dat<-rbind(dat1,dat2)
  dat<-dat[order(dat$id),]

  return(dat)
}

