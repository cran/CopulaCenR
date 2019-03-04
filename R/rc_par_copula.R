#' Copula regression models with parametric margins for bivariate right-censored data
#'
#' @description Fits a copula model with parametric margins for bivariate right-censored data.
#'
#' @name rc_par_copula
#' @aliases rc_par_copula
#' @param data a data frame; must have id (subject id), ind (1,2 for two margins), obs_time, status (0 for right-censoring, 1 for event).
#' @param var_list the list of covariates to be fitted into the model.
#' @param copula specify the copula family.
#' @param m.dist specify the marginal baseline distribution.
#' @param n.cons number of pieces, only for m.dist = "Piecewise". Default is 4.
#' @param method optimization method; default is "Newton" (from nlm).
#' @param iter number of iterations; default is 500.
#' @param stepsize size of optimization step; default is 1e-6.
#' @importFrom corpcor pseudoinverse
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival cluster
#' @importFrom stats as.formula
#' @importFrom stats cor
#' @importFrom stats optim
#' @importFrom stats pchisq
#' @importFrom stats quantile
#' @importFrom stats coef
#' @importFrom pracma grad
#' @importFrom pracma hessian
#' @source
#' Tao Sun, Yi Liu, Richard J. Cook, Wei Chen and Ying Ding (2018).
#' Copula-based Score Test for Bivariate Time-to-event Data, with Application to a Genetic Study of AMD Progression.
#' \emph{Lifetime Data Analysis} doi:10.1007/s10985-018-09459-5. \cr
#' Tao Sun and Ying Ding (2019).
#' Copula-based Semiparametric Transformation Model for Bivariate Data Under General Interval Censoring.
#' http://arxiv.org/abs/1901.01918.
#' @export
#'
#' @details The supported copula models are "Clayton", "Gumbel", "Frank", "AMH", "Joe", "Copula2".
#' The "Copula2" model is a two-parameter copula model that incorporates Clayton and Gumbel as special cases.
#' The supported marginal distributions are "Weibull" (proportional hazards), "Gompertz" (proportional hazards),
#' "Piecewise" (proportional hazards) and "Loglogistic" (proportional odds).
#' We assume the same baseline parameters between two margins. \cr
#'
#' The input data must be a data frame with columns id (subject id), ind (1,2 for two margins; each id must have both ind = 1 and 2), obs_time, status (0 for right-censoring, 1 for event) and covariates. \cr
#'
#' Optimization methods can be "Newton" (from nlm) and all methods from optim, such as "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN". "Brent" is not allowed due to more than one estimators.
#'
#' @return a CopulaCenR object summarizing the model. Can be used as an input to general S3 methods including
#' summary, plot, lines, coef, logLik, predict.
#'
#' @examples
#' # fit a Clayton-Weibull model
#' data(DRS)
#' clayton_wb <- rc_par_copula(data = DRS, var_list = "treat",
#'                             copula = "Clayton",
#'                             m.dist = "Weibull", method = "Newton",
#'                             iter = 500, stepsize = 1e-06)
#' summary(clayton_wb)



rc_par_copula <- function(data, var_list, copula="Clayton", m.dist="Weibull", n.cons = 4, method = "Newton", iter=500, stepsize=1e-6){


  # first screen the inputs: copula, m.dist, method #
  if (!is.data.frame(data)) stop('data must be a data frame')
  if ((!"id" %in% colnames(data)) | (!"ind" %in% colnames(data)) | (!"obs_time" %in% colnames(data)) | (!"status" %in% colnames(data))) stop('data must have id, ind, obs_time and status')
  if (!copula %in% c("Clayton","Gumbel","Copula2","Frank","Joe","AMH"))	stop('copula must be one of "Clayton","Gumbel","Copula2","Frank","Joe","AMH"')
  if (!m.dist %in% c("Weibull","Loglogistic","Gompertz","Piecewise"))	stop('m.dist must be one of "Weibull","Loglogistic","Gompertz","Piecewise"')
  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')


  # data pre-processing #
  data_2 <- data_preprocess_rc(data, var_list)
  indata1 <- data_2$indata1
  indata2 <- data_2$indata2

  tmp1 <- get_covariates_rc(indata1, var_list)
  tmp2 <- get_covariates_rc(indata2, var_list)
  x1 <- as.matrix(tmp1$x,nrow = data_2$n)
  x2 <- as.matrix(tmp2$x,nrow = data_2$n)
  var_list <- tmp1$var_list # new var_list after creating dummy variables
  p <- dim(x1)[2]
  x <- data.frame(id=c(indata1$id,indata2$id), obs_time = c(indata1$obs_time,indata2$obs_time), status = c(indata1$status,indata2$status), rbind(x1,x2))
  quantiles = NULL # for piecewise only


  ###################################
  ############ Step 1a ##############
  ###################################

  if (m.dist == "Weibull") {

    M <- survreg(as.formula(paste0("Surv(obs_time,status)~",paste0(var_list,collapse = "+"),"+cluster(id)")), data=x, dist="weibull")
    lambda_ini <- exp(M$coef[1]) # as in wikipedia
    k_ini <- 1/M$scale # k
    beta_ini <- -1*coef(M)[-1]*k_ini

  }

  if (m.dist == "Gompertz") {

    M <- flexsurvreg(as.formula(paste0("Surv(obs_time,status)~",paste0(var_list,collapse = "+"))), data=x, dist="gompertz")
    lambda_ini <- (M$coefficients[2]) # a
    k_ini <- (M$coefficients[1]) # b
    beta_ini <- M$coefficients[3:(p+2)]

  }

  if (m.dist == "Loglogistic") {

    M <- survreg(as.formula(paste0("Surv(obs_time,status)~",paste0(var_list,collapse = "+"),"+cluster(id)")), data=x, dist="loglogistic")
    lambda_ini <- exp(M$coef[1]) # as in wikipedia
    k_ini <- 1/M$scale # k
    beta_ini <- -1*coef(M)[-1]*k_ini

  }

  if (m.dist == "Piecewise") {

    quantiles<-c(0,quantile(data[data[,"status"]==1,"obs_time"],seq(1,(n.cons-1),1)/n.cons),max(data[,"obs_time"]))
    # starting values approximated by Weibull
    M <- survreg(as.formula(paste0("Surv(obs_time,status)~",paste0(var_list,collapse = "+"),"+cluster(id)")), data=x, dist="weibull", scale=2)
    lambda_ini<- rep(exp(-summary(M)$coefficients[1]),n.cons)
    beta_ini<- -summary(M)$coefficients[2:(length(var_list)+1)]

  }

  ###################################
  ############ Step 1b ##############
  ###################################
  if (copula == "AMH") {
    eta_ini<- 0.5
  }

  else if (copula == "Copula2") {
    eta_ini <- c(0.5, 1)
  }

  else {
    eta_ini <- 2
  }

  if (m.dist == "Piecewise") {
    fit0<-nlm(rc_copula_log_lik_eta, p=eta_ini, p2 = c(lambda_ini,beta_ini),
              quantiles=quantiles,x1=x1, x2=x2,indata1=indata1,indata2=indata2,
              iterlim=500, steptol = 1e-6, copula = copula, m.dist = m.dist)
  }

  if (m.dist != "Piecewise") {
    fit0<-nlm(rc_copula_log_lik_eta, p=eta_ini, p2 = c(lambda_ini,k_ini,beta_ini),
              x1=x1, x2=x2,indata1=indata1,indata2=indata2,
              iterlim=500, steptol = 1e-6, copula = copula, m.dist = m.dist)
  }

  # step 1b estimates
  eta_ini<-fit0$estimate

  # AMH shall be between 0 and 1
  if (copula == "AMH" & eta_ini[1] > 1) {eta_ini = 0.5}


  ###################################
  ############ Step 2 ###############
  ###################################

  if (method == "Newton") {
    if (m.dist != "Piecewise") {
      model_step2 <- nlm(rc_copula_log_lik, p=c(lambda_ini,k_ini,beta_ini,eta_ini),
                         x1=x1, x2=x2,indata1=indata1,indata2=indata2,
                         iterlim = iter, steptol = stepsize, hessian = T,
                         copula = copula, m.dist = m.dist)
      inv_info = pseudoinverse(model_step2$hessian) # used in SunJG 2017
      se = sqrt(diag(inv_info))
      beta = model_step2$estimate # contains lambda, k, beta and eta
      llk = -1 * model_step2$minimum
      AIC = 2*length(beta) - 2*llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat,1,lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(tmp_name1, var_list, tmp_name2)

      colnames(summary) = c("estimate","SE","stat","pvalue")
      code = model_step2$code
      output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m.dist=m.dist,
                     indata1=indata1, indata2=indata2, var_list=var_list, estimates=model_step2$estimate, x1=x1, x2=x2, inv_info = inv_info)
    }

    else if (m.dist == "Piecewise") {
      model_step2 <- nlm(rc_copula_log_lik, p=c(lambda_ini,eta_ini,beta_ini),
                         quantiles=quantiles, x1=x1, x2=x2,indata1=indata1,indata2=indata2,
                         iterlim = iter, steptol = stepsize, hessian = T,
                         copula = copula, m.dist = m.dist)
      inv_info = pseudoinverse(model_step2$hessian) # used in SunJG 2017
      se = sqrt(diag(inv_info))
      beta = model_step2$estimate # contains lambda, k, beta and eta
      llk = -1 * model_step2$minimum
      AIC = 2*length(beta) - 2*llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat,1,lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)
      rownames(summary) = if (copula != "Copula2") c(paste0("pc",1:n.cons),"eta",var_list) else c(paste0("pc",1:n.cons),"alpha","kappa",var_list)
      colnames(summary) = c("estimate","SE","stat","pvalue")
      code = model_step2$code
      output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m.dist=m.dist,
                     indata1=indata1, indata2=indata2, var_list=var_list, estimates=model_step2$estimate, x1=x1, x2=x2,
                     n.cons=n.cons, quantiles=quantiles, inv_info = inv_info)
    }
  }


  if (method != "Newton") {
    if (m.dist != "Piecewise") {
      model_step2 <- optim(c(lambda_ini,k_ini,beta_ini,eta_ini), rc_copula_log_lik,
                         x1=x1, x2=x2,indata1=indata1,indata2=indata2,
                         hessian = T, method = method,  control = list(maxit = 1000),
                         copula = copula, m.dist = m.dist)
      inv_info = pseudoinverse(model_step2$hessian) # used in SunJG 2017
      se = sqrt(diag(inv_info))
      beta = model_step2$par # contains lambda, k, beta and eta
      llk = -1 * model_step2$value
      AIC = 2*length(beta) - 2*llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat,1,lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(tmp_name1, var_list, tmp_name2)

      colnames(summary) = c("estimate","SE","stat","pvalue")
      code = model_step2$convergence
      output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m.dist=m.dist,
                     indata1=indata1, indata2=indata2, var_list=var_list, estimates=model_step2$estimate, x1=x1, x2=x2, inv_info = inv_info)
    }

    else if (m.dist == "Piecewise") {
      model_step2 <- optim(c(lambda_ini,eta_ini,beta_ini),rc_copula_log_lik,
                           quantiles=quantiles,x1=x1, x2=x2,indata1=indata1,indata2=indata2,
                           hessian = T, method = method,  control = list(maxit = 1000),
                           copula = copula, m.dist = m.dist)

      inv_info = pseudoinverse(model_step2$hessian) # used in SunJG 2017
      se = sqrt(diag(inv_info))
      beta = model_step2$par # contains lambda, k, beta and eta
      llk = -1 * model_step2$value
      AIC = 2*length(beta) - 2*llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat,1,lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)
      rownames(summary) = if (copula != "Copula2") c(paste0("pc",1:n.cons),"eta",var_list) else c(paste0("pc",1:n.cons),"alpha","kappa",var_list)
      colnames(summary) = c("estimate","SE","stat","pvalue")
      code = model_step2$convergence
      output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m.dist=m.dist,
                     indata1=indata1, indata2=indata2, var_list=var_list, estimates=model_step2$par, x1=x1, x2=x2,
                     n.cons=n.cons, quantiles=quantiles,inv_info = inv_info)
    }
  }

  class(output) <- "CopulaCenR"
  return(output)

}

