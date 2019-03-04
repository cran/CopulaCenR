#' Copula regression models with parametric margins for bivariate interval-censored data
#'
#' @description Fits a copula model with parametric margins for bivariate interval-censored data.
#'
#' @name ic_par_copula
#' @aliases ic_par_copula
#' @param data a data frame; must have id (subject id), ind (1,2 for two units in each subject),
#' Left (0 if left-censoring), Right (Inf if right-censoring), status (0 for right-censoring, 1 for interval-censoring or left-censoring), and covariates by column.
#' @param var_list the list of covariates to be fitted into the copula model.
#' @param copula Types of copula model.
#' @param m.dist baseline marginal distribution.
#' @param method optimization method; default is "Newton" (from nlm).
#' @param iter number of iterations; default is 300.
#' @param stepsize size of optimization step; default is 1e-5.
#' @param hes default = T for hessian calculation.
#' @importFrom corpcor pseudoinverse
#' @importFrom stats pchisq
#' @importFrom stats coef
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival cluster
#' @importFrom flexsurv flexsurvreg
#' @importFrom icenReg ic_par
#' @export
#'
#' @source
#' Tao Sun, Yi Liu, Richard J. Cook, Wei Chen and Ying Ding (2018).
#' Copula-based Score Test for Bivariate Time-to-event Data, with Application to a Genetic Study of AMD Progression.
#' \emph{Lifetime Data Analysis} doi:10.1007/s10985-018-09459-5. \cr
#' Tao Sun and Ying Ding (2019).
#' Copula-based Semiparametric Transformation Model for Bivariate Data Under General Interval Censoring.
#' http://arxiv.org/abs/1901.01918.
#'
#' @details The supported copula models are "Clayton", "Gumbel", "Frank", "AMH", "Joe", "Copula2".
#' The "Copula2" model is a two-parameter copula model that incorporates Clayton and Gumbel as special cases.
#' The supported marginal distributions are "Weibull" (proportional hazards), "Gompertz" (proportional hazards) and "Loglogistic" (proportional odds).
#' We assume the same baseline parameters between two margins. \cr
#'
#' The input data must be a data frame. with columns id (sample id), ind (1,2 for the two units from the same id),
#' Left (0 if left-censoring), Right (Inf if right-censoring), status (0 for right-censoring, 1 for interval-censoring or left-censoring), and covariates. Do not allow Left == Right. \cr
#'
#' Optimization methods can be "Newton" (from nlm) and all methods from optim, such as "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN". "Brent" is not allowed due to more than one estimators.
#'
#' @return a CopulaCenR object summarizing the model. Can be used as an input to general S3 methods including
#' summary, plot, lines, coef, logLik, predict.
#'
#' @examples
#' # fit a Copula2-Loglogistic model
#' data(AREDS)
#' copula2_loglog <- ic_par_copula(data = AREDS, copula = "Copula2",
#'                   m.dist = "Loglogistic",
#'                   var_list = c("ENROLLAGE","rs2284665","SevScaleBL"),
#'                   method = "Newton", iter = 300, stepsize = 1e-6)
#' summary(copula2_loglog)




ic_par_copula <- function(data, var_list, copula, m.dist = "Loglogistic", method = "Newton", iter=300, stepsize=1e-5, hes = TRUE){

  # first screen the inputs: copula, m.dist, method #
  if (!is.data.frame(data)) stop('data must be a data frame')
  if ((!"id" %in% colnames(data)) | (!"ind" %in% colnames(data)) | (!"Left" %in% colnames(data)) | (!"Right" %in% colnames(data)) | (!"status" %in% colnames(data))) stop('data must have id, ind, Left, Right and status')
  if (!copula %in% c("Clayton","Gumbel","Copula2","Frank","Joe","AMH"))	stop('copula must be one of "Clayton","Gumbel","Copula2","Frank","Joe","AMH"')
  if (!m.dist %in% c("Weibull","Loglogistic","Gompertz"))	stop('m.dist must be one of "Weibull","Loglogistic","Gompertz"')
  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')

  # data pre-processing #
  data_processed <- data_process(data, var_list)
  indata1 <- data_processed$indata1
  indata2 <- data_processed$indata2
  t1_left <- data_processed$t1_left
  t1_right <- data_processed$t1_right
  t2_left <- data_processed$t2_left
  t2_right <- data_processed$t2_right
  n <- data_processed$n
  p <- data_processed$p
  x1 <- data_processed$x1
  x2 <- data_processed$x2
  x <- data_processed$x
  var_list <- data_processed$var_list



  ###################################
  ############ Step 1a ##############
  ###################################

  if (m.dist == "Weibull") {

    weibull_reg<-ic_par(cbind(Left, Right) ~ . , data = x, dist = 'weibull', model = 'ph')
    lambda_ini <-exp(weibull_reg$coefficients[2]) # scale; equ to b in qweibull
    k_ini <-exp(weibull_reg$coefficients[1]) #shape; equ to a in qweibull
    beta_ini <- coef(weibull_reg)[3:(p+2)]
    names(lambda_ini) = NULL
    names(k_ini) = NULL
    names(beta_ini) = NULL

  }

  if (m.dist == "Gompertz") {

    tmp <-  Surv(x$Left, x$Right, type = "interval2")
    weibull_reg <- flexsurvreg(as.formula(paste0("tmp~",paste0(var_list,collapse = "+"))), data=x, dist="gompertz")
    lambda_ini <-  weibull_reg$res["shape","est"] # a
    k_ini <-  weibull_reg$res["rate","est"] # b
    beta_ini <- weibull_reg$coefficients[3:(p+2)]
    names(lambda_ini) = NULL
    names(k_ini) = NULL
    names(beta_ini) = NULL

  }

  if (m.dist == "Loglogistic") {

    loglogistic_reg<-ic_par(cbind(Left, Right) ~ ., data = x, dist = 'loglogistic', model = 'po')
    lambda_ini<-exp(loglogistic_reg$coefficients[1])  # scale; equ to b in qweibull
    k_ini<-exp(loglogistic_reg$coefficients[2]) #shape; equ to a in qweibull
    beta_ini <- -coef(loglogistic_reg)[3:(p+2)]
    names(lambda_ini) = NULL
    names(k_ini) = NULL
    names(beta_ini) = NULL

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



  model_step1b <- nlm(ic_copula_log_lik_param_eta, eta_ini, hessian = FALSE,
                    lambda=lambda_ini, k=k_ini, beta=beta_ini,
                    x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                    copula = copula, m.dist = m.dist)


  # step 1b estimates
  eta_ini<-model_step1b$estimate

  # AMH shall be between 0 and 1
  if (copula == "AMH" & eta_ini[1] > 1) {eta_ini = 0.5}


  ###################################
  ############ Step 2 ###############
  ###################################
  if (method == "Newton") {

    model_step2<-nlm(ic_copula_log_lik_param, c(lambda_ini,k_ini,beta_ini,eta_ini),
                     p, x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                     hessian = hes, iterlim = iter ,steptol = stepsize,
                     copula = copula, m.dist = m.dist)

    if (isTRUE(hes)) {
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
    }

    if (!isTRUE(hes)) {

      inv_info = NULL
      beta = model_step2$estimate # contains lambda, k, beta and eta
      llk = -1 * model_step2$minimum
      AIC = 2*length(beta) - 2*llk
      summary = cbind(beta)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(tmp_name1, var_list, tmp_name2)

      colnames(summary) = c("estimate")
      code = model_step2$code
    }

    output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m.dist=m.dist,  indata1=indata1, indata2=indata2, var_list=var_list, estimates=model_step2$estimate, x1=x1, x2=x2,inv_info = inv_info)
  }


  if (method != "Newton") {

    model_step2 <- optim(par=c(lambda_ini,k_ini,beta_ini,eta_ini), ic_copula_log_lik_param,
                         method = method, hessian = hes, control=list(type=1),
                         p=p, x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                         copula = copula, m.dist = m.dist)


    if (isTRUE(hes)) {
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
    }

    if (!isTRUE(hes)) {

      inv_info = NULL
      beta = model_step2$par # contains lambda, k, beta and eta
      llk = -1 * model_step2$value
      AIC = 2*length(beta) - 2*llk
      summary = cbind(beta)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(tmp_name1, var_list, tmp_name2)

      colnames(summary) = c("estimate")
      code = model_step2$convergence
    }

    output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m.dist=m.dist,  indata1=indata1, indata2=indata2, var_list=var_list, estimates=model_step2$par, x1=x1, x2=x2,inv_info = inv_info)
  }

  class(output) <- "CopulaCenR"
  return(output)
}





