#' Copula regression models with semiparametric margins for bivariate interval-censored data
#'
#' @description Fits a copula model with semiparametric margins for bivariate interval-censored data.
#'
#' @name ic_sp_copula
#' @aliases ic_sp_copula
#' @param data a data frame; must have id (subject id), ind (1,2 for two units in each subject),
#' Left (0 if left-censoring), Right (Inf if right-censoring), status (0 for right-censoring, 1 for interval-censoring or left-censoring), and covariates by column.
#' @param var_list the list of covariates to be fitted into the copula model.
#' @param copula Types of copula model.
#' @param r postive transformation parameter for the semi-parametric linear transformation model.
#' @param m integer, degree of Berstein polynomials (assume same phi for two margins, leading to m+1 Berstein coefficients parameters).
#' @param l the left bound for all Left and Right endpoints of observed finite intervals; default is 0.
#' @param u the right bound for all Left and Right endpoints of observed finite intervalsl; has to be a finite value
#' @param method optimization method; default is "Newton" using nlm; Alternatives are methods in the general optim function, such as "CG" (Conjugate Gradient).
#' @param iter number of iterations allowed for method = "Newton"; default is 300.
#' @param stepsize size of optimization step for method = "Newton"; default is 1e-5.
#' @param hes default = T for hessian calculation; if LRT is desired, can set hes = FALSE to save time.
#' @importFrom corpcor pseudoinverse
#' @importFrom stats pchisq
#' @importFrom stats optim
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
#' The marginal semiparametric transformation models are built based on Bernstein polynomials.
#' The parameter r (a positive number) decides the model assumption, with r = 1 for proportional hazards and 3 for proportional odds.
#' The parameter m (a positive integer) is the degree of Bernstein polynomials, with default m = 4.
#' Usually, m and r are selected based on the AIC value, which is provided after fitting the copula model. \cr
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
#' # fit a Copula2-Semiparametric model
#' data(AREDS)
#' copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2", l = 0, u = 15,
#'               m = 4, r = 3, var_list = c("ENROLLAGE","rs2284665","SevScaleBL"),
#'               iter = 300, stepsize = 1e-6, method = "Newton")
#'summary(copula2_sp)



ic_sp_copula <- function(data, var_list, l=0, u, copula = "Copula2", m = 3, r = 3, method = "Newton", iter=300, stepsize=1e-5, hes = TRUE){

  # first screen the inputs: copula, m.dist, method #
  if (!is.data.frame(data)) stop('data must be a data frame')
  if ((!"id" %in% colnames(data)) | (!"ind" %in% colnames(data)) | (!"Left" %in% colnames(data)) | (!"Right" %in% colnames(data)) | (!"status" %in% colnames(data))) stop('data must have id, ind, Left, Right and status')
  if (!copula %in% c("Clayton","Gumbel","Copula2","Frank","Joe","AMH"))	stop('copula must be one of "Clayton","Gumbel","Copula2","Frank","Joe","AMH"')
  if (l<0 | u<=0) stop('l must be >= 0 and u greater than non-infinite values of Left and Right')
  if (r <= 0) stop('r must be a positive number')
  if (m != round(m)) stop('m must be a positive integer')
  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')

  # data pre-processing #
  data_processed <- data_process_sieve(data, l, u, var_list, m)
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
  var_list <- data_processed$var_list

  # BP
  bl1 <- data_processed$bl1
  br1 <- data_processed$br1
  bl2 <- data_processed$bl2
  br2 <- data_processed$br2


  ###################################
  ############ Step 1a ##############
  ###################################
  model_step1a <- nlm(estimate_sieve_step1a,rep(0,(p+m+1)),hessian = FALSE,iterlim = 200,
                      p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r) # p beta, m+1 polynomails
  beta <- model_step1a$estimate[1:p]
  phi <- model_step1a$estimate[(p+1):(p+1+m)] # 4-7
  ep<-cumsum(exp(phi))



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


  if (method == "Newton") {
    model_step1b <- nlm(ic_copula_log_lik_sieve_eta,eta_ini,hessian = FALSE,
                        beta=beta, ep=ep, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                        copula = copula, print.level = 0)
    eta_ini<-model_step1b$estimate
  }


  if (method != "Newton") {
    model_step1b <- optim(par=eta_ini, ic_copula_log_lik_sieve_eta,
                          method = "CG", control=list(type=1), hessian = FALSE,
                         beta=beta, ep=ep, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                         copula = copula)
    eta_ini<-model_step1b$par
  }

  # AMH shall be between 0 and 1
  if (copula == "AMH" & eta_ini[1] > 1) {eta_ini = 0.5}

  ###################################
  ############ Step 2 ###############
  ###################################

  if (method == "Newton") {
    model_step2 <- nlm(ic_copula_log_lik_sieve, c(model_step1a$estimate,eta_ini),
                       hessian = hes, iterlim = iter ,steptol = stepsize,
                       p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                       copula = copula, print.level = 0)

    if (isTRUE(hes)) {
      inv_info = pseudoinverse(model_step2$hessian) # used in SunJG 2017
      dih = diag(inv_info)
      dih[dih<0] = 0
      dih = sqrt(dih)
      # dih = sqrt(diag(inv_info))
      se <- if (copula != "Copula2") dih[c(1:p,length(dih))] else dih[c(1:p,length(dih)-1,length(dih))]
      beta <- if (copula != "Copula2") model_step2$estimate[c(1:p,length(dih))] else model_step2$estimate[c(1:p,length(dih)-1,length(dih))]
      llk = -1 * model_step2$minimum
      AIC = 2*length(model_step2$estimate) - 2*llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat,1,lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(var_list, tmp_name2)

      colnames(summary) = c("estimate","SE","stat","pvalue")
      code = model_step2$code
    }

    if (!isTRUE(hes)) {
      inv_info = NULL
      beta <- if (copula != "Copula2")  model_step2$estimate[c(1:p,length(model_step2$estimate))] else model_step2$estimate[c(1:p,length(model_step2$estimate)-1,length(model_step2$estimate))]
      llk = -1 * model_step2$minimum
      AIC = 2*length(model_step2$estimate) - 2*llk
      summary = cbind(beta)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(var_list, tmp_name2)

      colnames(summary) = c("estimate")
      code = model_step2$code
    }
    output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m=m, r=r, indata1=indata1, indata2=indata2, var_list=var_list, l=l, u=u, bl1=bl1, br1=br1, bl2=bl2, br2=br2, estimates=model_step2$estimate, x1=x1, x2=x2,inv_info = inv_info)
  }



  if (method != "Newton") {

    model_step2 <- optim(par=c(model_step1a$estimate,eta_ini), ic_copula_log_lik_sieve,
                        method = method, hessian = hes, control=list(type=1),
                        p=p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                        copula = copula)

    if (isTRUE(hes)) {
      inv_info = pseudoinverse(model_step2$hessian) # used in SunJG 2017
      dih = diag(inv_info)
      dih[dih<0] = 0
      dih = sqrt(dih)
      # dih = sqrt(diag(inv_info))
      se <- if (copula != "Copula2") dih[c(1:p,length(dih))] else dih[c(1:p,length(dih)-1,length(dih))]
      beta <- if (copula != "Copula2") model_step2$par[c(1:p,length(dih))] else model_step2$par[c(1:p,length(dih)-1,length(dih))]
      llk = -1 * model_step2$value
      AIC = 2*length(model_step2$par) - 2*llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat,1,lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(var_list, tmp_name2)

      colnames(summary) = c("estimate","SE","stat","pvalue")
      code = model_step2$convergence
    }

    if (!isTRUE(hes)) {
      inv_info = NULL
      beta <- if (copula != "Copula2")  model_step2$par[c(1:p,length(model_step2$par))] else model_step2$par[c(1:p,length(model_step2$par)-1,length(model_step2$par))]
      llk = -1 * model_step2$value
      AIC = 2*length(model_step2$par) - 2*llk
      summary = cbind(beta)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(var_list, tmp_name2)

      colnames(summary) = c("estimate")
      code = model_step2$convergence
    }
    output <- list(code=code, summary=summary, llk=llk, AIC=AIC, copula=copula, m=m, r=r, indata1=indata1, indata2=indata2, var_list=var_list, l=l, u=u, bl1=bl1, br1=br1, bl2=bl2, br2=br2, estimates=model_step2$par, x1=x1, x2=x2, inv_info = inv_info)
  }

  class(output) <- "CopulaCenR"
  return(output)
}

