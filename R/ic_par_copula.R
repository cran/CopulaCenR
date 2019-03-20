#' Copula regression models with parametric margins for bivariate interval-censored data
#'
#' @description Fits a copula model with parametric margins for bivariate interval-censored data.
#'
#' @name ic_par_copula
#' @aliases ic_par_copula
#' @param data a data frame; must have \code{id} (subject id), \code{ind} (1,2 for two units in each subject),
#' \code{Left} (0 if left-censoring), \code{Right} (Inf if right-censoring), \code{status} (0 for right-censoring,
#' 1 for interval-censoring or left-censoring), and \code{covariates} by column.
#' @param var_list the list of covariates to be fitted into the copula model.
#' @param copula Types of copula model.
#' @param m.dist baseline marginal distribution.
#' @param method optimization method (see ?optim); default is "BFGS";
#' also can be "Newton" (see ?nlm).
#' @param iter number of iterations when method is \code{"Newton"};
#' default is 300.
#' @param stepsize size of optimization step when method is \code{"Newton"};
#' default is 1e-5.
#' @param control a list of control parameters for methods other than \code{"Newton"};
#' see \code{?optim}.
#' @param hes default is \code{TRUE} for hessian calculation.
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
#' Copula-based Score Test for Bivariate Time-to-event Data,
#' with Application to a Genetic Study of AMD Progression.
#' \emph{Lifetime Data Analysis} doi:10.1007/s10985-018-09459-5. \cr
#' Tao Sun and Ying Ding (2019).
#' Copula-based Semiparametric Transformation Model for Bivariate Data
#' Under General Interval Censoring.
#' http://arxiv.org/abs/1901.01918.
#'
#' @details The input data must be a data frame. with columns \code{id} (sample id),
#' \code{ind} (1,2 for the two units from the same id),
#' \code{Left} (0 if left-censoring), \code{Right} (Inf if right-censoring),
#' \code{status} (0 for right-censoring, 1 for interval-censoring or left-censoring),
#' and \code{covariates}. The function does not allow \code{Left} == \code{Right}. \cr
#'
#'
#' The supported copula models are \code{"Clayton"}, \code{"Gumbel"}, \code{"Frank"},
#' \code{"AMH"}, \code{"Joe"} and \code{"Copula2"}.
#' The \code{"Copula2"} model is a two-parameter copula model that incorporates \code{Clayton}
#' and \code{Gumbel} as special cases.
#' The parametric generator functions of copula functions are list below:
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
#' with \eqn{\eta \geq 1} and Kendall's \eqn{\tau = 1 - 4 \sum_{k=1}^{\infty} \frac{1}{k(\eta k+2)\{\eta(k-1)+2\}}}.
#'
#' The Two-parameter copula (Copula2) has a generator \deqn{\phi_{\eta}(t) = \{1/(1+t^{\alpha})\}^{\kappa},}
#' with \eqn{\alpha \in (0,1], \kappa > 0} and Kendall's \eqn{\tau = 1-2\alpha\kappa/(2\kappa+1)}. \cr
#'
#'
#' The supported marginal distributions are \code{"Weibull"} (proportional hazards),
#' \code{"Gompertz"} (proportional hazards) and \code{"Loglogistic"} (proportional odds).
#' These marginal distributions are listed below
#' and we assume the same baseline parameters between two margins. \cr
#'
#' The Weibull (PH) survival distribution is \deqn{\exp \{-(t/\lambda)^k  e^{Z^{\top}\beta}\},}
#' with \eqn{\lambda > 0} as scale and \eqn{k > 0} as shape.
#'
#' The Gompertz (PH) survival distribution is \deqn{\exp \{-\frac{b}{a}(e^{at}-1) e^{Z^{\top}\beta}\},}
#' with \eqn{a > 0} as shape and \eqn{b > 0} as rate.
#'
#' The Loglogistic (PO) survival distribution is \deqn{\{1+(t/\lambda)^{k} e^{Z^{\top}\beta} \}^{-1},}
#' with \eqn{\lambda > 0} as scale and \eqn{k > 0} as shape. \cr
#'
#'
#' Optimization methods can be all methods (except \code{"Brent"}) from \code{optim}, such as
#' \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"}.
#' Users can also use \code{"Newton"} (from \code{nlm}).
#'
#' @return a \code{CopulaCenR} object summarizing the model.
#' Can be used as an input to general \code{S3} methods including
#' \code{summary}, \code{print}, \code{plot}, \code{lines},
#' \code{coef}, \code{logLik}, \code{AIC},
#' \code{BIC}, \code{fitted}, \code{predict}.
#'
#' @examples
#' # fit a Copula2-Weibull model
#' data(AREDS)
#' copula2_wb <- ic_par_copula(data = AREDS, copula = "Copula2",
#'                   m.dist = "Weibull",
#'                   var_list = c("ENROLLAGE","rs2284665"))
#' summary(copula2_wb)




ic_par_copula <- function(data, var_list, copula, m.dist = "Weibull",
                          method = "BFGS", iter=300, stepsize=1e-5,
                          hes = TRUE, control = list()){

  # first screen the inputs: copula, m.dist, method #
  if (!is.data.frame(data)) {
    stop('data must be a data frame')
  }

  if ((!"id" %in% colnames(data)) |
      (!"ind" %in% colnames(data)) |
      (!"Left" %in% colnames(data)) |
      (!"Right" %in% colnames(data)) |
      (!"status" %in% colnames(data))) {
    stop('data must have id, ind, Left, Right and status')
  }

  if (!copula %in% c("Clayton","Gumbel","Copula2","Frank","Joe","AMH"))	{
    stop('copula must be one of "Clayton","Gumbel","Copula2","Frank","Joe","AMH"')
  }

  if (!m.dist %in% c("Weibull","Loglogistic","Gompertz"))	{
    stop('m.dist must be one of "Weibull","Loglogistic","Gompertz"')
  }

  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	{
    stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')
  }


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

    weibull_reg <- ic_par(cbind(Left, Right) ~ . , data = x, dist = 'weibull', model = 'ph')
    lambda_ini <- exp(weibull_reg$coefficients[2])
    k_ini <- exp(weibull_reg$coefficients[1])
    beta_ini <- coef(weibull_reg)[3:(p+2)]
    names(lambda_ini) <- NULL
    names(k_ini) <- NULL
    names(beta_ini) <- NULL

  }

  else if (m.dist == "Gompertz") {

    tmp <-  Surv(x$Left, x$Right, type = "interval2")
    weibull_reg <- flexsurvreg(as.formula(paste0("tmp~", paste0(var_list, collapse = "+"))),
                               data = x, dist = "gompertz")
    lambda_ini <-  weibull_reg$res["shape","est"] # a
    k_ini <-  weibull_reg$res["rate","est"] # b
    beta_ini <- weibull_reg$coefficients[3:(p+2)]
    names(lambda_ini) <- NULL
    names(k_ini) <- NULL
    names(beta_ini) <- NULL

  }

  else if (m.dist == "Loglogistic") {

    loglogistic_reg<-ic_par(cbind(Left, Right) ~ ., data = x, dist = 'loglogistic', model = 'po')
    lambda_ini<-exp(loglogistic_reg$coefficients[1])  # scale; equ to b in qweibull
    k_ini<-exp(loglogistic_reg$coefficients[2]) #shape; equ to a in qweibull
    beta_ini <- -coef(loglogistic_reg)[3:(p+2)]
    names(lambda_ini) <- NULL
    names(k_ini) <- NULL
    names(beta_ini) <- NULL

  }



  ###################################
  ############ Step 1b ##############
  ###################################

  if (copula == "AMH") {
    eta_ini <- 0
  }

  else if (copula == "Copula2") {
    eta_ini <- c(0, 0)
  }

  else {
    eta_ini <- 0
  }


  if (method == "Newton") {
    model_step1b <- nlm(ic_copula_log_lik_param_eta, eta_ini, hessian = FALSE,
                        lambda = lambda_ini, k = k_ini, beta = beta_ini,
                        x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right,
                        t2_left = t2_left, t2_right = t2_right, indata1 = indata1,
                        indata2 = indata2, copula = copula, m.dist = m.dist)
    eta_ini <- exp(model_step1b$estimate)

  } else {
    model_step1b <- optim(par = eta_ini, ic_copula_log_lik_param_eta, hessian = FALSE,
                          method = method, control = control,
                          lambda = lambda_ini, k = k_ini, beta = beta_ini,
                          x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right,
                          t2_left = t2_left, t2_right = t2_right, indata1 = indata1,
                          indata2 = indata2, copula = copula, m.dist = m.dist)
    eta_ini <- exp(model_step1b$par)
  }


  # Quality check for step 1b estimates
  # AMH shall be between 0 and 1
  if (copula == "AMH" & eta_ini[1] > 1) {eta_ini <- 0.5}
  # Gumbel shall be >= 1
  if (copula == "Gumbel" & eta_ini[1] < 1) {eta_ini <- 1}
  # Joe shall be >= 1
  if (copula == "Joe" & eta_ini[1] < 1) {eta_ini <- 1}
  # Copula2 alpha shall be between 0 and 1
  if (copula == "Copula2" & eta_ini[1] > 1) {eta_ini[1] <- 0.5}


  ###################################
  ############ Step 2 ###############
  ###################################
  if (method == "Newton") {

    model_step2<-nlm(ic_copula_log_lik_param, c(lambda_ini, k_ini, beta_ini, eta_ini),
                     p, x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right,
                     t2_left = t2_left, t2_right = t2_right, indata1 = indata1,
                     indata2 = indata2, hessian = hes, iterlim = iter,
                     steptol = stepsize, copula = copula, m.dist = m.dist)

    if (isTRUE(hes)) {
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      se <- sqrt(dih)
      beta <- model_step2$estimate # contains lambda, k, beta and eta
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(beta) - 2 * llk
      stat <- (beta-0)^2/se^2
      pvalue <- pchisq(stat,1,lower.tail = F)
      summary <- cbind(beta, se, stat, pvalue)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(tmp_name1, var_list, tmp_name2)

      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$code
    }

    if (!isTRUE(hes)) {

      inv_info <- NULL
      beta <- model_step2$estimate
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(beta) - 2 * llk
      summary <- cbind(beta)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(tmp_name1, var_list, tmp_name2)
      colnames(summary) <- c("estimate")
      code <- model_step2$code
    }

    output <- list(code = code, summary = summary, llk = llk, AIC = AIC, copula = copula,
                   m.dist = m.dist,  indata1 = indata1, indata2 = indata2, var_list = var_list,
                   estimates = model_step2$estimate, x1 = x1, x2 = x2,inv_info = inv_info)
  }


  if (method != "Newton") {

    model_step2 <- optim(par = c(lambda_ini, k_ini, beta_ini, eta_ini),
                         ic_copula_log_lik_param,
                         method = method, hessian = hes, control = control,
                         p = p, x1 = x1, x2 = x2, t1_left = t1_left,
                         t1_right = t1_right, t2_left = t2_left,
                         t2_right = t2_right, indata1 = indata1,
                         indata2 = indata2, copula = copula, m.dist = m.dist)

    if (isTRUE(hes)) {
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      se <- sqrt(dih)
      beta <- model_step2$par # contains lambda, k, beta and eta
      llk <- -1 * model_step2$value
      AIC <- 2 * length(beta) - 2 * llk
      stat <- (beta - 0)^2/se^2
      pvalue <- pchisq(stat,1,lower.tail = F)
      summary <- cbind(beta, se, stat, pvalue)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(tmp_name1, var_list, tmp_name2)
      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$convergence
    }

    if (!isTRUE(hes)) {

      inv_info <- NULL
      beta <- model_step2$par # contains lambda, k, beta and eta
      llk <- -1 * model_step2$value
      AIC <- 2 * length(beta) - 2 * llk
      summary <- cbind(beta)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(tmp_name1, var_list, tmp_name2)

      colnames(summary) <- c("estimate")
      code <- model_step2$convergence
    }

    output <- list(code = code, summary = summary, llk = llk, AIC = AIC, copula = copula,
                   m.dist = m.dist,  indata1 = indata1, indata2 = indata2, var_list = var_list,
                   estimates = model_step2$par, x1 = x1, x2 = x2,inv_info = inv_info)
  }

  class(output) <- "CopulaCenR"
  return(output)
}





