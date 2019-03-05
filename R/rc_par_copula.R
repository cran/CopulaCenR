#' Copula regression models with parametric margins for bivariate right-censored data
#'
#' @description Fits a copula model with parametric margins for bivariate right-censored data.
#'
#' @name rc_par_copula
#' @aliases rc_par_copula
#' @param data a data frame; must have \code{id} (subject id), \code{ind} (1,2 for two margins),
#' \code{obs_time}, \code{status} (0 for right-censoring, 1 for event).
#' @param var_list the list of covariates to be fitted into the model.
#' @param copula specify the copula family.
#' @param m.dist specify the marginal baseline distribution.
#' @param n.cons number of pieces, only for \code{m.dist = "Piecewise"}.
#' Default is 4.
#' @param method optimization method (see \code{?optim}); default is \code{"BFGS"};
#' also can be \code{"Newton"} (see \code{?nlm}).
#' @param iter number of iterations when \code{method = "Newton"};
#' default is 500.
#' @param stepsize size of optimization step when \code{method = "Newton"};
#' default is 1e-6.
#' @param control a list of control parameters for methods other than \code{"Newton"};
#' see \code{?optim}.
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
#' Copula-based Score Test for Bivariate Time-to-event Data,
#' with Application to a Genetic Study of AMD Progression.
#' \emph{Lifetime Data Analysis} doi:10.1007/s10985-018-09459-5. \cr
#' Tao Sun and Ying Ding (2019).
#' Copula-based Semiparametric Transformation Model for Bivariate Data
#' Under General Interval Censoring.
#' http://arxiv.org/abs/1901.01918.
#' @export
#'
#' @details The input data must be a data frame with columns \code{id} (subject id),
#' \code{ind} (1,2 for two margins; each id must have both \code{ind = 1 and 2}),
#' \code{obs_time}, \code{status} (0 for right-censoring, 1 for event)
#' and \code{covariates}. \cr
#'
#'
#' The supported copula models are \code{"Clayton"}, \code{"Gumbel"}, \code{"Frank"},
#' \code{"AMH"}, \code{"Joe"} and \code{"Copula2"}.
#' The \code{"Copula2"} model is a two-parameter copula model that incorporates
#' \code{Clayton} and \code{Gumbel} as special cases.
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
#'
#' The supported marginal distributions are \code{"Weibull"} (proportional hazards),
#' \code{"Gompertz"} (proportional hazards), \code{"Piecewise"} (proportional hazards)
#' and \code{"Loglogistic"} (proportional odds).
#' These marginal distributions follow the standard parameterization in Wikipedia.
#' We also assume the same baseline parameters between two margins. \cr
#'
#' The Weibull (PH) survival distribution is \deqn{\exp \{-(t/\lambda)^k  e^{Z^{\top}\beta}\},}
#' with \eqn{\lambda > 0} as scale and \eqn{k > 0} as shape.
#'
#' The Gompertz (PH) survival distribution is \deqn{\exp \{-\frac{b}{a}(e^{at}-1) e^{Z^{\top}\beta}\},}
#' with \eqn{a > 0} as shape and \eqn{b > 0} as rate.
#'
#' The Piecewise constant (PH) survival distribution is \deqn{\exp \{-\int_{0}^{t} \lambda_{0}(s) ds \ e^{Z^{\top}\beta}\},}
#' with \eqn{ \lambda_{0}(s) = \rho_{k}} for \eqn{s \in (a_{k-1},a_k]}.
#'
#' The Loglogistic (PO) survival distribution is \deqn{\{1+(t/\lambda)^{k} e^{Z^{\top}\beta} \}^{-1},}
#' with \eqn{\lambda > 0} as scale and \eqn{k > 0} as shape. \cr
#'
#'
#' Optimization methods can be all methods (except \code{"Brent"}) from \code{optim},
#' such as \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"}.
#' Users can also use \code{"Newton"} (from \code{nlm}).
#'
#' @return a \code{CopulaCenR} object summarizing the model.
#' Can be used as an input to general \code{S3} methods including
#' \code{summary}, \code{print}, \code{plot}, \code{lines},
#' \code{coef}, \code{logLik], \code{AIC},
#' \code{BIC}, \code{fitted}, \code{predict}.
#'
#' @examples
#' # fit a Clayton-Weibull model
#' data(DRS)
#' clayton_wb <- rc_par_copula(data = DRS, var_list = "treat",
#'                             copula = "Clayton",
#'                             m.dist = "Weibull")
#' summary(clayton_wb)



rc_par_copula <- function(data, var_list, copula="Clayton", m.dist="Weibull",
                          n.cons = 4, method = "BFGS", iter = 500, stepsize = 1e-6,
                          control = list()) {


  # first screen the inputs: copula, m.dist, method #
  if (!is.data.frame(data)) {
    stop('data must be a data frame')
  }

  if ((!"id" %in% colnames(data)) |
      (!"ind" %in% colnames(data)) |
      (!"obs_time" %in% colnames(data)) |
      (!"status" %in% colnames(data))) {
    stop('data must have id, ind, obs_time and status')
  }

  if (!copula %in% c("Clayton","Gumbel","Copula2","Frank","Joe","AMH"))	{
    stop('copula must be one of "Clayton","Gumbel","Copula2","Frank","Joe","AMH"')
  }

  if (!m.dist %in% c("Weibull","Loglogistic","Gompertz","Piecewise"))	{
    stop('m.dist must be one of "Weibull","Loglogistic","Gompertz","Piecewise"')
  }

  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	{
    stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')
  }

  # data pre-processing #
  data_2 <- data_preprocess_rc(data, var_list)
  indata1 <- data_2$indata1
  indata2 <- data_2$indata2

  tmp1 <- get_covariates_rc(indata1, var_list)
  tmp2 <- get_covariates_rc(indata2, var_list)
  x1 <- as.matrix(tmp1$x, nrow = data_2$n)
  x2 <- as.matrix(tmp2$x, nrow = data_2$n)
  var_list <- tmp1$var_list # new var_list after creating dummy variables
  p <- dim(x1)[2]
  x <- data.frame(id = c(indata1$id, indata2$id),
                  obs_time = c(indata1$obs_time, indata2$obs_time),
                  status = c(indata1$status, indata2$status),
                  rbind(x1, x2))
  quantiles <- NULL # for piecewise only


  ###################################
  ############ Step 1a ##############
  ###################################

  if (m.dist == "Weibull") {

    M <- survreg(as.formula(paste0("Surv(obs_time,status)~",
                                   paste0(var_list,collapse = "+"),
                                   "+cluster(id)")),
                 data = x, dist = "weibull")
    lambda_ini <- exp(M$coef[1]) # as in wikipedia
    k_ini <- 1/M$scale # k
    beta_ini <- -1 * coef(M)[-1] * k_ini

  }

  if (m.dist == "Gompertz") {

    M <- flexsurvreg(as.formula(paste0("Surv(obs_time,status)~",
                                       paste0(var_list, collapse = "+"))),
                     data = x, dist = "gompertz")
    lambda_ini <- (M$coefficients[2]) # a
    k_ini <- (M$coefficients[1]) # b
    beta_ini <- M$coefficients[3:(p+2)]

  }

  if (m.dist == "Loglogistic") {

    M <- survreg(as.formula(paste0("Surv(obs_time,status)~",
                                   paste0(var_list, collapse = "+"),
                                   "+cluster(id)")),
                 data = x, dist = "loglogistic")
    lambda_ini <- exp(M$coef[1]) # as in wikipedia
    k_ini <- 1/M$scale # k
    beta_ini <- -1 * coef(M)[-1] * k_ini

  }

  if (m.dist == "Piecewise") {

    quantiles <- c(0, quantile(data[data[,"status"]==1,"obs_time"],
                               seq(1, (n.cons-1), 1)/n.cons),
                   max(data[,"obs_time"])
                   )
    # starting values approximated by Weibull
    M <- survreg(as.formula(paste0("Surv(obs_time,status)~",
                                   paste0(var_list, collapse = "+"),
                                   "+cluster(id)")),
                 data = x, dist = "weibull", scale = 2)
    lambda_ini<- rep(exp(-summary(M)$coefficients[1]), n.cons)
    beta_ini <- -summary(M)$coefficients[2:(length(var_list) + 1)]

  }

  ###################################
  ############ Step 1b ##############
  ###################################
  if (copula == "AMH") {
    eta_ini<- 0.5
  }

  else if (copula == "Copula2") {
    eta_ini <- c(1, 1)
  }

  else {
    eta_ini <- 2
  }

  if (method == "Newton") {

    if (m.dist == "Piecewise") {
      fit0 <- nlm(rc_copula_log_lik_eta, p=eta_ini, p2 = c(lambda_ini,beta_ini),
                quantiles = quantiles, x1 = x1, x2 = x2, indata1 = indata1,
                indata2 = indata2, iterlim = iter, steptol = stepsize,
                copula = copula, m.dist = m.dist)
    } else {
      fit0 <- nlm(rc_copula_log_lik_eta, p = eta_ini,
                  p2 = c(lambda_ini,k_ini,beta_ini),
                  x1 = x1, x2 = x2, indata1 = indata1,indata2 = indata2,
                  iterlim = iter, steptol = stepsize, copula = copula,
                  m.dist = m.dist)
    }

    eta_ini <- fit0$estimate

  } else {

    if (m.dist == "Piecewise") {
      fit0 <- optim(par = eta_ini, rc_copula_log_lik_eta,
                    p2 = c(lambda_ini, beta_ini),
                    method = method,  control = control, hessian = FALSE,
                    quantiles = quantiles,x1 = x1, x2 = x2,indata1 = indata1,
                    indata2 = indata2, copula = copula, m.dist = m.dist)
    } else {
      fit0 <- optim(par = eta_ini, rc_copula_log_lik_eta,
                    p2 = c(lambda_ini,k_ini,beta_ini),
                    method = method,  control = control, hessian = FALSE,
                    x1 = x1, x2 = x2,indata1 = indata1,indata2 = indata2,
                    copula = copula, m.dist = m.dist)
    }
    eta_ini <- fit0$par
  }



  # AMH shall be between 0 and 1
  if (copula == "AMH" & eta_ini[1] > 1) {eta_ini <- 0.5}


  ###################################
  ############ Step 2 ###############
  ###################################

  if (method == "Newton") {
    if (m.dist != "Piecewise") {
      model_step2 <- nlm(rc_copula_log_lik, p = c(lambda_ini,k_ini,beta_ini,eta_ini),
                         x1 = x1, x2 = x2,indata1 = indata1,indata2 = indata2,
                         iterlim = iter, steptol = stepsize, hessian = T,
                         copula = copula, m.dist = m.dist)
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      se <- sqrt(dih)
      beta <- model_step2$estimate # contains lambda, k, beta and eta
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(beta) - 2 * llk
      stat <- (beta - 0)^2/se^2
      pvalue <- pchisq(stat,1,lower.tail=F)
      summary <- cbind(beta, se, stat, pvalue)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) = c(tmp_name1, var_list, tmp_name2)

      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$code
      output <- list(code = code, summary = summary, llk = llk, AIC = AIC,
                     copula = copula, m.dist = m.dist, indata1 = indata1,
                     indata2 = indata2, var_list = var_list,
                     estimates = model_step2$estimate, x1 = x1, x2 = x2,
                     inv_info = inv_info)
    }

    else if (m.dist == "Piecewise") {
      model_step2 <- nlm(rc_copula_log_lik, p = c(lambda_ini,eta_ini,beta_ini),
                         quantiles = quantiles, x1 = x1, x2 = x2,
                         indata1 = indata1,indata2 = indata2,
                         iterlim = iter, steptol = stepsize, hessian = T,
                         copula = copula, m.dist = m.dist)
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      se <- sqrt(dih)
      beta <- model_step2$estimate # contains lambda, k, beta and eta
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(beta) - 2 * llk
      stat <- (beta - 0)^2/se^2
      pvalue <- pchisq(stat, 1, lower.tail=F)
      summary <- cbind(beta, se, stat, pvalue)
      rownames(summary) <- if (copula != "Copula2") c(paste0("pc",1:n.cons),"eta",var_list) else c(paste0("pc",1:n.cons),"alpha","kappa",var_list)
      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$code
      output <- list(code = code, summary = summary, llk = llk, AIC = AIC,
                     copula = copula, m.dist = m.dist, indata1 = indata1,
                     indata2 = indata2, var_list = var_list,
                     estimates = model_step2$estimate, x1 = x1, x2 = x2,
                     n.cons = n.cons, quantiles = quantiles, inv_info = inv_info)
    }
  }


  if (method != "Newton") {
    if (m.dist != "Piecewise") {
      model_step2 <- optim(c(lambda_ini,k_ini,beta_ini,eta_ini), rc_copula_log_lik,
                         x1 = x1, x2 = x2,indata1 = indata1, indata2 = indata2,
                         hessian = T, method = method,  control = control,
                         copula = copula, m.dist = m.dist)
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      se <- sqrt(dih)
      beta <- model_step2$par # contains lambda, k, beta and eta
      llk <- -1 * model_step2$value
      AIC <- 2 * length(beta) - 2 * llk
      stat <- (beta - 0)^2/se^2
      pvalue <- pchisq(stat, 1, lower.tail=F)
      summary <- cbind(beta, se, stat, pvalue)

      tmp_name1 <- if (m.dist != "Gompertz") c("lambda","k") else c("a","b")
      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(tmp_name1, var_list, tmp_name2)

      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$convergence
      output <- list(code = code, summary = summary, llk = llk, AIC = AIC,
                     copula = copula, m.dist = m.dist, indata1 = indata1,
                     indata2 = indata2, var_list = var_list,
                     estimates = model_step2$par, x1 = x1, x2 = x2,
                     inv_info = inv_info)
    }

    else if (m.dist == "Piecewise") {
      model_step2 <- optim(c(lambda_ini,eta_ini,beta_ini), rc_copula_log_lik,
                           quantiles = quantiles, x1 = x1, x2 = x2,
                           indata1 = indata1, indata2 = indata2,
                           hessian = T, method = method, control = control,
                           copula = copula, m.dist = m.dist)

      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      se <- sqrt(dih)
      beta <- model_step2$par # contains lambda, k, beta and eta
      llk <- -1 * model_step2$value
      AIC <- 2 * length(beta) - 2 * llk
      stat <- (beta - 0)^2/se^2
      pvalue <- pchisq(stat, 1, lower.tail=F)
      summary <- cbind(beta, se, stat, pvalue)
      rownames(summary) <- if (copula != "Copula2") c(paste0("pc",1:n.cons),"eta",var_list) else c(paste0("pc",1:n.cons),"alpha","kappa",var_list)
      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$convergence
      output <- list(code = code, summary = summary, llk = llk, AIC = AIC,
                     copula = copula, m.dist = m.dist,
                     indata1 = indata1, indata2 = indata2, var_list = var_list,
                     estimates = model_step2$par, x1 = x1, x2 = x2,
                     n.cons = n.cons, quantiles = quantiles,inv_info = inv_info)
    }
  }

  class(output) <- "CopulaCenR"
  return(output)

}

