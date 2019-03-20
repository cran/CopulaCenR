#' Copula regression models with semiparametric margins for bivariate interval-censored data
#'
#' @description Fits a copula model with semiparametric margins for bivariate interval-censored data.
#'
#' @name ic_spTran_copula
#' @aliases ic_spTran_copula
#' @param data a data frame; must have \code{id} (subject id), \code{ind} (1,2 for two units in each subject),
#' \code{Left} (0 if left-censoring), \code{Right} (Inf if right-censoring), \code{status} (0 for right-censoring,
#' 1 for interval-censoring or left-censoring), and \code{covariates} by column.
#' @param var_list the list of covariates to be fitted into the copula model.
#' @param copula Types of copula model.
#' @param r postive transformation parameter for the semiparametric transformation marginal model.
#' @param m integer, degree of Berstein polynomials for both margins; default is 3
#' @param l the left bound for all \code{Left} and \code{Right} endpoints of observed finite intervals;
#' default is 0.
#' @param u the right bound for all \code{Left} and \code{Right} endpoints of observed finite intervals;
#' has to be a finite value
#' @param method optimization method (see \code{?optim}); default is \code{"BFGS"};
#' also can be \code{"Newton"} (see \code{?nlm}).
#' @param iter number of iterations when \code{method = "Newton"};
#' default is 300.
#' @param stepsize size of optimization step when method is \code{"Newton"};
#' default is 1e-6.
#' @param control a list of control parameters for methods other than \code{"Newton"};
#' see ?optim.
#' @param hes default is \code{TRUE} for hessian calculation;
#' if LRT is desired, can set \code{hes = FALSE} to save time.
#' @importFrom corpcor pseudoinverse
#' @importFrom stats pchisq
#' @importFrom stats optim
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
#' The marginal semiparametric transformation models are built based on Bernstein polynomials, which is formulated below:
#'
#' \deqn{S(t|Z) = \exp[-G\{\Lambda(t) e^{Z^{\top}\beta}\}],} where \eqn{t} is time, \eqn{Z} is covariate,
#' \eqn{\beta} is coefficient and \eqn{\Lambda(t)} is an unspecified function with infinite dimensions.
#' We approximate \eqn{\Lambda(t)} in a sieve space constructed by Bernstein polynomials with degree \eqn{m}. By default, \eqn{m=3}.
#' In the end, all model parameters are estimated by the sieve estimators (Sun et.al., 2019).
#'
#' The \eqn{G(\cdot)} function is the transformation function with a parameter \eqn{r > 0}, which has a form of
#' \eqn{G(x) = \frac{(1+x)^r - 1}{r}}, when \eqn{0 < r \leq 2} and \eqn{G(x) = \frac{\log\{1 + (r-2)x\}}{r - 2}} when \eqn{r > 2}.
#' When \eqn{r = 1}, the marginal model becomes a proportional hazards model;
#' when \eqn{r = 3}, the marginal model becomes a proportional odds model.
#' In practice, \code{m} and \code{r} can be selected based on the AIC value. \cr
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
#' @examples
#' # fit a Copula2-Semiparametric model
#' data(AREDS)
#' copula2_sp <- ic_spTran_copula(data = AREDS, copula = "Copula2",
#'               l = 0, u = 15, m = 3, r = 3,
#'               var_list = c("ENROLLAGE","rs2284665","SevScaleBL"))
#'summary(copula2_sp)



ic_spTran_copula <- function(data, var_list, l=0, u, copula = "Copula2", m = 3, r = 3,
                         method = "BFGS", iter=300, stepsize=1e-6, hes = TRUE,
                         control = list()){

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

  if (l<0 | u <= max(data$Left[is.finite(data$Left)],data$Right[is.finite(data$Right)])) {
    stop('l must be >= 0 and u greater than non-infinite values of Left and Right')
  }

  if (r <= 0) {
    stop('r must be a positive number')
  }

  if (m != round(m)) {
    stop('m must be a positive integer')
  }

  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	{
    stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')
  }


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
  if (method == "Newton") {
    model_step1a <- nlm(estimate_sieve_step1a, rep(0, (p+m+1)), hessian = FALSE,
                        iterlim = iter, steptol = stepsize,
                        p, m = m, x1 = x1, x2 = x2, bl1 = bl1, br1 = br1, bl2 = bl2,
                        br2 = br2, indata1 = indata1, indata2 = indata2, r = r)
    beta <- model_step1a$estimate[1:p]
    phi <- model_step1a$estimate[(p+1):(p+1+m)]
    ep<-cumsum(exp(phi))
  }


  if (method != "Newton") {
    model_step1a <- optim(par = rep(0,(p+m+1)), estimate_sieve_step1a,
                          method = method, hessian = FALSE,
                          p = p, m = m, x1 = x1, x2 = x2, bl1 = bl1, br1 = br1,
                          bl2 = bl2, br2 = br2, indata1 = indata1, indata2 = indata2,
                          r = r, control = control)
    # model_step1a
    beta <- model_step1a$par[1:p]
    phi <- model_step1a$par[(p+1):(p+1+m)]
    ep <- cumsum(exp(phi))
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
    model_step1b <- nlm(ic_copula_log_lik_sieve_eta, eta_ini, hessian = FALSE,
                        beta = beta, ep = ep, x1 = x1, x2 = x2, bl1 = bl1,
                        br1 = br1, bl2 = bl2, br2 = br2, indata1 = indata1,
                        indata2 = indata2, r = r, copula = copula)
    eta_ini <- exp(model_step1b$estimate)
  } else {
    model_step1b <- optim(par = eta_ini, ic_copula_log_lik_sieve_eta,
                          method = method,  control = control, hessian = FALSE,
                         beta = beta, ep = ep, x1 = x1, x2 = x2, bl1 = bl1,
                         br1 = br1, bl2 = bl2, br2 = br2, indata1 = indata1,
                         indata2 = indata2, r = r, copula = copula)
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
    model_step2 <- nlm(ic_copula_log_lik_sieve, c(model_step1a$estimate,eta_ini),
                       hessian = hes, iterlim = iter ,steptol = stepsize,
                       p, m = m, x1 = x1, x2 = x2, bl1 = bl1, br1 = br1, bl2 = bl2,
                       br2 = br2, indata1 = indata1, indata2 = indata2, r = r,
                       copula = copula)

    if (isTRUE(hes)) {
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] <- 0
      dih <- sqrt(dih)
      se <- if (copula != "Copula2") dih[c(1:p,length(dih))] else dih[c(1:p,length(dih)-1,length(dih))]
      beta <- if (copula != "Copula2") model_step2$estimate[c(1:p,length(dih))] else model_step2$estimate[c(1:p,length(dih)-1,length(dih))]
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(model_step2$estimate) - 2 * llk
      stat <- (beta-0)^2/se^2
      pvalue <- pchisq(stat, 1, lower.tail=F)
      summary <- cbind(beta, se, stat, pvalue)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(var_list, tmp_name2)

      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$code
    }

    if (!isTRUE(hes)) {
      inv_info = NULL
      beta <- if (copula != "Copula2")  model_step2$estimate[c(1:p,length(model_step2$estimate))] else model_step2$estimate[c(1:p,length(model_step2$estimate)-1,length(model_step2$estimate))]
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(model_step2$estimate) - 2 * llk
      summary = cbind(beta)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(var_list, tmp_name2)

      colnames(summary) <- c("estimate")
      code <- model_step2$code
    }
    output <- list(code = code, summary = summary, llk = llk, AIC = AIC, copula = copula,
                   m = m, r = r, indata1 = indata1, indata2 = indata2, var_list = var_list,
                   l = l, u = u, bl1 = bl1, br1 = br1, bl2 = bl2, br2 = br2,
                   estimates = model_step2$estimate, x1 = x1, x2 = x2,inv_info  =  inv_info)
  }



  if (method != "Newton") {

    model_step2 <- optim(par = c(model_step1a$par,eta_ini), ic_copula_log_lik_sieve,
                        method = method, hessian = hes, control = control,
                        p = p, m = m, x1 = x1, x2 = x2, bl1 = bl1, br1 = br1,
                        bl2 = bl2, br2 = br2, indata1 = indata1, indata2 = indata2,
                        r = r, copula = copula)

    if (isTRUE(hes)) {
      inv_info <- pseudoinverse(model_step2$hessian)
      dih <- diag(inv_info)
      dih[dih < 0] = 0
      dih <- sqrt(dih)
      se <- if (copula != "Copula2") dih[c(1:p,length(dih))] else dih[c(1:p,length(dih)-1,length(dih))]
      beta <- if (copula != "Copula2") model_step2$par[c(1:p,length(dih))] else model_step2$par[c(1:p,length(dih)-1,length(dih))]
      llk <- -1 * model_step2$value
      AIC <- 2 * length(model_step2$par) - 2 * llk
      stat = (beta-0)^2/se^2
      pvalue = pchisq(stat, 1, lower.tail=F)
      summary = cbind(beta, se, stat, pvalue)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(var_list, tmp_name2)

      colnames(summary) <- c("estimate","SE","stat","pvalue")
      code <- model_step2$convergence
    }

    if (!isTRUE(hes)) {
      inv_info = NULL
      beta <- if (copula != "Copula2")  model_step2$par[c(1:p,length(model_step2$par))] else model_step2$par[c(1:p,length(model_step2$par)-1,length(model_step2$par))]
      llk <- -1 * model_step2$value
      AIC <- 2 * length(model_step2$par) - 2 * llk
      summary <- cbind(beta)

      tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
      rownames(summary) <- c(var_list, tmp_name2)

      colnames(summary) <- c("estimate")
      code <- model_step2$convergence
    }
    output <- list(code = code, summary = summary, llk = llk, AIC = AIC, copula = copula,
                   m = m, r = r, indata1 = indata1, indata2 = indata2, var_list = var_list,
                   l = l, u = u, bl1 = bl1, br1 = br1, bl2 = bl2, br2 = br2,
                   estimates = model_step2$par, x1 = x1, x2 = x2, inv_info = inv_info)
  }

  class(output) <- "CopulaCenR"
  return(output)
}

