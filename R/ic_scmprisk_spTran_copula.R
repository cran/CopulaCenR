#' Copula regression models with semi-parametric transformation margins for semi-competing risk data under interval-censoring and left-truncation
#'
#' @description Fits a copula model with semi-parametric transformation margins for semi-competing risk data under interval-censoring and left-truncation.
#'
#' @name ic_scmprisk_spTran_copula
#' @aliases ic_scmprisk_spTran_copula
#' @param data a data frame; must have \code{id} (subject id), \code{Left} (0 if left-censoring for non-terminal event),
#' \code{Right} (Inf if right-censoring for non-terminal event), \code{status} (0 for right-censoring,
#' 1 for interval-censoring or left-censoring of non-terminal event), \code{timeD} (observed terminal event),
#' \code{statusD} (for terminal event), \code{A} (left truncation time, 0 if none), and \code{covariates} by column.
#' @param var_list the list of covariates to be fitted into the copula model.
#' @param copula Types of copula model, only Copula2 is supported at this stage.
#' @param r1 for non-terminal event, postive transformation parameter for the semiparametric transformation marginal model.
#' @param m1 for non-terminal event, integer, degree of Berstein polynomials for both margins; default is 3
#' @param l1 for non-terminal event, the left bound for all \code{Left} and \code{Right} endpoints of observed finite intervals;
#' default is 0.
#' @param u1 for non-terminal event, the right bound for all \code{Left} and \code{Right} endpoints of observed finite intervals;
#' has to be a finite value
#' @param r2 for terminal event, postive transformation parameter for the semiparametric transformation marginal model.
#' @param m2 for terminal event, integer, degree of Berstein polynomials for both margins; default is 3
#' @param l2 for terminal event, the left bound for all \code{Left} and \code{Right} endpoints of observed finite intervals;
#' default is 0.
#' @param u2 for terminal event, the right bound for all \code{Left} and \code{Right} endpoints of observed finite intervals;
#' has to be a finite value
#' @param method optimization method (see ?optim); default is "BFGS";
#' also can be "Newton" (see ?nlm).
#' @param iter number of iterations when method is \code{"Newton"};
#' default is 300.
#' @param stepsize size of optimization step when method is \code{"Newton"};
#' default is 1e-5.
#' @param control a list of control parameters for methods other than \code{"Newton"};
#' see \code{?optim}.
#' @param eta_ini a vector of initial values for copula parameters, default is NULL
#' @importFrom corpcor pseudoinverse
#' @importFrom stats pchisq
#' @importFrom stats optim
#' @export
#'
#' @source
#' Tao Sun, Yunlong Li, Zhengyan Xiao, Ying Ding, Xiaojun Wang (2022).
#' Semiparametric copula method for semi-competing risks data
#' subject to interval censoring and left truncation:
#' Application to disability in elderly. Statistical Methods in Medical Research (Accepted).
#'
#' @details The input data must be a data frame. with columns \code{id} (sample id),
#' \code{Left} (0 if left-censoring), \code{Right} (Inf if right-censoring),
#' \code{status} (0 for right-censoring, 1 for interval-censoring or left-censoring),
#' \code{timeD} (for terminal event), \code{statusD},\code{A} (0 if no left truncation),
#' and \code{covariates}. The function does not allow \code{Left} == \code{Right}. \cr
#'
#'
#' The supported copula model in this version is \code{"Copula2"}.
#' The \code{"Copula2"} model is a two-parameter copula model that incorporates \code{Clayton}
#' and \code{Gumbel} as special cases.
#' The parametric generator functions of copula functions are list below:
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
#' In the end, all model parameters are estimated by the sieve estimators (Sun and Ding, In Press).
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
#' \code{summary}, \code{print}, \code{coef},
#' \code{logLik}, \code{AIC}, \code{BIC}.
#' @examples
#' # fit a Copula2-Semiparametric model
#' data("data_scmprisk")
#' copula2_sp <- ic_scmprisk_spTran_copula(data = data_scmprisk,
#'               var_list = c("x1"), copula = "Copula2",
#'               l1=0, u1 = 21, m1 = 3, r1 = 1,
#'               l2=0, u2 = 21, m2 = 3, r2 = 1,
#'               )
#'summary(copula2_sp)



ic_scmprisk_spTran_copula <- function(data, var_list, copula = "Copula2",
                                      l1=0, u1, m1 = 3, r1 = 1,
                                      l2=0, u2, m2 = 3, r2 = 1,
                                      method = "BFGS", iter=1000, stepsize=1e-5,
                                      control = list(), eta_ini = NULL){


  # first screen the inputs: copula, m.dist, method #
  if (!is.data.frame(data)) {
    stop('data must be a data frame')
  }

  if ((!"id" %in% colnames(data)) |
      (!"Left" %in% colnames(data)) |
      (!"Right" %in% colnames(data)) |
      (!"status" %in% colnames(data)) |
      (!"timeD" %in% colnames(data)) |
      (!"statusD" %in% colnames(data)) |
      (!"A" %in% colnames(data))) {
    stop('data must have id, Left, Right, status, timeD, statusD and A')
  }

  # if (!copula %in% c("Clayton","Gumbel","Copula2","Frank"))	{
  #   stop('copula must be one of "Clayton","Gumbel","Copula2","Frank"')
  # }

  if (!copula %in% c("Copula2"))	{
    stop('copula must be "Copula2"')
  }

  if (l1<0 | u1 <= max(data$Left[is.finite(data$Left)],data$Right[is.finite(data$Right)])) {
    stop('l1 must be >= 0 and u1 greater than non-infinite values of Left and Right')
  }

  if (l2<0 | u2 <= max(data$timeD[is.finite(data$timeD)])) {
    stop('l2 must be >= 0 and u2 greater than non-infinite values of timeD')
  }

  if (r1 <= 0) {
    stop('r1 must be a positive number')
  }

  if (r2 <= 0) {
    stop('r2 must be a positive number')
  }

  if (m1 != round(m1)) {
    stop('m1 must be a positive integer')
  }

  if (m2 != round(m2)) {
    stop('m2 must be a positive integer')
  }

  if (!method %in% c("Newton","Nelder-Mead","BFGS","CG","SANN"))	{
    stop('m.dist must be one of "Newton","Nelder-Mead","BFGS","CG","SANN"')
  }

  # data pre-processing #
  data$A1 = data$A2 = data$A
  data_processed <- data_process_scmprisk_ic_sp_LT_A1A2(data, var_list, l1, u1, m1, l2, u2, m2)
  indata1 <- data_processed$indata1
  indata2 <- data_processed$indata2
  t1_left <- data_processed$t1_left
  t1_right <- data_processed$t1_right
  t2 <- data_processed$t2
  A1 <- data_processed$A1
  A2 <- data_processed$A2
  n <- data_processed$n
  p <- data_processed$p
  x1 <- data_processed$x1
  x2 <- data_processed$x2
  var_list <- data_processed$var_list
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]

  # BP
  bl1 <- data_processed$bl1
  br1 <- data_processed$br1
  b2 <- data_processed$b2
  b1_A <- data_processed$b1_A
  b2_A <- data_processed$b2_A

  # BP derivatives
  bl1_d <- data_processed$bl1_d
  br1_d <- data_processed$br1_d
  b2_d <- data_processed$b2_d

  ###################################
  ############ Step 1a ##############
  ###################################

  ###### obtain consistent estimator for timeD under RC and LT (in fact, LT is not considered due to computational issues)#######
  x2_1a <- data.frame(timeD = c(indata2$timeD),
                      statusD = c(indata2$statusD),
                      A = A2,
                      x2)

  if (method == "Newton") {

    # omit LT
    model_step1a <- nlm(estimate_sieve_step1a_scmprisk_rc_LT, rep(0.1, (p2+m2+1)), hessian = F,
                        iterlim = iter, steptol = stepsize,
                        p2=p2, m2 = m2, x2 = x2, b2 = b2,  b2_d = b2_d, b2_A = b2_A,
                        indata2 = x2_1a, r2 = r2)

    # adjust LT
    model_step1a <- nlm(estimate_sieve_step1a_scmprisk_rc_LT_2, model_step1a$estimate, hessian = F,
                        iterlim = iter, steptol = stepsize,
                        p2=p2, m2 = m2, x2 = x2, b2 = b2,  b2_d = b2_d, b2_A = b2_A,
                        indata2 = x2_1a, r2 = r2)

    beta_ini_2 <- model_step1a$estimate[1:p2]
    phi_ini_2 <- model_step1a$estimate[(p2+1):(p2+1+m2)]
  }


  if (method != "Newton") {

    # omitting LT
    model_step1a <- optim(par = rep(0, (p2+m2+1)), estimate_sieve_step1a_scmprisk_rc_LT,
                          method = method, hessian = F,
                          p2 = p2, m2 = m2, x2 = x2,
                          b2 = b2, b2_d = b2_d, indata2 = x2_1a, b2_A = b2_A,
                          r2 = r2, control = control)

    # adjusting LT
    model_step1a <- optim(par = model_step1a$par, estimate_sieve_step1a_scmprisk_rc_LT_2,
                          method = method, hessian = F,
                          p2 = p2, m2 = m2, x2 = x2,
                          b2 = b2, b2_d = b2_d, indata2 = x2_1a, b2_A = b2_A,
                          r2 = r2, control = control)

    # model_step1a
    beta_ini_2 <- model_step1a$par[1:p2]
    phi_ini_2 <- model_step1a$par[(p2+1):(p2+1+m2)]
  }


  ###### obtain initial (inconsistent) estimator for event1 under IC and LT #######
  x1_1a <- data.frame(Left = c(indata1$Left),
                      Right = c(indata1$Right),
                      status = indata1$status,
                      A = indata1$A1,
                      (x1))

  if (method == "Newton") {
    model_step1a <- nlm(estimate_sieve_step1a_scmprisk_ic_LT_1, rep(0, (p1+m1+1)), hessian = FALSE,
                        iterlim = iter, steptol = stepsize,
                        p1, m1 = m1, x1 = x1, bl1 = bl1, br1 = br1, b1_A = b1_A,
                        indata1 = x1_1a, r1 = r1)


    model_step1a <- nlm(estimate_sieve_step1a_scmprisk_ic_LT_2, model_step1a$estimate, hessian = FALSE,
                        iterlim = iter, steptol = stepsize,
                        p1, m1 = m1, x1 = x1, bl1 = bl1, br1 = br1, b1_A = b1_A,
                        indata1 = x1_1a, r1 = r1)

    beta_ini_1 <- model_step1a$estimate[1:p1]
    phi_ini_1 <- model_step1a$estimate[(p1+1):(p1+1+m1)]
  }


  if (method != "Newton") {
    model_step1a <- optim(par = rep(0, (p1+m1+1)), estimate_sieve_step1a_scmprisk_ic_LT_1,
                          method = method, hessian = FALSE,
                          p1 = p1, m1 = m1, x1 = x1, bl1 = bl1, br1 = br1, b1_A = b1_A,
                          indata1 = x1_1a,
                          r1 = r1, control = control)

    model_step1a <- optim(par = model_step1a$par, estimate_sieve_step1a_scmprisk_ic_LT_2,
                          method = method, hessian = FALSE,
                          p1 = p1, m1 = m1, x1 = x1, bl1 = bl1, br1 = br1, b1_A = b1_A,
                          indata1 = x1_1a,
                          r1 = r1, control = control)

    # model_step1a
    beta_ini_1 <- model_step1a$par[1:p1]
    phi_ini_1 <- model_step1a$par[(p1+1):(p1+1+m1)]
  }


  ###################################
  ############ Step 1b ##############
  ###################################

  if (is.null(eta_ini)) {
    if (copula == "AMH") {
      eta_ini <- 1
    }

    else if (copula == "Copula2") {
      eta_ini <- c(log(0.5/0.5), log(1))
    }

    else {
      eta_ini <- 1
    }
  }


  if (method == "Newton") {

    fit0 <- nlm(ic_scmprisk_copula_log_lik_sieve_pseudo_LT_1_copula2, p = c(eta_ini, phi_ini_1, beta_ini_1),
                fitted = c(phi_ini_2,beta_ini_2),
                x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right, t2=t2, indata1 = indata1,indata2 = indata2,
                bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                b1_A = b1_A, b2_A = b2_A,
                iterlim = iter, steptol = stepsize, copula = copula)

    fit0 <- nlm(ic_scmprisk_copula_log_lik_sieve_pseudo_LT_copula2, p = fit0$estimate,
                fitted = c(phi_ini_2,beta_ini_2),
                x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right, t2=t2, indata1 = indata1,indata2 = indata2,
                bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                b1_A = b1_A, b2_A = b2_A,
                iterlim = iter, steptol = stepsize, copula = copula)

    if (copula == "Copula2") {
      p_ini <- c((fit0$estimate[1]), (fit0$estimate[2]), fit0$estimate[3:length(fit0$estimate)]) # anti-log
    } else {
      p_ini <- c((fit0$estimate[1]), fit0$estimate[2:length(fit0$estimate)]) # anti-log
    }


  } else {

    fit0 <- optim(par = c(eta_ini, phi_ini_1, beta_ini_1),
                  ic_scmprisk_copula_log_lik_sieve_pseudo_LT_1_copula2,
                  fitted = c(phi_ini_2, beta_ini_2),
                  method = method,  control = control, hessian = FALSE,
                  x1 = x1, x2 = x2,indata1 = indata1,indata2 = indata2,
                  t1_left = t1_left, t1_right = t1_right, t2=t2,
                  bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                  b1_A = b1_A, b2_A = b2_A,
                  copula = copula)

    fit0 <- optim(par = fit0$par,
                  ic_scmprisk_copula_log_lik_sieve_pseudo_LT_copula2,
                  fitted = c(phi_ini_2, beta_ini_2),
                  method = method,  control = control, hessian = FALSE,
                  x1 = x1, x2 = x2,indata1 = indata1,indata2 = indata2,
                  t1_left = t1_left, t1_right = t1_right, t2=t2,
                  bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                  b1_A = b1_A, b2_A = b2_A,
                  copula = copula)

    if (copula == "Copula2") {
      p_ini <- c((fit0$par[1]), (fit0$par[2]), fit0$par[3:length(fit0$par)]) # anti-log
    } else {
      p_ini <- c((fit0$par[1]), fit0$par[2:length(fit0$par)]) # anti-log
    }

  }


  ###################################
  ############ Step 2 ###############
  ###################################
  if (method == "Newton") {

    model_step2 <- nlm(ic_scmprisk_copula_log_lik_sieve_LT_1_copula2, p = c(p_ini, phi_ini_2, beta_ini_2), # eta, margin1,margin2
                       x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right, t2=t2, indata1 = indata1,indata2 = indata2,
                       bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                       b1_A = b1_A, b2_A = b2_A,
                       iterlim = iter, steptol = stepsize, hessian = F,
                       copula = copula)

    # alpha = expit(alpha0), kappa = exp(alpha0)
    model_step2 <- nlm(ic_scmprisk_copula_log_lik_sieve_LT_copula2, p = model_step2$estimate, # eta, margin1,margin2
                       x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right, t2=t2, indata1 = indata1,indata2 = indata2,
                       bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                       b1_A = b1_A, b2_A = b2_A,
                       iterlim = iter, steptol = stepsize, hessian = T,
                       copula = copula)
    inv_info <- pseudoinverse(model_step2$hessian)
    dih <- diag(inv_info)
    dih[dih < 0] <- 0
    se <- sqrt(dih)
    beta <- if (copula != "Copula2") c((model_step2$estimate[1]), model_step2$estimate[c((1+1+m1+1):(2+m1+p1), (2+m1+p1+1+m2+1):(2+m1+p1+1+m2+p2))]) else c(exp((model_step2$estimate[1]))/(1+exp(model_step2$estimate[1])), exp(model_step2$estimate[2]), model_step2$estimate[c((2+1+m1+1):(2+1+m1+p1), (2+1+m1+p1+1+m2+1):(2+1+m1+p1+1+m2+p2))])
    se <- if (copula != "Copula2") c(se[1], se[c((1+1+m1+1):(2+m1+p1), (2+m1+p1+1+m2+1):(2+m1+p1+1+m2+p2))]) else c(se[1]*exp(model_step2$estimate[1])/(1+exp(model_step2$estimate[1]))^2, se[2]*exp(model_step2$estimate[2]), se[c((2+1+m1+1):(2+1+m1+p1), (2+1+m1+p1+1+m2+1):(2+1+m1+p1+1+m2+p2))])
    llk <- -1 * model_step2$minimum
    AIC <- 2 * length(model_step2$estimate) - 2 * llk
    stat <- (beta - 0)^2/se^2
    pvalue <- pchisq(stat,1,lower.tail=F)
    summary <- cbind(beta, se, stat, pvalue)

    tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
    rownames(summary) = c(tmp_name2, paste0(c(var_list),"_1"), paste0(c(var_list),"_2"))
    colnames(summary) <- c("estimate","SE","stat","pvalue")
    code <- model_step2$code
    output <- list(code = code, summary = summary, llk = llk, AIC = AIC,
                   copula = copula, indata1 = indata1,
                   indata2 = indata2, var_list = var_list,
                   estimates = model_step2$estimate, x1 = x1, x2 = x2,
                   inv_info = inv_info,
                   bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, b1_A = b1_A, b2_A = b2_A,
                   m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2, data = data)

  }


  if (method != "Newton") {

    model_step2 <- optim(c(p_ini, phi_ini_2, beta_ini_2), ic_scmprisk_copula_log_lik_sieve_LT_1_copula2,
                         x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right, t2=t2, indata1 = indata1,indata2 = indata2,
                         bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                         b1_A = b1_A, b2_A = b2_A,
                         hessian = F, method = method,  control = control,
                         copula = copula)

    model_step2 <- optim(model_step2$par, ic_scmprisk_copula_log_lik_sieve_LT_copula2,
                         x1 = x1, x2 = x2, t1_left = t1_left, t1_right = t1_right, t2=t2, indata1 = indata1,indata2 = indata2,
                         bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2,
                         b1_A = b1_A, b2_A = b2_A,
                         hessian = T, method = method,  control = control,
                         copula = copula)

    inv_info <- pseudoinverse(model_step2$hessian)
    dih <- diag(inv_info)
    dih[dih < 0] <- 0
    se <- sqrt(dih)
    beta <- if (copula != "Copula2") c((model_step2$par[1]), model_step2$par[c((1+1+m1+1):(2+m1+p1), (2+m1+p1+1+m2+1):(2+m1+p1+1+m2+p2))]) else c(exp((model_step2$par[1]))/(1+exp(model_step2$par[1])), exp(model_step2$par[2]), model_step2$par[c((2+1+m1+1):(2+1+m1+p1), (2+1+m1+p1+1+m2+1):(2+1+m1+p1+1+m2+p2))])
    se <- if (copula != "Copula2") c(se[1], se[c((1+1+m1+1):(2+m1+p1), (2+m1+p1+1+m2+1):(2+m1+p1+1+m2+p2))]) else c(se[1]*exp(model_step2$par[1])/(1+exp(model_step2$par[1]))^2, se[2]*exp(model_step2$par[2]), se[c((2+1+m1+1):(2+1+m1+p1), (2+1+m1+p1+1+m2+1):(2+1+m1+p1+1+m2+p2))])
    llk <- -1 * model_step2$value
    AIC <- 2 * length(model_step2$par) - 2 * llk
    stat <- (beta - 0)^2/se^2
    pvalue <- pchisq(stat, 1, lower.tail=F)
    summary <- cbind(beta, se, stat, pvalue)

    tmp_name2 <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
    rownames(summary) = c(tmp_name2, paste0(c(var_list),"_1"), paste0(c(var_list),"_2"))
    colnames(summary) <- c("estimate","SE","stat","pvalue")

    code <- model_step2$convergence
    output <- list(code = code, summary = summary, llk = llk, AIC = AIC,
                   copula = copula, indata1 = indata1,
                   indata2 = indata2, var_list = var_list,
                   estimates = model_step2$par, x1 = x1, x2 = x2,
                   inv_info = inv_info,
                   bl1 = bl1, br1 = br1, b2 = b2, b2_d = b2_d, b1_A = b1_A, b2_A = b2_A,
                   m1 = m1, m2 = m2, r1 = r1, r2 = r2, p1 = p1, p2 = p2, data = data)
  }

  class(output) <- "CopulaCenR"
  return(output)
}


