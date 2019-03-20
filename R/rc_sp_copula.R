#' Copula regression models with Cox semiparametric margins for bivariate right-censored data
#'
#' @description Fits a copula model with Cox semiparametric margins for bivariate right-censored data.
#'
#' @name rc_spCox_copula
#' @aliases rc_spCox_copula
#' @param data a data frame; must have \code{id} (subject id), \code{ind} (1,2 for two margins),
#' \code{obs_time}, \code{status} (0 for right-censoring, 1 for event).
#' @param var_list the list of covariates to be fitted into the model.
#' @param copula specify the copula family.
#' @param method optimization method (see \code{?optim}); default is \code{"BFGS"};
#' also can be \code{"Newton"} (see \code{?nlm}).
#' @param iter number of iterations when \code{method = "Newton"};
#' default is 500.
#' @param stepsize size of optimization step when \code{method = "Newton"};
#' default is 1e-6.
#' @param control a list of control parameters for methods other than \code{"Newton"};
#' see \code{?optim}.
#' @param B number of bootstraps for estimating standard errors with default 100;
#' @param seed the bootstrap seed; default is 1
#' @importFrom caret createResample
#' @importFrom corpcor pseudoinverse
#' @importFrom survival survreg
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom survival cluster
#' @importFrom stats as.formula
#' @importFrom stats cor
#' @importFrom stats optim
#' @importFrom stats pchisq
#' @importFrom stats quantile
#' @importFrom stats coef
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
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
#' The marginal distribution is a Cox semiparametric proportional hazards model.
#' The copula parameter and coefficient standard errors are estimated from bootstrap. \cr
#'
#' Optimization methods can be all methods (except \code{"Brent"}) from \code{optim},
#' such as \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"}.
#' Users can also use \code{"Newton"} (from \code{nlm}).
#'
#' @return a \code{CopulaCenR} object summarizing the model.
#' Can be used as an input to general \code{S3} methods including
#' \code{summary}, \code{print}, \code{plot}, \code{lines},
#' \code{coef}, \code{logLik}, \code{AIC},
#' \code{BIC}, \code{fitted}, \code{predict}.
#'
#' @examples
#' # fit a Clayton-Cox model
#' data(DRS)
#' clayton_cox <- rc_spCox_copula(data = DRS, var_list = "treat",
#'                             copula = "Clayton", B = 2)
#' summary(clayton_cox)



rc_spCox_copula <- function(data, var_list, copula="Clayton",
                         method = "BFGS", iter = 500, stepsize = 1e-6,
                         control = list(), B = 100,
                         seed = 1) {


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

  ###################################
  ############ Step 1a ##############
  ###################################
  M <- coxph(as.formula(paste0("Surv(obs_time, status)~",
                               paste0(var_list,collapse = "+"),
                               "+cluster(id)")),
             data = x, ties = "breslow")
  beta_ini <- as.numeric(coef(M))

  ###################################
  ############ Step 1b ##############
  ###################################
  if (copula == "AMH") {
    eta_ini<- 0
  }

  else if (copula == "Copula2") {
    eta_ini <- c(0, 0)
  }

  else {
    eta_ini <- 0
  }

  if (method == "Newton") {

    fit0 <- nlm(rc_copula_log_lik_cox_eta, p = eta_ini,
                p2 = c(beta_ini),
                x1 = x1, x2 = x2, indata1 = indata1,indata2 = indata2,
                iterlim = iter, steptol = stepsize, copula = copula)

    eta_ini <- exp(fit0$estimate) # anti-log

  } else {

    fit0 <- optim(par = eta_ini, rc_copula_log_lik_cox_eta,
                  p2 = c(beta_ini),
                  method = method,  control = control, hessian = FALSE,
                  x1 = x1, x2 = x2, indata1 = indata1,indata2 = indata2,
                  copula = copula)
    eta_ini <- exp(fit0$par) # anti-log
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
      model_step2 <- nlm(rc_copula_log_lik_cox, p = c(eta_ini, beta_ini),
                         x1 = x1, x2 = x2,indata1 = indata1,indata2 = indata2,
                         iterlim = iter, steptol = stepsize, hessian = F,
                         copula = copula)
      beta <- model_step2$estimate
      llk <- -1 * model_step2$minimum
      AIC <- 2 * length(beta) - 2 * llk
      code <- model_step2$code
    }


    if (method != "Newton") {
      model_step2 <- optim(c(eta_ini, beta_ini), rc_copula_log_lik_cox,
                           x1 = x1, x2 = x2, indata1 = indata1, indata2 = indata2,
                           hessian = F, method = method, control = control,
                           copula = copula)
      beta <- model_step2$par
      llk <- -1 * model_step2$value
      AIC <- 2 * length(beta) - 2 * llk
      code <- model_step2$convergence
    }


    # partitiioning samples by adjusting censoring
    set.seed(seed)
    status.id <- aggregate(data$status, by = list(data$id), FUN = sum)$x
    groups <- createResample(status.id, times = B, list = F)
    tmp <- c()
    pb = txtProgressBar(title = "progress bar",
                        min = 0, max = B, initial = 0,
                        style = 3)  # progress bar

    for (K in 1:B) {

      tryCatch({


      # new data
      b = groups[, K]
      x1.b = x1[b, ]
      x2.b = x2[b, ]
      indata1.b = indata1[b, ]
      indata2.b = indata2[b, ]

      if (method == "Newton") {
        model_step2.b <- nlm(rc_copula_log_lik_cox, p = c(eta_ini, beta_ini),
                             x1 = x1.b, x2 = x2.b, indata1 = indata1.b,
                             indata2 = indata2.b,
                             iterlim = iter, steptol = stepsize, hessian = F,
                             copula = copula)
        beta.b <- model_step2.b$estimate
      }


      if (method != "Newton") {
        model_step2.b <- optim(c(eta_ini, beta_ini), rc_copula_log_lik_cox,
                               x1 = x1.b, x2 = x2.b, indata1 = indata1.b,
                               indata2 = indata2.b,
                               hessian = F, method = method, control = control,
                               copula = copula)
        beta.b <- model_step2.b$par
      }

    }, error=function(e){cat("Error",conditionMessage(e), "\n")})

      tmp <- cbind(tmp, beta.b)
      setTxtProgressBar(pb = pb, value = K, label = paste0(round(K/B*100, 0), "% done"))

    }
    close(pb)

    # calculate s.e
    se <- as.numeric(apply(tmp, 1, function(x) {sd(x, na.rm = T)}))
    stat <- (beta - 0)^2/se^2
    pvalue <- pchisq(stat, 1, lower.tail = F)

    # summary
    summary <- cbind(beta, se, stat, pvalue)
    tmp_name <- if (copula != "Copula2") c("eta") else c("alpha","kappa")
    rownames(summary) = c(tmp_name, var_list)
    colnames(summary) <- c("estimate","SE","stat","pvalue")
    summary <- summary[c(var_list, tmp_name), ]

    # output
    if (method == "Newton") {
      output <- list(code = code, llk = llk, AIC = AIC, summary = summary,
                     copula = copula, indata1 = indata1,
                     indata2 = indata2, var_list = var_list,
                     estimates = model_step2$estimate, x1 = x1, x2 = x2,
                     cox = TRUE)
    } else {
      output <- list(code = code, llk = llk, AIC = AIC, summary = summary,
                     copula = copula, indata1 = indata1,
                     indata2 = indata2, var_list = var_list,
                     estimates = model_step2$par, x1 = x1, x2 = x2,
                     cox = TRUE)
    }


  class(output) <- "CopulaCenR"
  return(output)

}

