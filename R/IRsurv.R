#' An information ratio-based goodness-of-fit test for copula models on censored data
#'
#' Fits an Information ratio (IR)-based goodness-of-fit test for copula models under various censoring types.
#'
#' @name IRsurv
#' @aliases IRsurv
#' @param data the input data; see examples for details.
#' @param censoring types of censoring, including "rc", "ic", "rec_bivariate", "rec_multivariate".
#' @param copula specify the copula family to be tested; default is "clayton"; others include "copula2", "gumbel", "frank", "gaussian", and "d-vine". "d-vine" is only for censoring = "rec_multivariate".
#' @param fams specify the unconditional copulas by following the style of the VineCopula package when copula = "d-vine". Only d = 4 is supported at this stage. The conditional copulas are set as "frank" by default.
#' @param R number of Bootstraps; default is 200.
#' @param parallel indicator of parallel computing; can be "no" or "multicore"; default is "no".
#' @param ncpus number of cpus to be assigned for parallel computing; default is 1.
#' @import boot
#' @importFrom copula archmCopula
#' @importFrom copula rCopula
#' @importFrom copula normalCopula
#' @importFrom copula plackettCopula
#' @importFrom copula dCopula
#' @importFrom copula cCopula
#' @importFrom copula pCopula
#' @importFrom copula iTau
#' @importFrom copBasic PLcop
#' @importFrom copBasic derCOP
#' @importFrom copBasic derCOP2
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom corpcor pseudoinverse
#' @importFrom pracma grad
#' @importFrom pracma hessian
#' @importFrom stats as.formula
#' @importFrom stats as.formula
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom stats na.omit
#' @importFrom stats reshape
#' @importFrom stats rexp
#' @importFrom stats runif
#' @import VineCopula
#' @importFrom flexsurv flexsurvreg
#' @import icenReg
#' @return the p value of the IR test
#' @export
#' @source
#' Tao Sun, Yu Cheng, Ying Ding (2022). An information Ratio-based
#' Goodness-of-Fit Test for Copula Models on Censored Data. Biometrics (Accepted).
#' @examples
#'\dontrun{
#' # Goodness of fit under right censoring
#' data("data_sim_RC")
#' test_rc <- IRsurv(data = data_sim_RC, censoring = "rc", copula = "clayton", R = 200)
#' test_rc
#' # Goodness of fit under interval censoring
#' data("data_sim_ic")
#' test_ic <- IRsurv(data = data_sim_ic, censoring = "ic", copula = "clayton", R = 200)
#' test_ic
#' # Goodness of fit under bivariate recurrent events
#' data("data_sim_rec")
#' test_rec_bi <- IRsurv(data = data_sim_rec, censoring = "rec_bivariate",
#'                       copula = "clayton", R = 200)
#' test_rec_bi
#' # Goodness of fit of D-vine copula under multivariate recurrent events
#' data("data_sim_multi_rec")
#' test_rec_mv <- IRsurv(data = data_sim_multi_rec, censoring = "rec_multivariate",
#'                       fams = c(3,3,3), R = 200)
#' test_rec_mv
#'}


IRsurv <- function(data, censoring = "rc", copula = "clayton", fams = c(3,3,3), R = 200, parallel = "no", ncpus = 1){



  if (censoring == "rc") {
    data_sim_RC = data
    fit1 <- Marginal_KM(data = data_sim_RC)
    fit2 <- MLE_copula_np_PMLE_eta(data = data_sim_RC, copula = copula)
    fit <- list(param=c(fit2), S1=fit1$S1, S2=fit1$S2, copula = copula)
    bs.res <- boot(R = R , data = data_sim_RC, statistic = copula_param_bootstrap_IR, copula = copula,
                   sim = "parametric", ran.gen = gen_copula_RC, mle = fit,
                   parallel = "multicore", ncpus = ncpus)
  }


  if (censoring == "ic") {
    data_sim_ic = data
    fit1 <- Marginal_Turnbull(data_sim_ic)
    fit2 <- MLE_ic_copula_Turnbull_eta(data_sim_ic, copula = copula)
    fit <- list(param=c(fit2), u1_left=fit1$u1_left, u1_right=fit1$u1_right, u2_left=fit1$u2_left, u2_right=fit1$u2_right, copula = copula)
    bs.res <- boot(R = R , data = data_sim_ic, statistic = ic_copula_param_bootstrap_IR, copula = copula,
                   sim = "parametric", ran.gen = gen_copula_IC, mle = fit,
                   parallel = "multicore", ncpus = ncpus)
  }

  if (censoring == "rec_bivariate") {
    data_sim_rec = data
    fit1 <- modified_KM(data = data_sim_rec)
    fit2 <- MLE_copula_np_PMLE_eta_rec(data = data_sim_rec, copula = copula)
    fit <- list(param=c(fit2), S1=fit1$copData[,1], S2=fit1$copData[,2], status1 = fit1$status[,1], status2 = fit1$status[,2], copula = copula)
    bs.res <- boot(R = R , data = data_sim_rec, statistic = copula_param_bootstrap_IR_rec, copula = copula,
                   sim = "parametric", ran.gen = gen_copula_RC_rec, mle = fit,
                   parallel = "multicore", ncpus = 1) #
  }


  if (censoring == "rec_multivariate") {
    data_sim_multi_rec = data
    fit1 <- modified_KM_vine(data = data_sim_multi_rec)
    fit2 <- MLE_copula_np_PMLE_eta_rec_vine(data = data_sim_multi_rec, fams = fams)
    fit <- list(param=c(fit2), S1=fit1$copData[,1], S2=fit1$copData[,2],
                S3=fit1$copData[,3], S4=fit1$copData[,4],
                status1 = fit1$status[,1], status2 = fit1$status[,2],
                status3 = fit1$status[,3], status4 = fit1$status[,4],
                fams = fams)
    bs.res <- boot(R = R , data = data_sim_multi_rec, statistic = copula_param_bootstrap_IR_rec_vine, fams = fams,
                   sim = "parametric", ran.gen = gen_copula_RC_rec_vine, mle = fit,
                   parallel = "multicore", ncpus = ncpus)
  }

  summary <- data.frame(empir.pvalue = mean(abs(bs.res$t) > abs(bs.res$t0)))
  return(summary)

}

