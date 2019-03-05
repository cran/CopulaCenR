#' Generalized score test for covariate effect(s)
#'
#' Generalized score test on covariate effect(s) under a fitted copula model.
#'
#' @name score_copula
#' @aliases score_copula
#' @param object The output object from the main functions
#' (\code{rc_par_copula}, \code{ic_sp_copula}, \code{ic_par_copula})
#' under the null hypothesis
#' @param var_score the list of covariates to be tested by the score test
#' @importFrom corpcor pseudoinverse
#' @importFrom stats nlm
#' @importFrom pracma grad
#' @importFrom pracma hessian
#' @return the score statistics, p value
#' @export
#'
#' @examples
#' # Score test for "rs2284665" in AREDS data
#' # fit a Copula2-semiparametric model under NULL
#' data(AREDS)
#' copula2_sp_null <- ic_sp_copula(data = AREDS, copula = "Copula2",
#'                    l = 0, u = 15, m = 3, r = 3,
#'                    var_list = c("ENROLLAGE","SevScaleBL"))
#' score_copula(object = copula2_sp_null, var_score = "rs2284665")


score_copula <- function(object, var_score){


  # IC, transformation model
  if (is.numeric(object$m)){

    output <- ic_sp_copula_score(object, var_score)

  }

  # IC, parametric margins
  else if (!is.numeric(object$m) & !("obs_time" %in% colnames(object$indata1)) ) {

    output <- ic_par_copula_score(object, var_score)

  }

  # RC, parametric margins
  else if (!is.numeric(object$m) & ("obs_time" %in% colnames(object$indata1)) ) {

    output <- rc_par_copula_score(object, var_score)

  }

  return(output)
}
