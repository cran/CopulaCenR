#' Likelihood-ratio test for covariate effect(s) in copula models
#'
#' This function (lrt_copula) is used to perform the likelihood ratio test (LRT) between two nested copula models
#'
#' @name lrt_copula
#' @aliases lrt_copula
#' @param model1 The output of the larger model
#' @param model2 The output of the smaller model
#' @importFrom stats pchisq
#' @return the LRT statistics, p value
#' @export
#'
#' @examples
#' #' # Likelihood-ratio test for "rs2284665" in AREDS data
#' data(AREDS)
#' # Fit null model without "rs2284665"
#' copula2_sp_null <- ic_spTran_copula(data = AREDS, copula = "Copula2",
#'                    l = 0, u = 15, m = 3, r = 3,
#'                    var_list = c("SevScaleBL"))
#' # Fit full model
#' copula2_sp <- ic_spTran_copula(data = AREDS, copula = "Copula2",
#'               l = 0, u = 15, m = 3, r = 3,
#'               var_list = c("rs2284665","SevScaleBL"))
#' lrt_copula(model1 = copula2_sp, model2 = copula2_sp_null)


lrt_copula <- function(model1, model2){

  llk.1 <- model1$llk
  p.1 <- length(model1$var_list)

  llk.2 <- model2$llk
  p.2 <- length(model2$var_list)

  stat = llk.1 - llk.2
  df = p.1 - p.2
  pvalue = pchisq(stat,df,lower.tail=F)
  output = c(stat, pvalue)
  names(output) = c("stat", "pvalue")

  return(output)
}







