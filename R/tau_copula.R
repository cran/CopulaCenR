#' Calculate Kendall's tau
#'
#' This function (tau_copula) is to obtain Kendall's tau from copula parameter(s)
#'
#' @name tau_copula
#' @aliases tau_copula
#' @param eta copula parameter(s); if Coupla2, input as a vector of (alpha, kappa)
#' @param copula specify the type of copula model
#' @importFrom copula archmCopula
#' @importFrom copula tau
#' @return Kendall's tau
#' @export


tau_copula <- function(eta, copula){

  if (tolower(copula) == "copula2"){
    alpha <- eta[1]
    kappa <- eta[2]
    output <- 1-2*alpha*kappa/(1+2*kappa)
  }

  if (tolower(copula) != "copula2") {
    output <- tau(archmCopula(tolower(copula), param = eta, dim = 2))
  }

  return(output)
}
