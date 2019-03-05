#' Printing outputs of a CopulaCenR object
#' @name print.CopulaCenR
#' @aliases print.CopulaCenR
#' @param x a CopulaCenR object
#' @param ... further arguments
#' @importFrom stats printCoefmat
#' @export
print.CopulaCenR <- function(x,...) {
  if (!is.null(x$m.dist)){
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  ",x$m.dist,"\n")
    cat("\n")
    if (dim(x$summary)[2] > 1){
      printCoefmat(x$summary, P.values = T, has.Pvalue = T)
      cat("(The Wald tests are testing whether each coefficient is 0)","\n")
    } else {
      printCoefmat(x$summary, P.values = F, has.Pvalue = F)
    }
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    if (x$code == 0) {cat("Convergence is completed successfully","\n")}
  } else {
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  semiparametric","\n")
    cat("\n")
    if (dim(x$summary)[2] > 1){
      printCoefmat(x$summary, P.values = T, has.Pvalue = T)
      cat("(The Wald tests are testing whether each coefficient is 0)","\n")
    } else {
      printCoefmat(x$summary, P.values = F, has.Pvalue = F)
    }
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    if (x$code == 0) {cat("Convergence is completed successfully","\n")}
  }
}


#' Summarizing outputs of a CopulaCenR object
#' @name summary.CopulaCenR
#' @aliases summary.CopulaCenR
#' @param object a CopulaCenR object
#' @param ... further arguments
#' @export
summary.CopulaCenR <- function(object,...) {
  res <- list(copula=object$copula, m.dist=object$m.dist, summary=object$summary,
              llk=object$llk, AIC=object$AIC, code=object$code)
  class(res) <- "summary.CopulaCenR"
  res
}

#' Print the summary of a CopulaCenR object
#' @name print.summary.CopulaCenR
#' @aliases print.summary.CopulaCenR
#' @param x a summary.CopulaCenR object
#' @param ... further arguments
#' @importFrom stats printCoefmat
#' @export
print.summary.CopulaCenR <- function(x,...) {
  if (!is.null(x$m.dist)){
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  ",x$m.dist,"\n")
    cat("\n")
    if (dim(x$summary)[2] > 1){
      printCoefmat(x$summary, P.values = T, has.Pvalue = T)
      cat("(The Wald tests are testing whether each coefficient is 0)","\n")
    } else {
      printCoefmat(x$summary, P.values = F, has.Pvalue = F)
    }
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    if (x$code == 0) {cat("Convergence is completed successfully","\n")}
  } else {
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  semiparametric","\n")
    cat("\n")
    if (dim(x$summary)[2] > 1){
      printCoefmat(x$summary, P.values = T, has.Pvalue = T)
      cat("(The Wald tests are testing whether each coefficient is 0)","\n")
    } else {
      printCoefmat(x$summary, P.values = F, has.Pvalue = F)
    }
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    if (x$code == 0) {cat("Convergence is completed successfully","\n")}
  }
}

#' the coefficient estimates of a CopulaCenR object
#' @name coef.CopulaCenR
#' @aliases coef.CopulaCenR
#' @param object a CopulaCenR object
#' @param ... further arguments
#' @export
coef.CopulaCenR <- function(object,...) {

  res = object$summary[,1]
  res

}


#' the log-likelihood of a CopulaCenR object
#' @name logLik.CopulaCenR
#' @aliases logLik.CopulaCenR
#' @param object a CopulaCenR object
#' @param ... further arguments
#' @importFrom stats logLik
#' @export
logLik.CopulaCenR <- function(object,...) {

  res <- object$llk
  res

}


#' the AIC of a CopulaCenR object
#' @name AIC.CopulaCenR
#' @aliases AIC.CopulaCenR
#' @param object a CopulaCenR object
#' @param ... further arguments
#' @param k numeric, with k = 2 for AIC
#' @importFrom stats AIC
#' @export
AIC.CopulaCenR <- function(object, ..., k = 2) {

  res <- object$AIC
  return(res)

}


#' the BIC of a CopulaCenR object
#' @name BIC.CopulaCenR
#' @aliases BIC.CopulaCenR
#' @param object a CopulaCenR object
#' @param ... further arguments
#' @importFrom stats BIC
#' @export
BIC.CopulaCenR <- function(object, ...) {

  # log(n)*k - 2*llk
  n <- nrow(object$indata1)
  k <- length(object$estimates)
  res <- -2 * object$llk + log(n) * k
  return(res)

}

