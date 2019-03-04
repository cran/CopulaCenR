#' Printing outputs of an CopulaCenR object
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
    printCoefmat(x$summary, P.values = T, has.Pvalue = T)
    cat("Note: The Wald tests are testing whether each coefficient is 0","\n")
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    cat("Final AIC:  ",x$AIC,"\n")
    cat("Convergence:  ",x$code,"\n")
  } else {
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  semiparametric","\n")
    cat("\n")
    printCoefmat(x$summary, P.values = T, has.Pvalue = T)
    cat("Note: The Wald tests are testing whether each coefficient is 0","\n")
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    cat("Final AIC:  ",x$AIC,"\n")
    cat("Convergence:  ",x$code,"\n")
  }
}


#' Summarizing outputs of an CopulaCenR object
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

#' Print the summary of an CopulaCenR object
#' @name print.summary.CopulaCenR
#' @aliases print.summary.CopulaCenR
#' @param x a CopulaCenR object
#' @param ... further arguments
#' @importFrom stats printCoefmat
#' @export
print.summary.CopulaCenR <- function(x,...) {
  if (!is.null(x$m.dist)){
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  ",x$m.dist,"\n")
    cat("\n")
    printCoefmat(x$summary, P.values = T, has.Pvalue = T)
    cat("Note: The Wald tests are testing whether each coefficient is 0","\n")
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    cat("Final AIC:  ",x$AIC,"\n")
    cat("Convergence:  ",x$code,"\n")
  } else {
    cat("Copula:  ",x$copula,"\n")
    cat("Margin:  semiparametric","\n")
    cat("\n")
    printCoefmat(x$summary, P.values = T, has.Pvalue = T)
    cat("Note: The Wald tests are testing whether each coefficient is 0","\n")
    cat("\n")
    cat("Final llk:  ",x$llk,"\n")
    cat("Final AIC:  ",x$AIC,"\n")
    cat("Convergence:  ",x$code,"\n")
  }
}

#' the coefficient estimates of an CopulaCenR object
#' @name coef.CopulaCenR
#' @aliases coef.CopulaCenR
#' @param object a CopulaCenR object
#' @param ... further arguments
#' @export
coef.CopulaCenR <- function(object,...) {

  res = object$summary[,1]
  res

}


#' the log-likelihood of an CopulaCenR object
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

