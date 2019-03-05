#' Calculate Kendall's tau
#'
#' To obtain Kendall's tau from copula parameter(s)
#'
#' @name tau_copula
#' @aliases tau_copula
#' @param eta copula parameter(s);
#' if \code{copula = "Coupla2"}, input \eqn{\alpha} and \eqn{\kappa}
#' @param copula specify the type of copula model
#' @importFrom copula archmCopula
#' @importFrom copula tau
#' @return Kendall's \eqn{\tau}
#' @export
#'
#' @details
#' The supported copula models are \code{"Clayton"}, \code{"Gumbel"}, \code{"Frank"},
#' \code{"AMH"}, \code{"Joe"} and \code{"Copula2"}.
#' The \code{"Copula2"} model is a two-parameter copula model that incorporates
#' \code{Clayton} and \code{Gumbel} as special cases. \cr
#'
#'
#' The Kendall's \eqn{\tau} formulas are list below:
#'
#' The Clayton copula Kendall's \eqn{\tau = \eta/(2+\eta)}.
#'
#' The Gumbel copula Kendall's \eqn{\tau = 1 - 1/\eta}.
#'
#' The Frank copula Kendall's \eqn{\tau = 1+4\{D_1(\eta)-1\}/\eta},
#' in which \eqn{D_1(\eta) = \frac{1}{\eta} \int_{0}^{\eta} \frac{t}{e^t-1}dt}.
#'
#' The AMH copula Kendall's \eqn{\tau =  1-2\{(1-\eta)^2 \log (1-\eta) + \eta\}/(3\eta^2)}.
#'
#' The Joe copula Kendall's \eqn{\tau = 1 - 4 \sum_{k=1}^{\infty} \frac{1}{k(\eta k+2)\{\eta(k-1)+2\}}}.
#'
#' The Two-parameter copula (\code{Copula2}) Kendall's \eqn{\tau = 1-2\alpha\kappa/(2\kappa+1)}. \cr
#'
#' @source
#' Ali MM, Mikhail NN, Haq MS (1978).
#' A Class of Bivariate Distributions Including the Bi- variate Logistic.
#' \emph{Journal of Multivariate Analysis} doi:10.1016/0047-259X(78)90063-5. \cr
#' Clayton DG (1978).
#' A Model for Association in Bivariate Life Tables and Application in
#' Epidemiological Studies of Familial Tendency in Chronic Disease Incidence.
#' \emph{Biometrika} doi:10.2307/2335289. \cr
#' Gumbel EJ (1960).
#' Bivariate Exponential Distributions.
#' \emph{Journal of the American Statistical Association}
#' doi:10.2307/2281591. \cr
#' Joe H (1993).
#' Parametric Families of Multivariate Distributions with Given Margins.
#' \emph{Journal of Multivariate Analysis}
#' doi:10.1006/jmva.1993.1061. \cr
#' Joe H (1997).
#' Multivariate Models and Dependence Concepts.
#' \emph{Chapman & Hall, London}. \cr
#' Frank MJ (1979).
#' On the Simultaneous Associativity of \eqn{F(x, y)}
#' and \eqn{x + y - F(x, y)}.
#' \emph{Aequationes Mathematicae}. \cr

tau_copula <- function(eta, copula){

  if (tolower(copula) == "copula2"){
    alpha <- eta[1]
    kappa <- eta[2]
    output <- 1 - 2 * alpha * kappa/(1 + 2 * kappa)
  }

  if (tolower(copula) != "copula2") {
    output <- tau(archmCopula(tolower(copula), param = eta, dim = 2))
  }

  return(output)
}
