#' Predictions from CopulaCenR regression models
#'
#' Predict survival distribution for a new observation from ic_sp_copula, ic_par_copula and rc_par_copula.
#'
#' @name predict.CopulaCenR
#' @aliases predict.CopulaCenR
#' @param object a CopulaCenR object from ic_sp_copula, ic_par_copula and rc_par_copula
#' @param class one of "joint", "conditional" or "marginal"
#' @param evalPoints number of time points to be evaluated; default is 50
#' @param evalTimes a vector of times to be evaluated within the observed time range; default is NULL; will override evalPoints if non-NULL; if class is "conditional", the evaluation times are evalTimes + cond_time
#' @param newdata a data frame with colname names id, ind and covariate names
#' @param cond_margin for class = "conditional" only; indicator of the margin for which event has occurred (either 1 or 2); default is 2 for ind = 2
#' @param cond_time for class = "conditional" only; the time by which event has occurred in the margin indicated by cond_margin; must be smaller than the largest observed time
#' @param ... further arguments
#' @importFrom stats predict
#' @export
#'
#' @details
#'
#' The newdata must be a data frame with columns id (subject id), ind (1,2 for two margins) and covariates. \cr
#' The argument class determines the prediction output: "joint" for joint survival probabilities,
#' "conditional" for conditional probabilities and "marginal" or marginal probabilities. The function
#' evaluates on a series of time points (given by evalPoints or evalTimes; evalTimes will override evalPoints). \cr
#'
#' If class = "conditional", one needs to specify the margin that has the event (by cond_margin) and time when the event has occurred (by cond_time).
#' For example, if cond_margin = 2 and cond_time = 5, then the function produces the conditional survival probability (after time 5) in margin 1 given that
#' margin 2 has got an event by time 5. This measurement is useful for predicting the second event given the first event has occurred. \cr
#'
#' @return if class is "marginal" or "conditional", return a vector of survival probabilities with time grids;
#' if class is "joint", returns  a matrix (named surv2, row for event 1, column for event 2) of joint survival probabilities based on the fitted model in object, grid1 (selected times for event 1) and grid2 (selected times for event 2)
#'
#' @examples
#' data(AREDS)
#' # fit a Copula2-Sieve model
#' copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2", l = 0, u = 15,
#'               m = 3, r = 3, var_list = c("ENROLLAGE","rs2284665","SevScaleBL"),
#'               iter = 300, stepsize = 1e-6, method = "Newton")
#' # Predicted probabilities for newdata
#' newdata = data.frame(id = rep(1, each=2), ind = c(1,2), SevScaleBL = c(6,8),
#'                      ENROLLAGE = c(67,67), rs2284665 = c(1,1))
#' joint <- predict(object = copula2_sp, class = "joint", newdata = newdata)
#' conditional <- predict(object = copula2_sp, class = "conditional",
#'                        newdata = newdata, cond_margin = 2, cond_time = 5)
#' marginal <- predict(object = copula2_sp, class = "marginal",
#'                     newdata = newdata)


predict.CopulaCenR <- function(object, class = "joint", newdata, evalPoints = 50, evalTimes = NULL, cond_time = NULL, cond_margin = 2,...) {


  # first screen the inputs #
  if (class(object) != "CopulaCenR") stop('object must be a CopulaCenR class object')
  if (!class %in% c("joint","conditional","marginal")) stop('class must be one of joint, conditional and marginal')
  if (!is.data.frame(newdata) | !"id" %in% colnames(newdata) | !"ind" %in% colnames(newdata)) stop('newdata must be a data frame with columns id, ind and var_list')

  newdata1 = newdata[newdata$ind==1, object$var_list]
  newdata2 = newdata[newdata$ind==2, object$var_list]

  newdata1 = as.numeric(newdata1)
  newdata2 = as.numeric(newdata2)


  if (is.null(evalTimes) & is.null(evalPoints)) {
    evalPoints = 50 # set default
  }

  if (is.null(evalTimes) & !is.null(evalPoints)) {
    # calculate grid length
    xmin_1 = min(c(object$indata1[,"Left"], object$indata1[,"Right"]))
    xmax_1 = max(c(object$indata1[,"Left"], object$indata1[!is.infinite(object$indata1[,"Right"]),"Right"]))
    xmin_2 = min(c(object$indata2[,"Left"], object$indata2[,"Right"]))
    xmax_2 = max(c(object$indata2[,"Left"], object$indata2[!is.infinite(object$indata2[,"Right"]),"Right"]))
    grid.length1 = (xmax_1 - xmin_1)/evalPoints
    grid.length2 = (xmax_2 - xmin_2)/evalPoints

    if (class == "marginal") {
      output1 <- m_copula(object = object, grid.length = grid.length1, newdata = newdata1)
      output2 <- m_copula(object = object, grid.length = grid.length2, newdata = newdata2)
      output <- list(grid1 = output1$grid, m1 = output1$m, grid2 = output2$grid, m2 = output2$m,...)
    }

    if (class == "conditional") {
      output <- cond_copula(object = object, grid.length = grid.length1, newdata1 = newdata1, newdata2 = newdata2, cond_time = cond_time, cond_margin = cond_margin)
    }

    if (class == "joint") {
      output <- surv2_copula(object = object, grid.length1 = grid.length1, grid.length2 = grid.length2, newdata1 = newdata1, newdata2 = newdata2)
    }

  }


  if (!is.null(evalTimes)) {

    grid.length1 = NULL
    grid.length2 = NULL

    if (class == "marginal") {
      output1 <- m_copula(object = object, grid.length = grid.length1, newdata = newdata1, evalTimes = evalTimes)
      output2 <- m_copula(object = object, grid.length = grid.length2, newdata = newdata2, evalTimes = evalTimes)
      output <- list(grid1 = output1$grid, m1 = output1$m, grid2 = output2$grid, m2 = output2$m,...)

    }

    if (class == "conditional") {
      output <- cond_copula(object = object, grid.length = grid.length1, newdata1 = newdata1, newdata2 = newdata2, cond_time = cond_time, cond_margin = cond_margin, evalTimes = evalTimes)
    }

    if (class == "joint") {
      output <- surv2_copula(object = object, grid.length1 = grid.length1, grid.length2 = grid.length2, newdata1 = newdata1, newdata2 = newdata2, evalTimes = evalTimes)
    }

  }



  return(output)
}


