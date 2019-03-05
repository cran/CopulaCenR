#' Plotting for CopulaCenR fits
#'
#' Plotting for CopulaCenR fits from \code{ic_sp_copula}, \code{ic_par_copula} and \code{rc_par_copula}.
#'
#' @name lines.CopulaCenR
#' @aliases lines.CopulaCenR
#' @param x an object of \code{ic_sp_copula, ic_par_copula, rc_par_copula}
#' @param y new data frame with colname names \code{id}, \code{ind} and \code{covariate}
#' @param newdata new data frame (ignored if \code{y} is included)
#' @param class one of "joint", "conditional" or "marginal"
#' @param plotly_object only for \code{class = "joint"}, an object of \code{plot.CopulaCenR}
#' @param evalPoints number of time points to be evaluated; default is 50
#' @param evalTimes1 a vector of times for margin 1 to be evaluated;
#' default is NULL; will override evalPoints if non-NULL
#' @param evalTimes2 a vector of times for margin 2 to be evaluated
#' @param plot_margin for \code{class = "marginal"} only; indicator of which margin to plot
#' (either 1 or 2); default is 1 for margin 1
#' @param cond_margin for \code{class = "conditional"} only; indicator of the margin
#' where event has occurred (either 1 or 2); default is 2 for margin 2
#' @param cond_time for \code{class = "conditional"} only; the time
#' by which event has occurred in the margin indicated by cond_margin;
#' must be smaller than the largest observed time
#' @param ... further arguments
#' @importFrom caret dummyVars
#' @importFrom plotly plot_ly
#' @importFrom plotly add_surface
#' @importFrom plotly layout
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom stats predict
#' @import magrittr
#' @export
#'
#' @details
#'
#' y must be a data frame with columns \code{id} (subject id),
#' \code{ind} (1,2 for two margins) and \code{covariates}. \cr
#' The argument \code{class} determines the plot:
#' \code{"joint"} for joint survival probabilities,
#' \code{"conditional"} for conditional probabilities and
#' \code{"marginal"} for marginal probabilities.
#'
#' The function evaluates on a series of time points
#' (given by \code{evalPoints} or \code{evalTimes};
#' \code{evalTimes} will override \code{evalPoints}).
#' By default, the time points are automatically
#' selected by specifying the number of points (\code{evalPoints = 50}).
#' Users can also provide the specific time points through \code{evalTimes1} and
#' \code{evalTimes2} for the two margins, respectively.
#' When \code{class} \code{= "conditional"}, only \code{evalTimes1} is needed
#' and the evaluation times are actually \code{evalTimes1} plus \code{cond_time}. \cr
#'
#' If \code{class = "conditional"}, one needs to specify the margin
#' that has the event (by \code{cond_margin})
#' and time when the event has occurred (by \code{cond_time}).
#' For example, if \code{cond_margin = 2} and \code{cond_time = 5},
#' then the function produces the conditional survival probability
#' (after time 5) in margin 1 given that margin 2 has got an event by time 5.
#' This measurement is useful for predicting the second event
#' given the first event has occurred. See the example for details.\cr
#'
#' If \code{class = "marginal"}, one needs to specify which margin to plot
#' through the argument \code{plot_margin}. See the example for details. \cr
#'
#' If \code{class = "joint"}, one needs to include a \code{plot_ly} object
#' (from \code{plot.CopulaCenR} with \code{class} = \code{"joint"})
#' through the argument \code{plotly_object}. See the example for details.\cr
#'
#' @return a 3D joint survival distribution plot if \code{class = "joint"};
#' a 2D survival distribution plot if \code{class} \code{= "marginal"}
#' or \code{"conditional"}.
#'
#' @examples
#' data(AREDS)
#' # fit a Copula2-Sieve model
#' copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2",
#'               l = 0, u = 15, m = 3, r = 3,
#'               var_list = c("ENROLLAGE","rs2284665","SevScaleBL"))
#' newdata = data.frame(id = rep(1:3, each=2), ind = rep(c(1,2),3),
#'                      SevScaleBL = rep(3,6), ENROLLAGE = rep(60,6),
#'                      rs2284665 = c(0,0,1,1,2,2))
#' # Plot marginal survival probabilities
#' plot(x = copula2_sp, class = "marginal",
#'      newdata = newdata[newdata$id==1,],
#'      plot_margin = 1, ylim = c(0.6,1),
#'      ylab = "Marginal Survival Probability")
#' lines(x = copula2_sp, class = "marginal",
#'       newdata = newdata[newdata$id==2,],
#'       plot_margin = 1, lty = 2)
#' legend("bottomleft", c("id: 1","id: 2"), lty = c(1,2))
#'
#' # Plot conditional survival probabilities
#' plot(x = copula2_sp, class = "conditional",
#'      newdata = newdata[newdata$id==1,],
#'      cond_margin = 2, cond_time = 5, ylim = c(0.25,1),
#'      xlab = "years", ylab = "Conditional Survival Probability")
#' lines(x = copula2_sp, class = "conditional",
#'       newdata = newdata[newdata$id==2,],
#'      cond_margin = 2, cond_time = 5, lty = 2)
#' legend("bottomleft", c("GG","GT"), lty = c(1,2))
#'
#' # Plot joint survival probabilities
#' plot3d <- plot(x = copula2_sp, class = "joint",
#'                newdata = newdata[newdata$id==1,])
#' plot3d <- lines(x = copula2_sp, class = "joint",
#'                 newdata = newdata[newdata$id==2,], plotly_object = plot3d)


lines.CopulaCenR <- function(x, y, class = "joint", newdata, evalPoints = 50,
                             evalTimes1 = NULL, evalTimes2 = NULL,
                             plot_margin = 1, cond_time = NULL,
                             cond_margin = 2, plotly_object = NULL,
                             ...) {

  object <- x
  if(missing(y)) y <- newdata
  newdata <- y

  if (is.null(evalTimes2)) {evalTimes2 <- evalTimes1}
  if (is.null(evalTimes1)) {evalTimes1 <- evalTimes2}

  # first screen the inputs #
  if (class(object) != "CopulaCenR") {
    stop('object must be a CopulaCenR class object')
  }

  if (!class %in% c("joint","conditional","marginal")) {
    stop('class must be one of joint, conditional and marginal')
  }

  if (!is.data.frame(newdata) |
      !"id" %in% colnames(newdata) |
      !"ind" %in% colnames(newdata)) {
    stop('newdata must be a data frame with columns id, ind and var_list')
  }


  # generate dummy variables
  tmp <- dummyVars(~ ., data = newdata, fullRank = T, sep = "")
  newdata <- data.frame(predict(tmp, newdata = newdata))

  if (class == "marginal") {
    output <- internal_predict.CopulaCenR(object = object, class = "marginal",
                                          evalPoints = evalPoints,
                                          evalTimes1 = evalTimes1,
                                          evalTimes2 = evalTimes2,
                                          newdata = newdata)
    if (plot_margin == 1) {
      lines(x = output$grid1, y = output$m1, ...)
    } else {
      lines(x = output$grid2, y = output$m2, ...)
    }
  }


  else if (class == "conditional") {
    output <- internal_predict.CopulaCenR(object = object, class = "conditional",
                                          evalPoints = evalPoints,
                                          evalTimes1 = evalTimes1,
                                          evalTimes2 = evalTimes2,
                                          newdata = newdata,
                                          cond_time = cond_time,
                                          cond_margin = cond_margin)
    if (cond_margin == 2) {
      lines(x = output$grid1, y = output$condition, ...)
    } else {
      lines(x = output$grid2, y = output$condition, ...)
    }
  }

  else if (class == "joint") {
    output <- internal_predict.CopulaCenR(object = object, class = "joint",
                                          evalPoints = evalPoints,
                                          evalTimes1 = evalTimes1,
                                          evalTimes2 = evalTimes2,
                                          newdata = newdata)
    add_surface(p = plotly_object,
                x = output$grid1, # along columns
                y = output$grid2, # along rows
                z = output$surv2,
                opacity = 0.9)
  }

}


