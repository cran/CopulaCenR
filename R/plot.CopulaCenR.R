#' Plotting for CopulaCenR fits
#'
#' Plotting for CopulaCenR fits from  \code{ic_spTran_copula}, \code{rc_spCox_copula},
#' \code{ic_par_copula} and \code{rc_par_copula}.
#'
#' @name plot.CopulaCenR
#' @aliases plot.CopulaCenR
#' @param x an object of \code{ic_spTran_copula} or \code{rc_spCox_copula}
#' or \code{ic_par_copula} or \code{rc_par_copula}
#' @param y new data frame with colname names \code{id}, \code{ind} and \code{covariate}
#' @param newdata new data frame (ignored if \code{y} is included)
#' @param class one of "joint", "conditional" or "marginal"
#' @param evalPoints number of time points to be evaluated in both margins;
#' default is 50
#' @param evalTimes1 a vector of times for margin 1 to be evaluated;
#' default is NULL;
#' will override evalPoints if non-NULL
#' @param evalTimes2 a vector of times for margin 2 to be evaluated
#' @param plot_margin for \code{class = "marginal"} only; indicator of which margin to plot
#' (either 1 or 2); default is 1 for margin 1
#' @param cond_margin for \code{class = "conditional"} only; indicator of the margin
#' where event has occurred (either 1 or 2); default is 2 for margin 2
#' @param cond_time for \code{class = "conditional"} only; the time
#' by which event has occurred in the margin indicated by \code{cond_margin};
#' must be smaller than the largest observed time
#' @param type type of plot with default \code{type = "l"}.
#' @param xlab a title for the x axis.
#' @param ylab a title for the x axis.
#' @param cex.main cex for main.
#' @param cex.lab cex for lab.
#' @param cex.axis cex for axis.
#' @param legend whether to show legend with default \code{legend = TRUE}.
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
#' @details
#'
#' y must be a data frame with columns \code{id} (subject id),
#' \code{ind} (1,2 for two margins) and \code{covariates}. \cr
#' The argument class determines the plot:
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
#' If \code{class} \code{= "conditional"}, one needs to specify the margin
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
#'
#' @return a 3D joint survival distribution plot if \code{class = "joint"};
#' a 2D survival distribution plot if \code{class} = \code{"marginal"} or \code{"conditional"}.
#'
#' @examples
#' data(AREDS)
#' # fit a Copula2-Sieve model
#' copula2_sp <- ic_spTran_copula(data = AREDS, copula = "Copula2",
#'               l = 0, u = 15, m = 3, r = 3,
#'               var_list = c("ENROLLAGE","rs2284665","SevScaleBL"))
#' newdata = data.frame(id = rep(1, each=2), ind = rep(c(1,2),1),
#'                      SevScaleBL = rep(3,2), ENROLLAGE = rep(60,2),
#'                      rs2284665 = c(0,0))
#' # Plot joint survival probabilities
#' plot(x = copula2_sp, class = "joint", newdata = newdata)
#'
#' # Plot conditional survival probabilities
#' plot(x = copula2_sp, class = "conditional", newdata = newdata,
#'      cond_margin = 2, cond_time = 5, ylim = c(0.25,1),
#'      ylab = "Conditional Survival Probability")
#'
#' # Plot marginal survival probabilities
#' plot(x = copula2_sp, class = "marginal", newdata = newdata,
#'      plot_margin = 1, ylim = c(0.6,1),
#'      ylab = "Marginal Survival Probability")


plot.CopulaCenR <- function(x, y, class = "joint", newdata, evalPoints = 50,
                            evalTimes1 = NULL, evalTimes2 = NULL,
                            plot_margin = 1, cond_time = NULL, cond_margin = 2,
                            type = "l", xlab = "years",
                            ylab = "survival probability",
                            cex.main = 1.4, cex.lab = 1.4, cex.axis = 1.4,
                            legend = TRUE, ...) {

  object <- x
  if (missing(y)) {y <- newdata}
  newdata <- y

  if (is.null(evalTimes2)) {
    evalTimes2 <- evalTimes1
  }

  if (is.null(evalTimes1)) {
    evalTimes1 <- evalTimes2
  }


  # first screen the inputs #
  # if (class(object) != "CopulaCenR") {
  if (isFALSE(inherits(object, what = "CopulaCenR"))) {
    stop('object must be a CopulaCenR class object')
  }

  if (!class %in% c("joint","conditional","marginal")) {
    stop('class must be one of joint, conditional and marginal')
  }

  if (!is.data.frame(newdata) |
      !"id" %in% colnames(newdata) |
      !"ind" %in% colnames(newdata)) {
    stop('newdata must be a data frame with columns id, ind and covariates')
  }


  # generate dummy variables
  tmp <- dummyVars(~ ., data = newdata, fullRank = T, sep = "")
  newdata <- data.frame(predict(tmp, newdata = newdata))

  # sort
  newdata <- newdata[order(newdata$id, newdata$ind), ]
  id <- unique(newdata$id)

  # plot a single line
  if (length(id) == 1) {

    if (class == "marginal") {
      output <- internal_predict.CopulaCenR(object = object, class = "marginal",
                                            evalPoints = evalPoints,
                                            evalTimes1 = evalTimes1,
                                            evalTimes2 = evalTimes2,
                                            newdata = newdata)
      if (plot_margin == 1) {
        plot(x = output$grid1, y = output$m1,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
      } else {
        plot(x = output$grid2, y = output$m2,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
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
        plot(x = output$grid1, y = output$condition,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
      } else {
        plot(x = output$grid2, y = output$condition,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
      }
    }

    else if (class == "joint") {
      output <- internal_predict.CopulaCenR(object = object, class = "joint",
                                            evalPoints = evalPoints,
                                            evalTimes1 = evalTimes1,
                                            evalTimes2 = evalTimes2,
                                            newdata = newdata)
      plot_ly(x = output$grid1, # along columns
              y = output$grid2, # along rows
              z = output$surv2,
              showscale = F, opacity = 0.9) %>%
        add_surface() %>%
        layout(title = "Joint Survival Probability",
               scene = list(xaxis = list(title = "time 1"),
                            yaxis = list(title = "time 2"),
                            zaxis = list(title = "", textangle = 90))
               )
    }

  } else { # more than one subject

    if (class == "marginal") {

      get_margin <- function(object, newdata, evalPoints, evalTimes1, evalTimes2) {
        internal_predict.CopulaCenR(object = object, class = "marginal",
                                    evalPoints = evalPoints,
                                    evalTimes1 = evalTimes1,
                                    evalTimes2 = evalTimes2,
                                    newdata = newdata)
      }
      output <- lapply(1:length(id),
                       function(x) {
                         tmp <- get_margin(object = object,
                                           evalPoints = evalPoints,
                                           newdata = newdata[newdata$id == id[x], ],
                                           evalTimes1 = evalTimes1,
                                           evalTimes2 = evalTimes1)
                         return(list(grid1 = tmp$grid1,
                                     m1 = tmp$m1,
                                     grid2 = tmp$grid2,
                                     m2=tmp$m2))
                         }
                       )

      if (plot_margin == 1) {
        plot(x = output[[1]]$grid1, y = output[[1]]$m1,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
        for(i in 2:length(output)) {
          lines(x = output[[i]]$grid1, y = output[[i]]$m1, lty = i)
        }

      } else {
        plot(x = output[[1]]$grid2, y = output[[1]]$m2,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
        for(i in 2:length(output)) {
          lines(x = output[[i]]$grid2, y = output[[i]]$m2, lty = i)
        }
      }
      if (isTRUE(legend)) {
        legend("bottomleft", legend = paste0("id: ",id),
               lty = 1:length(output), cex = 1.6)
        }
    }

    else if (class == "conditional") {
      get_conditional <- function(object, newdata, evalPoints,
                                  evalTimes1, evalTimes2,
                                  cond_time, cond_margin) {
        internal_predict.CopulaCenR(object = object, class = "conditional",
                                    evalPoints = evalPoints, evalTimes1 = evalTimes1,
                                    evalTimes2 = evalTimes2, newdata = newdata,
                                    cond_time = cond_time, cond_margin = cond_margin)
      }
      output <- lapply(1:length(id),
                       function(x) {
                         tmp <- get_conditional(object = object,
                                                evalPoints = evalPoints,
                                                newdata = newdata[newdata$id == id[x], ],
                                                evalTimes1 = evalTimes1,
                                                evalTimes2 = evalTimes1,
                                                cond_time = cond_time,
                                                cond_margin = cond_margin)
                         return(list(grid1 = tmp$grid1,
                                     grid2 = tmp$grid2,
                                     condition = tmp$condition))
                         }
                       )

      if (cond_margin == 2) {
        plot(x = output[[1]]$grid1, y = output[[1]]$condition,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
        for(i in 2:length(output)) {
          lines(x = output[[i]]$grid1, y = output[[i]]$condition, lty = i)
        }

      } else {
        plot(x = output[[1]]$grid2, y = output[[1]]$condition,
             type = type,
             xlab = xlab,
             ylab = ylab,
             cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
             ...)
        for(i in 2:length(output)) {
          lines(x = output[[i]]$grid2, y = output[[i]]$condition, lty = i)
        }
      }
      if (isTRUE(legend)) { legend("bottomleft", legend = paste0("id: ",id),
                                   lty = 1:length(output), cex = 1.6) }
    }

    else if (class == "joint") {
      get_joint <- function(object, newdata, evalPoints, evalTimes1, evalTimes2) {
        internal_predict.CopulaCenR(object = object, class = "joint",
                                    evalPoints = evalPoints, evalTimes1 = evalTimes1,
                                    evalTimes2 = evalTimes2, newdata = newdata)
      }
      output <- lapply(1:length(id),
                       function(x) {
                         tmp <- get_joint(object = object, evalPoints = evalPoints,
                                          newdata = newdata[newdata$id == id[x], ],
                                          evalTimes1 = evalTimes1, evalTimes2 = evalTimes1)
                         return(list(grid1 = tmp$grid1,
                                     grid2 = tmp$grid2,
                                     surv2 = tmp$surv2))
                         }
                       )

      plot3d <- plot_ly(x = output[[1]]$grid1, # along columns
                        y = output[[1]]$grid2, # along rows
                        z = output[[1]]$surv2,
                        showscale = F, opacity = 0.9) %>%
                        add_surface() %>%
                        layout(title = "Joint Survival Probability",
                               scene = list(xaxis = list(title = "time 1"),
                                            yaxis = list(title = "time 2"),
                                            zaxis = list(title = "", textangle = 90)))

      for(i in 2:length(output)) {
        plot3d <- add_surface(p = plot3d,
                              x = output[[i]]$grid1, # along columns
                              y = output[[i]]$grid2, # along rows
                              z = output[[i]]$surv2,
                              opacity = 0.9)
      }
      plot3d
   }
  }
}



