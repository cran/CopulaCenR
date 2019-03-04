#' Plotting for CopulaCenR fits
#'
#' Plotting for CopulaCenR fits from ic_sp_copula, ic_par_copula and rc_par_copula
#'
#' @name plot.CopulaCenR
#' @aliases plot.CopulaCenR
#' @param x an object of ic_sp_copula or ic_par_copula or rc_par_copula
#' @param y new data frame with colname names id, ind and covariate names
#' @param newdata new data frame (ignored if \code{y} is included)
#' @param class one of "joint", "conditional" or "marginal"
#' @param evalPoints number of time points to be evaluated; default is 50
#' @param evalTimes a vector of times to be evaluated within the observed time range; default is NULL; will override evalPoints if non-NULL; if class is "conditional", the evaluation times are evalTimes + cond_time
#' @param plot_margin for class = "marginal" only; indicator of which margin to plot (either 1 or 2); default is 1 for ind = 1
#' @param cond_margin for class = "conditional" only; indicator of the margin where event has occurred (either 1 or 2); default is 2 for ind = 2
#' @param cond_time for class = "conditional" only; the time by which event has occurred in the margin indicated by cond_margin; must be smaller than the largest observed time
#' @param ... further arguments
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
#' y must be a data frame with columns id (subject id), ind (1,2 for two margins) and covariates. \cr
#' The argument class determines the plot output: "joint" for joint survival probabilities,
#' "conditional" for conditional probabilities and "marginal" or marginal probabilities. The function
#' evaluates on a series of time points (given by evalPoints or evalTimes; evalTimes will override evalPoints). \cr
#'
#' If class = "conditional", one needs to specify the margin that has the event (by cond_margin) and time when the event has occurred (by cond_time).
#' For example, if cond_margin = 2 and cond_time = 5, then the function produces the conditional survival probability (after time 5) in margin 1 given that
#' margin 2 has got an event by time 5. This measurement is useful for predicting the second event given the first event has occurred. \cr
#'
#' If class = "marginal", one needs to specify which margin to plot through the argument plot_margin.
#'
#' @return a 3D joint survival distribution plot if class is "joint"; a 2D survival distribution plot if class is "marginal" or "conditional"
#'
#' @examples
#' data(AREDS)
#' # fit a Copula2-Sieve model
#' copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2", l = 0, u = 15,
#'               m = 3, r = 3, var_list = c("ENROLLAGE","rs2284665","SevScaleBL"),
#'               iter = 300, stepsize = 1e-6, method = "Newton")
#' newdata = data.frame(id = rep(1:3, each=2), ind = rep(c(1,2),3),
#'                      SevScaleBL = rep(3,6), ENROLLAGE = rep(60,6),
#'                      rs2284665 = c(0,0,1,1,2,2))
#' # Plot marginal survival probabilities
#' plot(x = copula2_sp, class = "marginal", newdata = newdata[newdata$id==1,],
#'      plot_margin = 1, ylim = c(0.6,1), type = "l", xlab = "years",
#'      ylab = "Marginal Survival Probability", cex.main = 1.3, cex.lab = 1.3,
#'      cex.axis = 1.3)
#' lines(x = copula2_sp, class = "marginal", newdata = newdata[newdata$id==2,],
#'       plot_margin = 1, lty = 2)
#' legend("bottomleft", c("GG","GT"), lty = c(1,2), lwd = c(1,1),
#'        col = c("black","black"), cex = 1.3)
#'
#' # Plot conditional survival probabilities
#' plot(x = copula2_sp, class = "conditional", newdata = newdata[newdata$id==1,],
#'      cond_margin = 2, cond_time = 5, ylim = c(0.25,1), type = "l",
#'      xlab = "years", ylab = "Conditional Survival Probability",
#'      cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3)
#' lines(x = copula2_sp, class = "conditional", newdata = newdata[newdata$id==2,],
#'      cond_margin = 2, cond_time = 5, lty = 2)
#' legend("bottomleft", c("GG","GT"), lty = c(1,2), lwd = c(1,1),
#'       col = c("black","black"), cex = 1.3)
#'
#' # Plot joint survival probabilities
#' plot3d <- plot(x = copula2_sp, class = "joint",
#'                newdata = newdata[newdata$id==1,])
#' plot3d <- lines(x = copula2_sp, class = "joint",
#'                 newdata = newdata[newdata$id==2,], plotly_object = plot3d)


plot.CopulaCenR <- function(x, y, class = "joint", newdata, evalPoints = 50, evalTimes = NULL, plot_margin = 1, cond_time = NULL, cond_margin = 2, ...) {

  object <- x
  if(missing(y)) y <- newdata
  newdata <- y

  # first screen the inputs #
  if (class(object) != "CopulaCenR") stop('object must be a CopulaCenR class object')
  if (!class %in% c("joint","conditional","marginal")) stop('class must be one of joint, conditional and marginal')
  if (!is.data.frame(newdata) | !"id" %in% colnames(newdata) | !"ind" %in% colnames(newdata)) stop('newdata must be a data frame with columns id, ind and var_list')


    if (class == "marginal") {
      output <- predict(object = object, class = "marginal", evalPoints = evalPoints, evalTimes = evalTimes, newdata = newdata)
      if (plot_margin == 1) {
        plot(x = output$grid1, y = output$m1, ...)
      } else {
        plot(x = output$grid2, y = output$m2, ...)
      }
    }

    else if (class == "conditional") {
      output <- predict(object = object, class = "conditional", evalPoints = evalPoints, evalTimes = evalTimes, newdata = newdata, cond_time = cond_time, cond_margin = cond_margin)
      if (cond_margin == 2) {
        plot(x = output$grid1, y = output$condition, ...)
      } else {
        plot(x = output$grid2, y = output$condition, ...)
      }
    }

    else if (class == "joint") {
      output <- predict(object = object, class = "joint", evalPoints = evalPoints, evalTimes = evalTimes, newdata = newdata)
      plot_ly(x = output$grid1, # along columns
              y = output$grid2, # along rows
              z = output$surv2,
              showscale = F, opacity = 0.9) %>%
        add_surface() %>%
        layout(title = "Joint Survival Probability", scene = list(xaxis = list(title = "time 1"), yaxis = list(title = "time 2"), zaxis = list(title = "",textangle = 90)))
    }


}


