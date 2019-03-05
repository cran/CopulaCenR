#' Fitted values from CopulaCenR regression models
#'
#' Fitted values based on models from \code{ic_sp_copula}, \code{ic_par_copula} and \code{rc_par_copula}.
#'
#' @name fitted.CopulaCenR
#' @aliases fitted.CopulaCenR
#' @param object a \code{CopulaCenR} object from \code{ic_sp_copula}, \code{ic_par_copula} and \code{rc_par_copula}
#' @param type \code{"lp"} for linear predictors or
#' \code{"survival"} for marginal and joint survival probabilities
#' @param ... further arguments
#' @importFrom stats predict
#' @export
#'
#' @details
#' When the argument \code{type = "lp"}, it gives a linear predictor for each margin
#' (i.e., log hazards ratio in the proportional hazards model,
#' log proportional odds in the proportional odds model). \cr
#'
#' When the argument \code{type = "survival"} and the fitted data is bivariate right-censored,
#' the marginal and joint survival values will be evaluated at the observed times.
#' For bivariate interval-censored, evaluation times are the interval middle points or
#' left bound (if right bound is infinity). \cr
#'
#' @return If \code{type = "lp"}, it returns a data frame with \code{id},
#' \code{lp1} (linear predictor for margin 1), \code{lp2}.
#' If \code{type = "survival"}, it returns a data frame with \code{id},
#' \code{t1} (evaluated times for the margin 1), \code{t2},
#' \code{S1} (predicted marginal survival probabilities for margin 1),
#' \code{S2} and \code{S12} (the predicted joint survival probabilities)
#'
#' @examples
#' data(AREDS)
#' # fit a Copula2-Sieve model
#' copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2",
#'               l = 0, u = 15, m = 3, r = 3,
#'               var_list = c("ENROLLAGE","rs2284665","SevScaleBL"))
#' output <- fitted(object = copula2_sp)


fitted.CopulaCenR <- function(object, type = "lp", ...) {


  # first screen the inputs #
  if (class(object) != "CopulaCenR") {
    stop('object must be a CopulaCenR class object')
  }

  # pre-process, sort
  y1 <- cbind(object$indata1$id, object$indata1$ind, object$indata1$obs_time,
              object$indata1$Left, object$indata1$Right, object$indata1$status)
  y2 <- cbind(object$indata2$id, object$indata2$ind, object$indata2$obs_time,
              object$indata2$Left, object$indata2$Right, object$indata2$status)
  newdata1 <- cbind(y1, object$x1)
  newdata2 <- cbind(y2, object$x2)
  newdata <- rbind(newdata1, newdata2)
  newdata <- data.frame(newdata)
  if (ncol(y1) == 5) {
    colnames(newdata) <- c("id","ind","Left","Right","status",object$var_list)
  }

  if (ncol(y1) == 4) {
    colnames(newdata) <- c("id","ind","obs_time","status",object$var_list)
  }

  newdata <- newdata[order(newdata$id, newdata$ind), ]

  # time to be evaluated
  newdata$new_time <- sapply(1:nrow(newdata),
                             function(x) {
                               tmp <- newdata[x, ]
                               tmp <- c(tmp$obs_time, tmp$Left, tmp$Right)
                               return(mean(tmp[is.finite(tmp)]))
                               }
                             )
  newdata$new_time[newdata$status == 0] <- sapply(which(newdata$status == 0),
                                                  function(x) {
                                                    tmp <- newdata[x, ]
                                                    tmp <- c(tmp$obs_time, tmp$Left)
                                                    return(mean(tmp[is.finite(tmp)]))
                                                    }
                                                  )
  t1 <- newdata$new_time[newdata$ind == 1]
  t2 <- newdata$new_time[newdata$ind == 2]
  id <- unique(newdata$id)

  # S1, S2, S12

  if (type == "survival") {
    get_margin <- function(object, newdata, evalTimes1, evalTimes2) {
      internal_predict.CopulaCenR(object = object, class = "marginal",
                                  newdata = newdata,
                                  evalTimes1 = evalTimes1,
                                  evalTimes2 = evalTimes2)
    }

    get_joint <- function(object, newdata, evalTimes1, evalTimes2) {
      internal_predict.CopulaCenR(object = object, class = "joint",
                                  newdata = newdata,
                                  evalTimes1 = evalTimes1,
                                  evalTimes2 = evalTimes2)
    }

    marginal <- lapply(1:length(id),
                       function(x) {
                         tmp <- get_margin(object = object,
                                           newdata = newdata[newdata$id == id[x], ],
                                           evalTimes1 = t1[x], evalTimes2 = t2[x])
                         return(c(tmp$m1, tmp$m2))
                         }
                       )

    marginal <- data.frame(matrix(unlist(marginal),
                                  nrow = length(marginal), byrow = T))
    joint <- unlist(lapply(1:length(id),
                           function(x) {
                             tmp <- get_joint(object = object,
                                              newdata = newdata[newdata$id == id[x], ],
                                              evalTimes1 = t1[x], evalTimes2 = t2[x])
                             return(tmp$surv2)
                             }
                           )
                    )

    output <- data.frame(id, t1, t2, marginal, joint)
    colnames(output) <- c("id", "t1", "t2", "S1", "S2", "S12")

  } else {
    # lp1, lp2
    beta <- object$summary[object$var_list, "estimate"]
    lp1 <- object$x1 %*% beta
    lp2 <- object$x2 %*% beta

    output <- data.frame(id, lp1, lp2)
    colnames(output) <- c("id", "lp1", "lp2")
  }

  return(output)
}


