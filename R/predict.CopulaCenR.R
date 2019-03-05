#' Predictions from CopulaCenR regression models
#'
#' Predictions for new observations based on \code{ic_sp_copula}, \code{ic_par_copula} and \code{rc_par_copula}
#'
#' @name predict.CopulaCenR
#' @aliases predict.CopulaCenR
#' @param object a \code{CopulaCenR} object from \code{ic_sp_copula},
#' \code{ic_par_copula} and \code{rc_par_copula}
#' @param newdata a data frame (see details)
#' @param type \code{"lp"} for linear predictors or
#' \code{"survival"} for marginal and joint survival probabilities
#' @param ... further arguments
#' @importFrom stats predict
#' @importFrom caret dummyVars
#' @export
#'
#' @details
#'
#' For the \code{newdata}, when \code{type = "survival"}, it must be a data frame with columns
#' \code{id} (subject id), \code{ind} (1,2 for two margins),
#' \code{time} (to be evaluted) and \code{covariates};
#' when \code{type = "lp"}, the newdata needs to have \code{id},
#' \code{ind} and \code{covariates}, but \code{time} is not needed. \cr
#'
#' When the argument \code{type = "lp"}, it gives a linear predictor for
#' each margin (i.e., log hazards ratio in the proportional hazards model,
#' log proportional odds in the proportional odds model). \cr
#'
#' When the argument \code{type = "survival"}, the marginal and joint survival values
#' will be evaluated at the given time points in the \code{newdata}. \cr
#'
#' @return If \code{type = "lp"}, it returns a data frame with \code{id},
#' \code{lp1} (linear predictor for margin 1), \code{lp2}.
#' If \code{type = "survival"}, it returns a data frame with \code{id},
#' \code{t1} (evaluated times for the margin 1), \code{t2},
#' \code{S1} (predicted marginal survival probabilities for margin 1),
#' \code{S2} and
#' \code{S12} (the predicted joint survival probabilities at \code{t1, t2})
#' @examples
#' data(AREDS)
#' # fit a Copula2-Sieve model
#' copula2_sp <- ic_sp_copula(data = AREDS, copula = "Copula2",
#'               l = 0, u = 15, m = 3, r = 3,
#'               var_list = c("ENROLLAGE","rs2284665","SevScaleBL"))
#' # Predicted probabilities for newdata
#' newdata = data.frame(id = rep(1:3, each=2), ind = rep(c(1,2),3),
#'                      time = c(2,3,5,6,7,8),
#'                     SevScaleBL = rep(3,6),
#'                     ENROLLAGE = rep(60,6),
#'                     rs2284665 = c(0,0,1,1,2,2))
#' output <- predict(object = copula2_sp, newdata = newdata)


predict.CopulaCenR <- function(object, newdata, type = "lp", ...) {


  # first screen the inputs #
  if (class(object) != "CopulaCenR") {
    stop('object must be a CopulaCenR class object')
  }

  if (type == "survival") {
    if (!is.data.frame(newdata) |
        !"id" %in% colnames(newdata) |
        !"ind" %in% colnames(newdata) |
        !"time" %in% colnames(newdata)) {
      stop('when type is "survival", newdata must be a data frame with columns id, ind, time and covariates')
    }

  } else {
    if (!is.data.frame(newdata) |
        !"id" %in% colnames(newdata) |
        !"ind" %in% colnames(newdata)) {
      stop('when type is "lp", newdata must be a data frame with columns id, ind and covariates')
    }
  }

  # generate dummy variables
  tmp <- dummyVars(~ ., data = newdata, fullRank = T, sep = "")
  newdata <- data.frame(predict(tmp, newdata = newdata))
  # sort by id and ind
  newdata <- newdata[order(newdata$id, newdata$ind), ]
  id <- unique(newdata$id)

  if (type == "survival") {

    t1 <- newdata$time[newdata$ind == 1]
    t2 <- newdata$time[newdata$ind == 2]

    get_margin <- function(object, newdata, evalTimes1, evalTimes2) {
      internal_predict.CopulaCenR(object = object,
                                  class = "marginal",
                                  newdata = newdata,
                                  evalTimes1 = evalTimes1,
                                  evalTimes2 = evalTimes2)
    }

    get_joint <- function(object, newdata, evalTimes1, evalTimes2) {
      internal_predict.CopulaCenR(object = object,
                                  class = "joint",
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
    marginal <- data.frame(matrix(unlist(marginal), nrow = length(marginal),
                                  byrow = T))
    joint <- unlist(lapply(1:length(id),
                           function(x) {
                             tmp <- get_joint(object = object,
                                              newdata = newdata[newdata$id == id[x], ],
                                              evalTimes1 = t1[x],
                                              evalTimes2 = t2[x])
                             return(tmp$surv2)
                             }
                           )
                    )

    output <- data.frame(id, t1, t2, marginal, joint)
    colnames(output) <- c("id", "t1", "t2", "S1", "S2", "S12")
  } else {
    # lp1, lp2
    beta <- object$summary[object$var_list, "estimate"]
    x1 <- newdata[newdata$ind == 1, object$var_list]
    x2 <- newdata[newdata$ind == 2, object$var_list]
    x1 <- as.matrix(x1, ncol = length(object$var_list))
    x2 <- as.matrix(x2, ncol = length(object$var_list))

    lp1 <- x1 %*% beta
    lp2 <- x2 %*% beta

    output <- data.frame(id, lp1, lp2)
    colnames(output) <- c("id", "lp1", "lp2")
  }


  return(output)
}


