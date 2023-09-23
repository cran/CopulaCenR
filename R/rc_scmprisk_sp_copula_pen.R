#' Penalized copula regression models with Cox semiparametric margins for semi-competing risk data under right-censoring.
#'
#' @description Fits a penalized copula model with Cox semiparametric margins for semi-competing risk data under right-censoring.
#'
#' @name rc_scmprisk_sp_copula_pen
#' @aliases rc_scmprisk_sp_copula_pen
#' @param data a data frame; must have \code{id} (subject id), \code{ind} (1,2 for two margins),
#' \code{obs_time}, \code{status} (0 for right-censoring, 1 for event).
#' @param var_list1 the list of covariates to be fitted into the copula model for non-terminal event.
#' @param var_list2 the list of covariates to be fitted into the copula model for terminal event.
#' @param m1 integer, degree of Bernstein polynomials for non-terminal event; default is 3.
#' @param m2 integer, degree of Bernstein polynomials for terminal event; default is 3.
#' @param initial a vector of initial values for all parameters, including phi1 and phi2
#' (two Bernstein polynomials parameters), beta1 and beta2 (regression coefficients in Cox margins),
#' \code{eta} (copula parameters).
#' @param a Tuning parameters in penalty function for non-terminal event.
#' @param b Tuning parameters in penalty function for terminal event.
#' @param pen1 Types of penalty function for non-terminal event, including \code{"NOPEN"}, \code{"RIDGE"},
#' \code{"BAR"}, \code{"LASSO"}, \code{"MCP"}, \code{"SCAD"}.
#' @param pen2 Types of penalty function for terminal event, including \code{"NOPEN"}, \code{"RIDGE"},
#' \code{"BAR"}, \code{"LASSO"}, \code{"MCP"}, \code{"SCAD"}.
#' @param copula Types of copula, including \code{"Clayton"}, \code{"Gumbel"}.
#' @importFrom stats optim
#' @importFrom stats nlm
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @return parameter estimates and BIC
#' @export
#' @source
#' Tao Sun, Weijie Liang, Gongzi Zhang, Danhui Yi, Ying Ding, Lihai Zhang (2023+). Penalized semiparametric
#' copula method for semi-competing risks data: Application to hip fracture in elderly. JRSSC.
#' @examples
#'\dontrun{
#' data("data_sim_scmprisk_vs")
#' var_list = paste0('x',seq(1,5,1))
#' phi1_ini = c(0,0,0,0)
#' phi2_ini = c(0,0,0,0)
#' beta1_ini = c(0,0,0,0,0)
#' beta2_ini = c(0,0,0,0,0)
#' eta_ini = 1
#' initial = c(phi1_ini,phi2_ini,beta1_ini, beta2_ini, eta_ini)
#' # obtain initial parameter estimates by ridge penalty
#' fit0 = rc_scmprisk_sp_copula_pen(data=data_sim_scmprisk_vs,
#'                                  var_list1 = var_list, var_list2 = var_list,
#'                                  m1=3, m2=3, initial=initial,
#'                                  a=0.1, b=0.1, pen1='RIDGE', pen2='RIDGE',
#'                                  copula='Clayton')
#' est_ini = c(fit0$Est_phi1,fit0$Est_phi2,fit0$Est_beta1,fit0$Est_beta2,fit0$Est_eta)
#'
#' fit = rc_scmprisk_sp_copula_pen(data=data_sim_scmprisk_vs,
#'                                 var_list1 = var_list, var_list2 = var_list,
#'                                 m1=3, m2=3, initial=est_ini,
#'                                 a=0.2, b=0.2, pen1='MCP', pen2='MCP',
#'                                 copula='Clayton')
#' fit$Est_beta1
#' fit$Est_beta2
#'}


rc_scmprisk_sp_copula_pen <- function(data, var_list1, var_list2, m1, m2, initial, a, b, pen1, pen2, copula){

  eps0 <- 1e-3  # threshold
  eps <- 1e-200
  indata1 <- data[data[,'ind']==1,]
  indata2 <- data[data[,'ind']==2,]
  n <- nrow(indata1)
  x1 <- as.matrix(indata1[,var_list1], nrow = dim(indata1)[1])
  x2 <- as.matrix(indata2[,var_list2], nrow = dim(indata2)[1])

  # Bernstein polynomial parameters
  l1 = 0
  u1 = max(indata1$obs_time)+1
  m1 = m1
  l2 = 0
  u2 = max(indata2$obs_time)+1
  m2 = m2

  data_processed1 <- data_process_scmprisk_rc_sp(indata1, var_list1, l=l1, u=u1, m=m1)
  data_processed2 <- data_process_scmprisk_rc_sp(indata2, var_list2, l=l2, u=u2, m=m2)

  # BP
  b1 <- data_processed1$b
  b2 <- data_processed2$b

  # BP derivatives
  b1_d <- data_processed1$b_d
  b2_d <- data_processed2$b_d

  # initial values
  phi1 <- initial[1:(m1+1)]
  phi2 <- initial[(m1+2):(m1+2+m2)]
  beta1 <- initial[(m1+3+m2):(m1+2+m2+ncol(x1))]
  beta2 <- initial[(m1+3+m2+ncol(x1)):(m1+2+m2+ncol(x1)+ncol(x2))]
  beta1_i <- beta1
  beta2_i <- beta2
  eta <- initial[length(initial)]

  act.set1 <- c(1:ncol(x1))
  act.set2 <- c(1:ncol(x2))
  if(pen1=='NOPEN') act.set1 <- which(abs(beta1)>eps0)
  if(pen2=='NOPEN') act.set2 <- which(abs(beta2)>eps0)

  # initial estimates \phi1 and \phi2
  fit1_0 <- optim(par = phi1,
                  scmprisk_log_lik_baseline1,
                  fitted = beta1,
                  method = "Nelder-Mead", control=c(maxit=100),
                  hessian = FALSE, x1 = x1, indata1 = indata1,
                  b1=b1, b1_d=b1_d, m1=m1)
  fit1_0$convergence
  phi1 <- fit1_0$par

  fit2_0 <- optim(par = phi2,
                  scmprisk_log_lik_baseline2,
                  fitted = beta2,
                  method = "Nelder-Mead", control=c(maxit=100),
                  hessian = FALSE, x2 = x2, indata2 = indata2,
                  b2=b2, b2_d=b2_d, m2=m2)
  fit2_0$convergence
  phi2 <- fit2_0$par

  iter<-0
  dd <- NULL
  repeat{

    # obtain \phi1 and \phi2
    fit <- nlm(p = c(phi1,phi2), rc_scmprisk_copula_loglik_sieve_baseline,
               fitted = c(log(eta),beta1,beta2),
               x1 = x1, x2 = x2, indata1 = indata1, indata2 = indata2,
               b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2,
               iterlim = 20, hessian = F,
               copula = copula)

    fit$code
    phi1_i <- fit$estimate[1:(m1+1)]
    phi2_i <- fit$estimate[(m1+2):(m1+2+m2)]

    # obtain \beta2
    if(pen2=='NOPEN'){

      foreach(dd=act.set2, .combine='c') %do% {

        dLc_b2 <- d_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b2 <- dd_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta2_i[dd] <- beta2_i[dd] - dLc_b2/ddLc_b2

      }
    }


    if(pen2=='RIDGE'){

      foreach(dd=act.set2, .combine='c') %do% {

        dLc_b2 <- d_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b2 <- dd_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta2_i[dd] <- beta2_i[dd] - (dLc_b2-2*b*beta2_i[dd])/(ddLc_b2-2*b)

      }
    }

    if(pen2=='BAR'){

      foreach(dd=act.set2, .combine='c') %do% {

        dLc_b2 <- d_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b2 <- dd_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta2_i[dd] <- beta2_i[dd] - (dLc_b2 - 2*b*beta2_i[dd]/(beta2_i[dd]^2+eps))/(ddLc_b2 - 2*b/(beta2_i[dd]^2+eps))
      }
    }

    if(pen2=='LASSO'){

      foreach(dd=act.set2, .combine='c') %do% {

        dLc_b2 <- d_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b2 <- dd_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta2_i[dd] <- (beta2_i[dd]*ddLc_b2-dLc_b2)/(ddLc_b2-n*b/(abs(beta2_i[dd])+eps))
      }
    }


    if(pen2=='MCP'){

      foreach(dd=act.set2, .combine='c') %do% {

        dLc_b2 <- d_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b2 <- dd_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta2_i[dd] <- (beta2_i[dd]*ddLc_b2-dLc_b2)/(ddLc_b2-n*dMCP2(abs(beta2_i[dd]),1.1,b)/(abs(beta2_i[dd])+eps))
      }
    }

    if(pen2=='SCAD'){

      foreach(dd=act.set2, .combine='c') %do% {

        dLc_b2 <- d_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b2 <- dd_logL_beta2(p=beta2_i,fitted=c(log(eta),phi1_i,phi2_i,beta1),x1=x1,x2=x2,indata1=indata1,indata2=indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta2_i[dd] <- (beta2_i[dd]*ddLc_b2-dLc_b2)/(ddLc_b2-n*dSCAD2(abs(beta2_i[dd]),3.7,b)/(abs(beta2_i[dd])+eps))

      }
    }


    act.set2 <- which(abs(beta2_i) > eps0)
    beta2_i[abs(beta2_i) < eps0] <- 0


    # obtain \beta1
    if(pen1=='NOPEN'){

      foreach(dd=act.set1, .combine='c') %do% {

        dLc_b1 <- d_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b1 <- dd_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta1_i[dd] <- beta1_i[dd] - dLc_b1/ddLc_b1

      }
    }


    if(pen1=='RIDGE'){

      foreach(dd=act.set1, .combine='c') %do% {

        dLc_b1 <- d_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b1 <- dd_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta1_i[dd] <- beta1_i[dd] - (dLc_b1 - 2*a*beta1_i[dd])/(ddLc_b1 - 2*a)

      }
    }

    if(pen1=='BAR'){

      foreach(dd=act.set1, .combine='c') %do% {

        dLc_b1 <- d_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b1 <- dd_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta1_i[dd] <- beta1_i[dd] - (dLc_b1 - 2*a*beta1_i[dd]/(beta1_i[dd]^2+eps))/(ddLc_b1 - 2*a/(beta1_i[dd]^2+eps))
      }
    }

    if(pen1=='LASSO'){

      foreach(dd=act.set1, .combine='c') %do% {

        dLc_b1 <- d_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b1 <- dd_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta1_i[dd] <- (beta1_i[dd]*ddLc_b1-dLc_b1)/(ddLc_b1-n*a/(abs(beta1_i[dd])+eps))
      }
    }


    if(pen1=='MCP'){

      foreach(dd=act.set1, .combine='c') %do% {

        dLc_b1 <- d_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b1 <- dd_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta1_i[dd] <- (beta1_i[dd]*ddLc_b1-dLc_b1)/(ddLc_b1-n*dMCP1(abs(beta1_i[dd]),1.1,a)/(abs(beta1_i[dd])+eps))
      }
    }

    if(pen1=='SCAD'){

      foreach(dd=act.set1, .combine='c') %do% {

        dLc_b1 <- d_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)
        ddLc_b1 <- dd_logL_beta1(p=beta1_i,fitted=c(log(eta),phi1_i,phi2_i,beta2_i),x1=x1,x2=x2,indata1=indata1,indata2=indata2, b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2, copula=copula,dd=dd)

        beta1_i[dd] <- (beta1_i[dd]*ddLc_b1-dLc_b1)/(ddLc_b1-n*dSCAD1(abs(beta1_i[dd]),3.7,a)/(abs(beta1_i[dd])+eps))

      }
    }


    act.set1 <- which(abs(beta1_i) > eps0)
    beta1_i[abs(beta1_i) < eps0] <- 0


    # obtain \eta
    phi <- nlm(rc_scmprisk_copula_loglik_sieve_eta, p = log(eta),
               fitted = c(phi1_i,phi2_i,beta1_i,beta2_i),
               x1 = x1, x2 = x2, indata1 = indata1, indata2 = indata2,
               b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2,
               iterlim = 1, stepmax = 1e-3, hessian = F,
               copula = copula)

    phi$code
    eta_i <- exp(phi$estimate)

    # evaluate the convergence
    theta=as.matrix(c(phi1, phi2, beta1, beta2, eta))
    theta_i=as.matrix(c(phi1_i, phi2_i, beta1_i, beta2_i, eta_i))

    if((mean(abs(theta_i-theta))<1e-3) | (iter>300)) (break) #stopping rule

    iter <- iter+1

    # update all parameters
    beta2 <- beta2_i
    beta1 <- beta1_i
    eta <- eta_i
    phi1 <- phi1_i
    phi2 <- phi2_i

    # cat(iter, '\n')
  }

  # After convergence, calculate log-likelihood and BIC

  # log-likelihood
  loglike <- rc_scmprisk_copula_loglik(beta = c(beta1_i,beta2_i), fitted=c(log(eta_i),phi1_i,phi2_i),x1 = x1, x2 = x2, indata1 = indata1,indata2 = indata2,b1=b1, b2=b2, b1_d=b1_d, b2_d=b2_d, m1=m1, m2=m2,copula = copula)  # ???????????????????????????

  # BIC
  df<-sum(beta2_i!=0)+sum(beta1_i!=0)
  bic <- (-2)*loglike+log(n)*df

  res<-list(iter, bic, loglike, phi2_i, phi1_i, beta2_i, beta1_i, eta_i)

  names(res)[[1]]<-"iter"
  names(res)[[2]]<-"BIC"
  names(res)[[3]]<-"loglike"
  names(res)[[4]]<-"Est_phi2"
  names(res)[[5]]<-"Est_phi1"
  names(res)[[6]]<-"Est_beta2"
  names(res)[[7]]<-"Est_beta1"
  names(res)[[8]]<-"Est_eta"

  return(res)
}
