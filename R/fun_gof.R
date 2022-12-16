# functions for RC

Marginal_KM <- function(data) {

  dat <- data
  dat1 <- dat[,c("event_time.L","status.L")]
  dat2 <- dat[,c("event_time.R","status.R")]
  colnames(dat1) = colnames(dat2) = c("event_time","status")

  KM1 <- survfit(Surv(event_time,status) ~ 1,data = dat1)
  KM2 <- survfit(Surv(event_time,status) ~ 1,data = dat2)
  S1 <- data.frame(x=KM1$time,y=KM1$surv,event=KM1$n.event)
  S2 <- data.frame(x=KM2$time,y=KM2$surv,event=KM2$n.event)

  S1$y[S1$y==0] <- min(S1$y[S1$y>0])
  S2$y[S2$y==0] <- min(S2$y[S2$y>0])
  S1$y <- (1-1/nrow(dat))*S1$y
  S2$y <- (1-1/nrow(dat))*S2$y

  S1 <- S1[match(dat1$event_time, S1$x),]
  S2 <- S2[match(dat2$event_time, S2$x),]

  return(list(S1=S1, S2=S2))

}


copula_log_lik_in_semi <- function(para, S1, S2, data1, data2, copula) {

  u1<-S1$y
  u2<-S2$y

  if (copula == "plackett") {
    eta <- para
    myCop <- plackettCopula(param=eta)
    c_val <- dCopula(cbind(u1, u2), myCop)
    c_u1_val <- derCOP(u1, u2, cop = PLcop, para = eta)
    c_u2_val <- derCOP2(u1, u2, cop = PLcop, para = eta)
    C_val <- pCopula(cbind(u1, u2), myCop)
  }

  if (copula == "gaussian") {
    eta <- para
    myCop <- normalCopula(param=eta, dim=2)
    c_val <- dCopula(cbind(u1, u2), myCop)
    c_u1_val <- cCopula(cbind(u1, u2), myCop)[,2]
    c_u2_val <- cCopula(cbind(u2, u1), myCop)[,2]
    C_val <- pCopula(cbind(u1, u2), myCop)
  }

  if (copula == "clayton") {
    eta <- para
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
  }

  if (copula == "gumbel") {
    eta <- para
    c_val<- gh_F(u1,u2,eta) * (1/(u1*u2)) * ((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2) * (log(u1)*log(u2))^(eta-1) +
      gh_F(u1,u2,eta) * (-1/(u1*u2)) * (1-eta) * ((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2) * (log(u1)*log(u2))^(eta-1)
    c_u1_val<- ((-log(u1))^(eta) + (-log(u2))^(eta))^(1/eta -1) * (-log(u1))^(eta-1) * (1/u1) * gh_F(u1,u2,eta)
    c_u2_val<- ((-log(u1))^(eta) + (-log(u2))^(eta))^(1/eta -1) * (-log(u2))^(eta-1) * (1/u2) * gh_F(u1,u2,eta)
    C_val<-gh_F(u1,u2,eta)
  }

  if (copula == "frank") {
    eta <- para
    c_val <- (-eta)*exp((-eta)*u1)*exp((-eta)*u2)/((exp((-eta))-1)+(exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)) -
      (-eta)*exp((-eta)*u1)*exp((-eta)*u2)*(exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)/((exp((-eta))-1)+(exp((-eta)*u1)-1)*(exp((-eta)*u2)-1))^2
    c_u1_val <- (exp((-eta)*u2)-1)*exp((-eta)*u1)/((1 + (exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)/(exp((-eta))-1))*(exp((-eta))-1))
    c_u2_val <- (exp((-eta)*u1)-1)*exp((-eta)*u2)/((1 + (exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)/(exp((-eta))-1))*(exp((-eta))-1))
    C_val <- (1/(-eta)) * log(1 + (exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)/(exp((-eta))-1))
  }


  term1 <- log(C_val)
  term1 <- ifelse((data1[,"status"] == 0) & (data2[,"status"] == 0), term1, 0)

  term2 <- log(c_u1_val)
  term2 <- ifelse((data1[,"status"] == 1) & (data2[,"status"] == 0), term2, 0)

  term3 <- log(c_u2_val)
  term3 <- ifelse((data1[,"status"] == 0) & (data2[,"status"] == 1), term3, 0)

  term4 <- log(c_val)
  term4 <- ifelse((data1[,"status"] == 1) & (data2[,"status"] == 1), term4, 0)

  logL<-(-1)*sum( term1 + term2 + term3 + term4 )
  return(logL)
}


MLE_copula_np_PMLE_eta <- function(data, copula) {

  dat <- data
  dat1 <- dat[,c("event_time.L","status.L")]
  dat2 <- dat[,c("event_time.R","status.R")]
  colnames(dat1) = colnames(dat2) = c("event_time","status")

  KM1 <- survfit(Surv(event_time,status) ~ 1,data = dat1)
  KM2 <- survfit(Surv(event_time,status) ~ 1,data = dat2)
  S1 <- approx(x = KM1$time, y = KM1$surv, xout = dat1$event_time, rule = 2)
  S2 <- approx(x = KM2$time, y = KM2$surv, xout = dat2$event_time, rule = 2)
  S1$y[S1$y==0] <- min(S1$y[S1$y>0])
  S2$y[S2$y==0] <- min(S2$y[S2$y>0])
  S1$y <- (1-1/nrow(dat))*S1$y
  S2$y <- (1-1/nrow(dat))*S2$y

  S1 <- data.frame(x=S1$x, y=S1$y)
  S2 <- data.frame(x=S2$x, y=S2$y)
  S1 <- S1[match(dat1$event_time, S1$x),]
  S2 <- S2[match(dat2$event_time, S2$x),]


  if (copula != "frank" & copula != "copula2" & copula != "gaussian" & copula != "custom") {

    tau0 <- cor(dat1$event_time[dat1$status==1 & dat2$status==1], dat2$event_time[dat1$status==1 & dat2$status==1], method="kendall")
    if (copula == "plackett") {
      eta0 <- copula::iTau(copula = plackettCopula(), tau = tau0)
    } else {
      eta0 <- copula::iTau(copula = copula::archmCopula(family = tolower(copula), dim = 2), tau = tau0)
    }

    PMLE.step2 <- nlm(copula_log_lik_in_semi, eta0, S1 = S1, S2 = S2, data1 = dat1, data2 = dat2, copula = copula)
    PMLE.step2 <- PMLE.step2$estimate
  }

  if (copula == "frank") {

    tau0 <- cor(dat1$event_time[dat1$status==1 & dat2$status==1], dat2$event_time[dat1$status==1 & dat2$status==1], method="kendall")
    eta0 <- copula::iTau(copula = copula::archmCopula(family = tolower(copula), dim = 2), tau = tau0)

    PMLE.step2 <- optim(eta0, copula_log_lik_in_semi, S1 = S1, S2 = S2, data1=dat1,data2=dat2, copula=copula,
                        method = "Brent", lower = 0, upper = 30)
    PMLE.step2 <- PMLE.step2$par
  }
  if (copula == "gaussian") {
    eta0 <- sin(cor(dat1$event_time[dat1$status==1 & dat2$status==1], dat2$event_time[dat1$status==1 & dat2$status==1], method="kendall")*pi/2)
    PMLE.step2 <- optim(eta0, copula_log_lik_in_semi, S1 = S1, S2 = S2, data1=dat1,data2=dat2, copula=copula,
                        method = "Brent", lower = -1, upper = 1)
    PMLE.step2 <- PMLE.step2$par
  }
  return(PMLE.step2)
}


copula_param_bootstrap_IR <- function(data, copula = "clayton") {

  tryCatch({

      dat <- data
      dat1 <- dat[,c("event_time.L","status.L")]
      dat2 <- dat[,c("event_time.R","status.R")]
      colnames(dat1) = colnames(dat2) = c("event_time","status")

      KM1 <- survfit(Surv(event_time,status) ~ 1,data = dat1)
      KM2 <- survfit(Surv(event_time,status) ~ 1,data = dat2)
      S1 <- approx(x = KM1$time, y = KM1$surv, xout = dat1$event_time, rule = 2)
      S2 <- approx(x = KM2$time, y = KM2$surv, xout = dat2$event_time, rule = 2)
      S1$y[S1$y==0] <- min(S1$y[S1$y>0])
      S2$y[S2$y==0] <- min(S2$y[S2$y>0])
      S1$y <- (1-1/nrow(dat))*S1$y
      S2$y <- (1-1/nrow(dat))*S2$y

      S1 <- data.frame(x=S1$x, y=S1$y)
      S2 <- data.frame(x=S2$x, y=S2$y)
      S1 <- S1[match(dat1$event_time, S1$x),]
      S2 <- S2[match(dat2$event_time, S2$x),]

      if (copula != "frank" & copula != "copula2" & copula != "gaussian" & copula != "custom") {

        tau0 <- cor(dat1$event_time[dat1$status==1 & dat2$status==1], dat2$event_time[dat1$status==1 & dat2$status==1], method="kendall")
        if (copula == "plackett") {
          eta0 <- iTau(copula = plackettCopula(), tau = tau0)
        } else {
          eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)
        }

        PMLE.step2 <- nlm(copula_log_lik_in_semi, eta0, S1 = S1, S2 = S2, data1 = dat1, data2 = dat2, copula = copula)
        PMLE.step2 <- PMLE.step2$estimate
      }

      if (copula == "frank") {

        tau0 <- cor(dat1$event_time[dat1$status==1 & dat2$status==1], dat2$event_time[dat1$status==1 & dat2$status==1], method="kendall")
        eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

        PMLE.step2 <- optim(eta0, copula_log_lik_in_semi, S1 = S1, S2 = S2, data1=dat1,data2=dat2, copula=copula,
                            method = "Brent", lower = 0, upper = 30)
        PMLE.step2 <- PMLE.step2$par
      }
      if (copula == "gaussian") {
        eta0 <- sin(cor(dat1$event_time[dat1$status==1 & dat2$status==1], dat2$event_time[dat1$status==1 & dat2$status==1], method="kendall")*pi/2)
        PMLE.step2 <- optim(eta0, copula_log_lik_in_semi, S1 = S1, S2 = S2, data1=dat1,data2=dat2, copula=copula,
                            method = "Brent", lower = -1, upper = 1)
        PMLE.step2 <- PMLE.step2$par
      }

      if (copula != "frank" & copula != "custom"){
        hess <- hessian(copula_log_lik_in_semi, x0=PMLE.step2, S1=S1, S2=S2, data1=dat1,data2=dat2, copula=copula,
                        h = 1e-3
        )
      }
      if (copula == "frank") {
        hess <- hessian(copula_log_lik_in_semi, x0=PMLE.step2, S1=S1, S2=S2, data1=dat1,data2=dat2, copula=copula,
                        h = 1e-3)
      }

      if (copula != "copula2" & copula != "custom") {
        score <- sapply(1:nrow(dat), function(x) {grad(copula_log_lik_in_semi, x0=PMLE.step2, S1=S1[x,], S2=S2[x,], data1=dat1[x,], data2=dat2[x,], copula=copula
        ) } )
        score <- matrix(score, nrow = nrow(dat))
        score2 <- sum(sapply(1:nrow(dat), function(x) {matrix(score[x,],nrow=1) %*% matrix(score[x,],ncol=1)})) # sum of grad %*% grad over all subjects
      }

      IR <- sum(diag(solve(hess) %*% score2))

    return(IR)
  }, error = function(err) {return(NA)}
  )
}


# clt_f<-function(u1,u2,eta)
# {
#   c_val<-(1+eta)*(u1*u2)^(-1-eta)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
#   return(c_val)
# }
#
# gh_F<-function(u1,u2,eta)
# {
#   exp(-((-log(u1))^eta + (-log(u2))^eta)^(1/eta))
# }


surv_2_time <- function(p, S){

  S = S[order(S$x),]
  if (p %in% S$y) {t = runif(1, min = min(S$x[S$y==p]), max = max(S$x[S$y==p])) }
  else if (sum(p < S$y) == 0) { t = min(S$x)}
  else if (sum(p < S$y) < nrow(S)) {

    t = min(S$x[S$y <= p])

  }
  else if (sum(p < S$y) == nrow(S)) { t = max(S$x) }

  return(t)

}



get_c <- function(n, R=1, surv) {

  bsample <- function(x, ...) x[sample.int(length(x), replace = TRUE, ...)]

  survival <- surv$surv
  time <- surv$time
  n1 <- length(time)
  if (survival[n1] > 0L) {
    survival <- c(survival, 0)
    time <- c(time, Inf)
  }
  probs <- diff(-c(1, survival))
  matrix(bsample(time, n*R, prob = probs), R, n)
}



gen_copula_RC <- function(data, fit) {


  copula <- fit$copula

  n <- nrow(data)
  param <- fit$param
  S1 <- fit$S1
  S2 <- fit$S2

  if (copula != "copula2") {
    eta <- param
  }

  dat <- data
  dat1 <- dat[,c("event_time.L","status.L")]
  dat2 <- dat[,c("event_time.R","status.R")]
  colnames(dat1) = colnames(dat2) = c("event_time","status")
  dat3 <- dat1
  dat3$event_time <- ifelse(dat3$event_time > dat2$event_time, dat3$event_time, dat2$event_time)
  dat3$status <- 1- dat1$status*dat2$status

  if (copula != "copula2" & copula != "gaussian" & copula != "custom" & copula != "plackett") {
    cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
    Cop <- rCopula(n, cl)

    Cop <- Cop[min(S1$y) < Cop[,1] & Cop[,1] < max(S1$y),]
    Cop <- Cop[min(S2$y) < Cop[,2] & Cop[,2] < max(S2$y),]

    while (nrow(Cop) < n) {
      cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
      Cop.2 <- rCopula(n, cl)
      tmp <- Cop.2[min(S1$y) < Cop.2[,1] & Cop.2[,1] < max(S1$y),]
      tmp <- tmp[min(S2$y) < tmp[,2] & tmp[,2] < max(S2$y),]

      Cop <- rbind(Cop, tmp)
    }

    u<-Cop[1:n,1]
    v<-Cop[1:n,2]
  }



  if (copula == "plackett") {
    cl <- plackettCopula(param = eta)
    Cop <- rCopula(n, cl)

    Cop <- Cop[min(S1$y) < Cop[,1] & Cop[,1] < max(S1$y),]
    Cop <- Cop[min(S2$y) < Cop[,2] & Cop[,2] < max(S2$y),]

    while (nrow(Cop) < n) {
      cl <- plackettCopula(param = eta)
      Cop.2 <- rCopula(n, cl)
      tmp <- Cop.2[min(S1$y) < Cop.2[,1] & Cop.2[,1] < max(S1$y),]
      tmp <- tmp[min(S2$y) < tmp[,2] & tmp[,2] < max(S2$y),]

      Cop <- rbind(Cop, tmp)
    }

    u<-Cop[1:n,1]
    v<-Cop[1:n,2]
  }

  if (copula == "gaussian") {
    cl <- normalCopula(param = eta, dim = 2)
    Cop <- rCopula(n, cl)

    Cop <- Cop[min(S1$y) < Cop[,1] & Cop[,1] < max(S1$y),]
    Cop <- Cop[min(S2$y) < Cop[,2] & Cop[,2] < max(S2$y),]

    while (nrow(Cop) < n) {
      cl <- normalCopula(param = eta, dim = 2)
      Cop.2 <- rCopula(n, cl)
      tmp <- Cop.2[min(S1$y) < Cop.2[,1] & Cop.2[,1] < max(S1$y),]
      tmp <- tmp[min(S2$y) < tmp[,2] & tmp[,2] < max(S2$y),]

      Cop <- rbind(Cop, tmp)
    }

    u<-Cop[1:n,1]
    v<-Cop[1:n,2]
  }


  t1 <- sapply(1:n, function(x) surv_2_time(p=u[x], S=S1))
  t2 <- sapply(1:n, function(x) surv_2_time(p=v[x], S=S2))


  G_cen <- survfit(Surv(event_time, status) ~ 1, data = dat3)
  c <- get_c(n = nrow(dat3), R=1, surv = G_cen)
  c <- as.numeric(c)

  status1<-ifelse(t1<c,1,0)
  status2<-ifelse(t2<c,1,0)
  t1<-ifelse(t1<c,t1,c)
  t2<-ifelse(t2<c,t2,c)

    dat1<-data.frame(id=seq(1:n),enum=seq(1:n),data.frame(time=t1,status=status1))
    dat2<-data.frame(id=seq(1:n),enum=seq(1:n),data.frame(time=t2,status=status2))

    dat<-rbind(dat1,dat2)
    dat<-dat[order(dat$id),]
    dat$enum<-rep(c(1,2),n)

    colnames(dat) <- c("id","ind","event_time","status")
    dat$EYE <- "L"
    dat$EYE[dat$ind==2] <- "R"

    dat <- reshape(dat[,c("id","EYE","event_time","status")], idvar = "id", timevar = "EYE", direction = "wide")

  return(dat)
}




# functions for Rec
# Authors: N. Barthel, C. Geerdens, C. Czado and P. Janssen
NonparaMargEst <- function(data){

  #
  # store data
  #
  gaptimes <- data$gaptimes
  status   <- data$status

  #
  # get total time for each cluster and corresponding censoring status
  #
  time_total   <- rowSums(gaptimes, na.rm = TRUE)
  status_total <- sapply(1:nrow(status), function(x) status[x, length(which(!is.na(status[x,])))])

  #
  # obtain estimate for survival function of total times based on the Nelson-Aalen estimate
  #
  SKM_total    <- summary(survfit(Surv(time_total, status_total) ~ 1, type = "fleming-harrington"))

  #
  # calculate Nelson-Aalen weights pertaining to the survival function of the total time
  #
  weights <- matrix(NA, nrow = length(SKM_total$time))
  weights[1] <- (1 - SKM_total$surv[1])/SKM_total$n.event[1]
  for(j in 2:length(SKM_total$time)){

    weights[j] = (SKM_total$surv[j-1] - SKM_total$surv[j])/SKM_total$n.event[j]

  }

  #
  # estimation of marginal survival functions for the gap times
  #
  init.digits <- 8
  survWeights <- rep(0, nrow(gaptimes))
  for(i in 1:nrow(gaptimes)){

    if(status_total[i] == 1){

      digits <- init.digits
      while(length(which(round(time_total[i], digits = digits) == round(SKM_total$time, digits = digits))) != 1){

        if(length(which(round(time_total[i], digits = digits) == round(SKM_total$time, digits = digits))) > 1){
          digits <- digits + 1
        }else{
          digits <- digits - 1
        }

      }

      survWeights[i] <- weights[round(time_total[i], digits = digits) == round(SKM_total$time, digits = digits)]

    }

  }

  u <- matrix(NA, nrow = nrow(gaptimes), ncol = ncol(gaptimes))
  for(j in 1:ncol(gaptimes)){

    u[1:sum(!is.na(gaptimes[,j])), j] <- sapply(gaptimes[1:sum(!is.na(gaptimes[,j])),j],
                                                function(t) {1 - sum(survWeights*as.numeric(gaptimes[,j] <= t), na.rm = TRUE)})
  }

  return(list(copData = u, status = status))

}


modified_KM <- function(data) {
  data2 <- list(gaptimes = cbind(data$gap1, data$gap2),
                status = cbind(data$status1, data$status2))
  np.est <- NonparaMargEst(data2)
  return(np.est)
}

modified_KM_vine <- function(data) {

  # data2 <- list(gaptimes = cbind(data$gap1, data$gap2, data$gap3, data$gap4),
  #              status = cbind(data$status1, data$status2, data$status3, data$status4))
  np.est <- NonparaMargEst(data)
  return(np.est)
}

MLE_copula_np_PMLE_eta_rec <- function(data, copula) {

  data2 <- list(gaptimes = cbind(data$gap1, data$gap2),
                status = cbind(data$status1, data$status2))
  np.est <- NonparaMargEst(data2)
  S1 <- np.est$copData[,1]
  S2 <- np.est$copData[,2]
  status1 <- np.est$status[,1]
  status2 <- np.est$status[,2]


  if (copula != "frank" & copula != "copula2" & copula != "gaussian" & copula != "custom") {

    tau0 <- cor(data$gap1[data$status1==1 & data$status2==1], data$gap2[data$status1==1 & data$status2==1], method="kendall")
    if (copula == "plackett") {
      eta0 <- iTau(copula = plackettCopula(), tau = tau0)
    } else {
      eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)
    }
    PMLE.step2 <- nlm(copula_log_lik_rec_semi, eta0, S1 = S1, S2 = S2, status1 = status1, status2 = status2, copula = copula)
    PMLE.step2 <- PMLE.step2$estimate
  }

  if (copula == "frank") {
    tau0 <- cor(data$gap1[data$status1==1 & data$status2==1], data$gap2[data$status1==1 & data$status2==1], method="kendall")
    eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

    PMLE.step2 <- optim(eta0, copula_log_lik_rec_semi, S1 = S1, S2 = S2, status1 = status1, status2 = status2, copula=copula,
                        method = "Brent", lower = 0, upper = 30)
    PMLE.step2 <- PMLE.step2$par
  }
  if (copula == "gaussian") {
    eta0 <- sin(cor(data$gap1[data$status1==1 & data$status2==1], data$gap2[data$status1==1 & data$status2==1], method="kendall")*pi/2)
    PMLE.step2 <- optim(eta0, copula_log_lik_rec_semi, S1 = S1, S2 = S2, status1 = status1, status2 = status2, copula=copula,
                        method = "Brent", lower = -1, upper = 1)
    PMLE.step2 <- PMLE.step2$par
  }

  return(PMLE.step2)

}

MLE_copula_np_PMLE_eta_rec_vine <- function(data, fams = c(3, 3, 3)) {

  np.est <- NonparaMargEst(data)

  data_s <- data.frame(gap1 = np.est$copData[,1], gap2 = np.est$copData[,2],
                       gap3 = np.est$copData[,3], gap4 = np.est$copData[,4],
                       status1 = np.est$status[,1], status2 = np.est$status[,2],
                       status3 = np.est$status[,3], status4 = np.est$status[,4])

  dim_cop     <- 4      # maximum cluster size
  fams           <- fams  # e.g. Clayton-Gumbel-Gumbel in tree T1
  paras.init     <- list()
  paras.init$RVM <- list()
  fam.init       <- matrix(c(0, 5, 5, fams[3], # 5 for Frank for T1 and T3
                             0, 0, 5, fams[2],
                             0, 0, 0, fams[1],
                             0, 0, 0,      0), nrow = dim_cop, ncol = dim_cop)
  paras.init$RVM$family <- fam.init
  paras.init$dim_cop <- dim_cop

  PMLE.step2    <- MLE_2stage_seq(data = data_s, paras.init = paras.init)

  return(PMLE.step2)

}


copula_log_lik_rec_semi <- function(para, S1, S2, status1, status2, copula) {

  u1<-S1
  u2<-S2

  if (copula == "plackett") {
    eta <- para
    myCop <- plackettCopula(param=eta)
    c_val <- dCopula(cbind(u1, u2), myCop)
    c_u1_val <- derCOP(u1, u2, cop = PLcop, para = eta)
  }

  if (copula == "gaussian") {
    eta <- para
    myCop <- normalCopula(param=eta, dim=2)
    c_val <- dCopula(cbind(u1, u2), myCop)
    c_u1_val <- cCopula(cbind(u1, u2), myCop)[,2]
  }

  if (copula == "clayton") {
    eta <- para
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
  }

  if (copula == "gumbel") {
    eta <- para
    c_val<- gh_F(u1,u2,eta) * (1/(u1*u2)) * ((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2) * (log(u1)*log(u2))^(eta-1) +
      gh_F(u1,u2,eta) * (-1/(u1*u2)) * (1-eta) * ((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2) * (log(u1)*log(u2))^(eta-1)
    c_u1_val<- ((-log(u1))^(eta) + (-log(u2))^(eta))^(1/eta -1) * (-log(u1))^(eta-1) * (1/u1) * gh_F(u1,u2,eta)
  }

  if (copula == "frank") {
    eta <- para
    c_val <- (-eta)*exp((-eta)*u1)*exp((-eta)*u2)/((exp((-eta))-1)+(exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)) -
      (-eta)*exp((-eta)*u1)*exp((-eta)*u2)*(exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)/((exp((-eta))-1)+(exp((-eta)*u1)-1)*(exp((-eta)*u2)-1))^2
    c_u1_val <- (exp((-eta)*u2)-1)*exp((-eta)*u1)/((1 + (exp((-eta)*u1)-1)*(exp((-eta)*u2)-1)/(exp((-eta))-1))*(exp((-eta))-1))
  }


  term2 <- log(c_u1_val)
  term2 <- ifelse(status1 == 1 & status2 == 0, term2, 0)
  term2[is.infinite(term2)] <- 0

  term4 <- log(c_val)
  term4 <- ifelse(status1 == 1 & status2 == 1, term4, 0)

  term5 <- log(u1)
  term5 <- ifelse(status1 == 0 & status2 == 0, term5, 0)


  logL<-(-1)*sum( term2 + term4 + term5)
  return(logL)
}



copula_param_bootstrap_IR_rec <- function(data, copula = "clayton") {

  tryCatch({

      data2 <- list(gaptimes = cbind(data$gap1, data$gap2),
                    status = cbind(data$status1, data$status2))
      np.est <- NonparaMargEst(data2)
      S1 <- np.est$copData[,1]
      S2 <- np.est$copData[,2]
      status1 <- np.est$status[,1]
      status2 <- np.est$status[,2]

      if (copula != "frank" & copula != "copula2" & copula != "gaussian" & copula != "custom") {

        tau0 <- cor(data$gap1[data$status1==1 & data$status2==1], data$gap2[data$status1==1 & data$status2==1], method="kendall")
        if (copula == "plackett") {
          eta0 <- iTau(copula = plackettCopula(), tau = tau0)
        } else {
          eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)
        }
        PMLE.step2 <- nlm(copula_log_lik_rec_semi, eta0, S1 = S1, S2 = S2, status1 = status1, status2 = status2, copula = copula)
        PMLE.step2 <- PMLE.step2$estimate
      }

      if (copula == "frank") {
        tau0 <- cor(data$gap1[data$status1==1 & data$status2==1], data$gap2[data$status1==1 & data$status2==1], method="kendall")
        eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

        PMLE.step2 <- optim(eta0, copula_log_lik_rec_semi, S1 = S1, S2 = S2, status1 = status1, status2 = status2, copula=copula,
                            method = "Brent", lower = 0, upper = 30)
        PMLE.step2 <- PMLE.step2$par
      }
      if (copula == "gaussian") {
        eta0 <- sin(cor(data$gap1[data$status1==1 & data$status2==1], data$gap2[data$status1==1 & data$status2==1], method="kendall")*pi/2)
        PMLE.step2 <- optim(eta0, copula_log_lik_rec_semi, S1 = S1, S2 = S2, status1 = status1, status2 = status2, copula=copula,
                            method = "Brent", lower = -1, upper = 1)
        PMLE.step2 <- PMLE.step2$par
      }

      if (copula != "frank" & copula != "custom"){
        hess <- hessian(copula_log_lik_rec_semi, x0=PMLE.step2, S1=S1, S2=S2, status1=status1, status2=status2, copula=copula,
                        h = 1e-3
        )
      }
      if (copula == "frank") {
        hess <- hessian(copula_log_lik_rec_semi, x0=PMLE.step2, S1=S1, S2=S2, status1=status1, status2=status2, copula=copula,
                        h = 1e-3)
      }

      if (copula != "copula2" & copula != "custom") {
        score <- sapply(1:nrow(data), function(x) {grad(copula_log_lik_rec_semi, x0=PMLE.step2, S1=S1[x], S2=S2[x], status1=status1[x], status2=status2[x], copula=copula
        ) } )
        score <- matrix(score, nrow = nrow(data))
        score2 <- sum(sapply(1:nrow(data), function(x) {matrix(score[x,],nrow=1) %*% matrix(score[x,],ncol=1)})) # sum of grad %*% grad over all subjects
      }

      IR <- sum(diag(solve(hess) %*% score2))


  }, error = function(err) {return(NA)}
  )

}


copula_param_bootstrap_IR_rec_vine <- function(data, fams = c(3, 3, 3)) {

  tryCatch({

    np.est <- NonparaMargEst(data)

    data_s <- data.frame(gap1 = np.est$copData[,1], gap2 = np.est$copData[,2],
                         gap3 = np.est$copData[,3], gap4 = np.est$copData[,4],
                         status1 = np.est$status[,1], status2 = np.est$status[,2],
                         status3 = np.est$status[,3], status4 = np.est$status[,4])

    dim_cop     <- 4      # maximum cluster size
    fams           <- fams  # e.g. Clayton-Gumbel-Gumbel (3-4-4) in tree T1 or Frank-Gumbel-Gumbel (5-4-4) in tree T1
    paras.init     <- list()
    paras.init$RVM <- list()
    fam.init       <- matrix(c(0, 5, 5, fams[3], # 5 for Frank for T1 and T3
                               0, 0, 5, fams[2],
                               0, 0, 0, fams[1],
                               0, 0, 0,      0), nrow = dim_cop, ncol = dim_cop)
    paras.init$RVM$family <- fam.init
    paras.init$dim_cop <- dim_cop

    PMLE.step2    <- MLE_2stage_seq(data = data_s, paras.init = paras.init)



    hess <- hessian(loglik2stage, x0=PMLE.step2,
                    data = data_s,
                    fam = c(paras.init$RVM$family[4,3:1],
                            paras.init$RVM$family[3,2:1],
                            paras.init$RVM$family[2,1]))
    hess = hess * (-1) # llk to neg-llk

    score <- sapply(1:nrow(data_s), function(x) {grad(loglik2stage, x0=PMLE.step2, data=data_s[x,],
                                                      fam = c(paras.init$RVM$family[4,3:1],
                                                              paras.init$RVM$family[3,2:1],
                                                              paras.init$RVM$family[2,1])) } )
    score <- t(score)
    tmp <- lapply(1:nrow(data_s), function(x) {matrix(score[x,],nrow=6) %*% matrix(score[x,],ncol=6)}) # Note: dim =2 for copula2; sum of grad %*% grad over all subjects
    score2 <- matrix(c(sum(sapply(1:nrow(data_s), function(x) tmp[[x]][1,1])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][2,1])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][3,1])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][4,1])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][5,1])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][6,1])),

                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][1,2])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][2,2])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][3,2])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][4,2])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][5,2])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][6,2])),

                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][1,3])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][2,3])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][3,3])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][4,3])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][5,3])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][6,3])),

                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][1,4])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][2,4])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][3,4])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][4,4])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][5,4])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][6,4])),

                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][1,5])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][2,5])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][3,5])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][4,5])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][5,5])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][6,5])),

                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][1,6])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][2,6])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][3,6])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][4,6])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][5,6])),
                       sum(sapply(1:nrow(data_s), function(x) tmp[[x]][6,6]))
    ),
    ncol = 6, nrow = 6)

    IR <- sum(diag(pseudoinverse(hess) %*% score2))
    return(IR)

  }, error = function(err) {return(NA)}
  )

}


gen_copula_RC_rec <- function(data, fit) {


  copula <- fit$copula

    n <- nrow(data)
    param <- fit$param
    S1 <- data.frame(x=data$gap1, y=fit$S1)
    S2 <- data.frame(x=data$gap2, y=fit$S2)

    if (copula != "copula2") {
      eta <- param
    }

    data$total <- data$gap1 + data$gap2
    data$status <- 1-data$status2

    if (copula != "copula2" & copula != "gaussian" & copula != "custom" & copula != "plackett") {
      cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
      Cop <- rCopula(n, cl)

      Cop <- Cop[min(S1$y) < Cop[,1] & Cop[,1] < max(S1$y),]
      Cop <- Cop[min(S2$y) < Cop[,2] & Cop[,2] < max(S2$y),]

      while (nrow(Cop) < n) {
        cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
        Cop.2 <- rCopula(n, cl)
        tmp <- Cop.2[min(S1$y) < Cop.2[,1] & Cop.2[,1] < max(S1$y),]
        tmp <- tmp[min(S2$y) < tmp[,2] & tmp[,2] < max(S2$y),]

        Cop <- rbind(Cop, tmp)
      }

      u<-Cop[1:n,1]
      v<-Cop[1:n,2]
    }

    if (copula == "plackett") {
      cl <- plackettCopula(param = eta)
      Cop <- rCopula(n, cl)

      Cop <- Cop[min(S1$y) < Cop[,1] & Cop[,1] < max(S1$y),]
      Cop <- Cop[min(S2$y) < Cop[,2] & Cop[,2] < max(S2$y),]

      while (nrow(Cop) < n) {
        cl <- plackettCopula(param = eta)
        Cop.2 <- rCopula(n, cl)
        tmp <- Cop.2[min(S1$y) < Cop.2[,1] & Cop.2[,1] < max(S1$y),]
        tmp <- tmp[min(S2$y) < tmp[,2] & tmp[,2] < max(S2$y),]

        Cop <- rbind(Cop, tmp)
      }

      u<-Cop[1:n,1]
      v<-Cop[1:n,2]
    }

    if (copula == "gaussian") {
      cl <- normalCopula(param = eta, dim = 2)
      Cop <- rCopula(n, cl)

      Cop <- Cop[min(S1$y) < Cop[,1] & Cop[,1] < max(S1$y),]
      Cop <- Cop[min(S2$y) < Cop[,2] & Cop[,2] < max(S2$y),]

      while (nrow(Cop) < n) {
        cl <- normalCopula(param = eta, dim = 2)
        Cop.2 <- rCopula(n, cl)
        tmp <- Cop.2[min(S1$y) < Cop.2[,1] & Cop.2[,1] < max(S1$y),]
        tmp <- tmp[min(S2$y) < tmp[,2] & tmp[,2] < max(S2$y),]

        Cop <- rbind(Cop, tmp)
      }

      u<-Cop[1:n,1]
      v<-Cop[1:n,2]
    }

    g1 <- sapply(1:n, function(x) surv_2_time(p=u[x], S=S1))
    g2 <- sapply(1:n, function(x) surv_2_time(p=v[x], S=S2))

    t1 <- g1
    t2 <- g1+g2

    G_cen <- survfit(Surv(total, status) ~ 1, data = data)
    c <- get_c(n = nrow(data), R=1, surv = G_cen)
    c <- as.numeric(c)
    d <- ifelse(t1 > c, 1, 2)

    y1 <- ifelse(g1<=c, g1, c)
    status1 <- ifelse(g1<=c, 1, 0)
    c2 <- (c-g1)*(g1<=c)
    y2 <- ifelse(g2<=c2, g2, c2)
    status2 <- ifelse(g2<=c2, 1, 0)

    dat <- data.frame(id = seq(1:length(y1)), gap1 = y1, status1 = status1, gap2 = y2, status2 = status2, d = d)

  return(dat)
}


gen_copula_RC_rec_vine <- function(data, fit) {


  fams <- fit$fams

  # KM margins
  n <- nrow(data$gaptimes)
  param <- fit$param
  S1 <- data.frame(x=data$gaptimes[,1], y=fit$S1)
  S2 <- data.frame(x=data$gaptimes[,2], y=fit$S2)
  S3 <- data.frame(x=data$gaptimes[,3], y=fit$S3)
  S4 <- data.frame(x=data$gaptimes[,4], y=fit$S4)

  S1 <- na.omit(S1)
  S2 <- na.omit(S2)
  S3 <- na.omit(S3)
  S4 <- na.omit(S4)

  # setup Vine data simulation
  dim_cop     <- 4      # maximum cluster size
  fam.init <- matrix(c(0, 5, 5, fams[3], # 5 for Frank for T1 and T3
                       0, 0, 5, fams[2],
                       0, 0, 0, fams[1],
                       0, 0, 0,      0), nrow = dim_cop, ncol = dim_cop)
  theta.init <- matrix(c(0, param[6], param[5], param[3], # 5 for Frank for T1 and T3
                         0,        0, param[4], param[2],
                         0,        0,        0, param[1],
                         0,        0,        0,       0), nrow = dim_cop, ncol = dim_cop)
  Matrix     <- matrix(c(4, 1, 2, 3,
                         0, 3, 1, 2,
                         0, 0, 2, 1,
                         0, 0, 0, 1), nrow = dim_cop, ncol = dim_cop) # Lower (or upper) triangular d x d matrix that defines the R-vine tree structure.
  par2 <-  matrix(0, dim_cop, dim_cop) # some copula needs second eta
  RVM  <- RVineMatrix(Matrix = Matrix, family = fam.init, par = theta.init, par2 = par2)

  # generate non-censored 'dim_cop'-dimensional vine copula data with 'ncluster' individuals based
  # on D-vine copula specification given in 'RVM'
  # set.seed(2468+1)
  U <- RVineSim(N = n, RVM = RVM) # nX4, gap1 to gap4

  # input original data for KM of censoring variable
  data$total_cen <- rowSums(data$gaptimes, na.rm = T)
  data$status_cen <- rowSums(1-data$status, na.rm = T)

  # from U back to g1 g2 g3 g4
  g1 <- sapply(1:n, function(x) surv_2_time(p=U[x,1], S=S1))
  g2 <- sapply(1:n, function(x) surv_2_time(p=U[x,2], S=S2))
  g3 <- sapply(1:n, function(x) surv_2_time(p=U[x,3], S=S3))
  g4 <- sapply(1:n, function(x) surv_2_time(p=U[x,4], S=S4))

  # gap times to calendar times
  t1 <- g1
  t2 <- g1+g2
  t3 <- g1+g2+g3
  t4 <- g1+g2+g3+g4

  event = cbind(t1, t2, t3, t4)


    G_cen <- survfit(Surv(data$total_cen, data$status_cen) ~ 1)
    cens <- get_c(n = nrow(data$gaptimes), R=1, surv = G_cen)
    cens <- as.numeric(cens)

  # - generate observed right-censored event time (T) data via comparison of event times and
  #   censoring times
  # - generate matrix with right-censoring status
  # - determine individual cluster size
  obsT_temp    <- matrix(NA, nrow = n, ncol = dim_cop)
  status_temp  <- matrix(NA, nrow = n, ncol = dim_cop)
  scluster     <- rep(NA, n)
  for(i in 1:n){

    if(length(which(event[i,] > cens[i])) == 0){

      scluster[i]      <- dim_cop
      obsT_temp[i,]    <- event[i, 1:dim_cop]
      status_temp[i,]  <- rep(1, dim_cop)

    }else{

      scluster[i] <- min(which(event[i,] > cens[i]))
      if(scluster[i] == 1){

        obsT_temp[i, 1]    <- cens[i]
        status_temp[i, 1]  <- 0

      }else{

        obsT_temp[i, 1:scluster[i]]    <- c(event[i, 1:(scluster[i]-1)], cens[i])
        status_temp[i, 1:scluster[i]]  <- c(rep(1, (scluster[i]-1)), 0)

      }

    }

  }


  obsT    <- obsT_temp[which(rowSums(status_temp, na.rm = TRUE) == (dim_cop - 1) &
                               rowSums(!is.na(status_temp)) == dim_cop),]
  status  <- status_temp[which(rowSums(status_temp, na.rm = TRUE) == (dim_cop - 1) &
                                 rowSums(!is.na(status_temp)) == dim_cop),]
  obsT    <- rbind(obsT, obsT_temp[which(rowSums(status_temp, na.rm = TRUE) == dim_cop &
                                           rowSums(!is.na(status_temp)) == dim_cop),])
  status  <- rbind(status, status_temp[which(rowSums(status_temp, na.rm = TRUE) == dim_cop &
                                               rowSums(!is.na(status_temp)) == dim_cop),])
  for(n.events in (max(scluster)-2):0){

    obsT   <- rbind(obsT, obsT_temp[which(rowSums(status_temp, na.rm = TRUE) == n.events &
                                            rowSums(!is.na(status_temp)) == (n.events + 1)),])
    status <- rbind(status, status_temp[which(rowSums(status_temp, na.rm = TRUE) == n.events &
                                                rowSums(!is.na(status_temp)) == (n.events + 1)),])

  }


  obsG <- sapply(2:max(scluster), function(x) obsT[, x] - obsT[, x-1])
  obsG <- cbind(obsT[,1], obsG)

  return(list(gaptimes = obsG, status = status))

}


Marginal_Turnbull <- function(data) {

  dat <- data	#
  dat1 <- dat[,c("Left.L","Right.L","status.L")] #
  dat2 <- dat[,c("Left.R","Right.R","status.R")] #
  colnames(dat1) = colnames(dat2) = c("Left","Right","status")

  TB1 <- ic_np(dat1[,c("Left","Right")])
  TB2 <- ic_np(dat2[,c("Left","Right")])

  u1_left <- approx(x = c(0,TB1$T_bull_Intervals[1,]), y = c(1,1,(1-cumsum(TB1$p_hat))[-length(TB1$p_hat)]), xout = dat1$Left, rule = 2)
  u1_right <- approx(x = c(0,TB1$T_bull_Intervals[2,]), y = c(1,(1-cumsum(TB1$p_hat))), xout = dat1$Right, rule = 2)
  u2_left <- approx(x = c(0,TB2$T_bull_Intervals[1,]), y = c(1,1,(1-cumsum(TB2$p_hat))[-length(TB2$p_hat)]), xout = dat2$Left, rule = 2)
  u2_right <- approx(x = c(0,TB2$T_bull_Intervals[2,]), y = c(1,(1-cumsum(TB2$p_hat))), xout = dat2$Right, rule = 2)

  return(list(u1_left=u1_left, u1_right=u1_right, u2_left=u2_left, u2_right=u2_right))

}


MLE_ic_copula_Turnbull_eta <- function(data, copula) {

  dat <- data
  dat1 <- dat[,c("Left.L","Right.L","status.L")]
  dat2 <- dat[,c("Left.R","Right.R","status.R")]
  colnames(dat1) = colnames(dat2) = c("Left","Right","status")

  TB1 <- ic_np(dat1[,c("Left","Right")])
  TB2 <- ic_np(dat2[,c("Left","Right")])

  u1_left <- approx(x = c(0,TB1$T_bull_Intervals[1,]), y = c(1,1,(1-cumsum(TB1$p_hat))[-length(TB1$p_hat)]), xout = dat1$Left, rule = 2)
  u1_right <- approx(x = c(0,TB1$T_bull_Intervals[2,]), y = c(1,(1-cumsum(TB1$p_hat))), xout = dat1$Right, rule = 2)
  u2_left <- approx(x = c(0,TB2$T_bull_Intervals[1,]), y = c(1,1,(1-cumsum(TB2$p_hat))[-length(TB2$p_hat)]), xout = dat2$Left, rule = 2)
  u2_right <- approx(x = c(0,TB2$T_bull_Intervals[2,]), y = c(1,(1-cumsum(TB2$p_hat))), xout = dat2$Right, rule = 2)

  if (copula == "frank") {

    tau0 <- cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")
    eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

    PMLE.step2 <- optim(eta0, ic_copula_log_lik_in, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right, data1 = dat1, data2 = dat2, copula = copula, hessian = F,
                        method = "Brent", lower = 0, upper = 30)
    PMLE.step2 <- PMLE.step2$par

  }

  if (copula == "gaussian") {

    eta0 <- sin(cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")*pi/2)
    PMLE.step2 <- optim(eta0, ic_copula_log_lik_in, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right, data1 = dat1, data2 = dat2, copula = copula, hessian = F,
                        method = "Brent", lower = -1, upper = 1)
    PMLE.step2 <- PMLE.step2$par

  }

  if (copula == "plackett") {

    tau0 <- cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")
    eta0 <- iTau(copula = plackettCopula(), tau = tau0)

    PMLE.step2 <- optim(eta0, ic_copula_log_lik_in, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right, data1 = dat1, data2 = dat2, copula = copula, hessian = F,
                        method = "Brent", lower = 0.001, upper = 100000)
    PMLE.step2 <- PMLE.step2$par
  }

  if (copula != "copula2" & copula != "frank" & copula != "gaussian" & copula != "custom" & copula != "plackett") {

    tau0 <- cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")
    eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

    PMLE.step2 <- nlm(ic_copula_log_lik_in, eta0, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right,
                      data1 = dat1, data2 = dat2, copula = copula, hessian = F)
    PMLE.step2 <- PMLE.step2$estimate

  }

  return(PMLE.step2)

}



ic_copula_param_bootstrap_IR <- function(data, copula) {

  tryCatch({


      dat <- data
      dat1 <- dat[,c("Left.L","Right.L","status.L")]
      dat2 <- dat[,c("Left.R","Right.R","status.R")]
      colnames(dat1) = colnames(dat2) = c("Left","Right","status")

      TB1 <- ic_np(dat1[,c("Left","Right")])
      TB2 <- ic_np(dat2[,c("Left","Right")])

      u1_left <- approx(x = c(0,TB1$T_bull_Intervals[1,]), y = c(1,1,(1-cumsum(TB1$p_hat))[-length(TB1$p_hat)]), xout = dat1$Left, rule = 2)
      u1_right <- approx(x = c(0,TB1$T_bull_Intervals[2,]), y = c(1,(1-cumsum(TB1$p_hat))), xout = dat1$Right, rule = 2)
      u2_left <- approx(x = c(0,TB2$T_bull_Intervals[1,]), y = c(1,1,(1-cumsum(TB2$p_hat))[-length(TB2$p_hat)]), xout = dat2$Left, rule = 2)
      u2_right <- approx(x = c(0,TB2$T_bull_Intervals[2,]), y = c(1,(1-cumsum(TB2$p_hat))), xout = dat2$Right, rule = 2)

      u1_left <- data.frame(x = u1_left$x, y = u1_left$y)
      u1_right <- data.frame(x = u1_right$x, y = u1_right$y)
      u2_left <- data.frame(x = u2_left$x, y = u2_left$y)
      u2_right <- data.frame(x = u2_right$x, y = u2_right$y)


    if (copula == "frank") {

      tau0 <- cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")
      eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

      PMLE.step2 <- optim(eta0, ic_copula_log_lik_in, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right, data1 = dat1, data2 = dat2, copula = copula, hessian = F,
                          method = "Brent", lower = 0, upper = 30)
      PMLE.step2 <- PMLE.step2$par

    }

    if (copula == "gaussian") {

      eta0 <- sin(cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")*pi/2)
      PMLE.step2 <- optim(eta0, ic_copula_log_lik_in, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right, data1 = dat1, data2 = dat2, copula = copula, hessian = F,
                          method = "Brent", lower = -1, upper = 1)
      PMLE.step2 <- PMLE.step2$par

    }


    if (copula == "plackett") {

      tau0 <- cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")
      eta0 <- iTau(copula = plackettCopula(), tau = tau0)

      PMLE.step2 <- optim(eta0, ic_copula_log_lik_in, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right, data1 = dat1, data2 = dat2, copula = copula, hessian = F,
                          method = "Brent", lower = 0, upper = 100000)
      PMLE.step2 <- PMLE.step2$par
    }

    if (copula != "copula2" & copula != "frank" & copula != "gaussian" & copula != "custom" & copula != "plackett") {

      tau0 <- cor(0.5*(dat1$Left+dat1$Right)[dat1$status==1 & dat2$status==1], 0.5*(dat2$Left+dat2$Right)[dat1$status==1 & dat2$status==1], method="kendall")
      eta0 <- iTau(copula = archmCopula(family = tolower(copula), dim = 2), tau = tau0)

      PMLE.step2 <- nlm(ic_copula_log_lik_in, eta0, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right,
                        data1 = dat1, data2 = dat2, copula = copula, hessian = F)
      PMLE.step2 <- PMLE.step2$estimate

    }


    if (copula != "frank" & copula != "copula2" & copula != "custom") {
      hess <- hessian(ic_copula_log_lik_in, x0=PMLE.step2, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right,
                                data1 = dat1,data2 = dat2, copula = copula)
    }

    if (copula == "frank") {
      hess <- hessian(ic_copula_log_lik_in, x0=PMLE.step2, u1_left = u1_left, u1_right = u1_right, u2_left = u2_left, u2_right = u2_right,
                                data1 = dat1,data2 = dat2, copula = copula)
    }


    if (copula == "frank") {

      score <- sapply(1:nrow(dat), function(x) {grad(ic_copula_log_lik_in, x0=PMLE.step2, u1_left = u1_left[x,], u1_right = u1_right[x,], u2_left = u2_left[x,], u2_right = u2_right[x,], data1=dat1[x,], data2=dat2[x,], copula = copula) } )
      score <- matrix(score, nrow = nrow(dat))
      score2 <- sum(sapply(1:nrow(dat), function(x) {matrix(score[x,],nrow=1) %*% matrix(score[x,],ncol=1)})) # sum of grad %*% grad over all subjects

    }


    if (copula != "copula2" & copula != "frank" & copula != "custom") {
      score <- sapply(1:nrow(dat), function(x) {grad(ic_copula_log_lik_in, x0=PMLE.step2, u1_left = u1_left[x,], u1_right = u1_right[x,], u2_left = u2_left[x,], u2_right = u2_right[x,], data1=dat1[x,], data2=dat2[x,], copula = copula) } )
      score <- matrix(score, nrow = nrow(dat))
      score2 <- sum(sapply(1:nrow(dat), function(x) {matrix(score[x,],nrow=1) %*% matrix(score[x,],ncol=1)})) # sum of grad %*% grad over all subjects

    }

    IR <- sum(diag(solve(hess) %*% score2))

    return(IR)

  }, error = function(err) {return(NA)}
  )

}



ic_copula_log_lik_in <- function(para, u1_left, u1_right, u2_left, u2_right, data1, data2, copula) {

  u1_left<-u1_left$y
  u1_right<-u1_right$y
  u2_left<-u2_left$y
  u2_right<-u2_right$y

  if (copula == "plackett"){

    eta <- para
    myCop <- plackettCopula(param=eta)

    C_val_1 <- pCopula(cbind(u1_left, u2_left), myCop)
    C_val_2 <- pCopula(cbind(u1_left, u2_right), myCop)
    C_val_3 <- pCopula(cbind(u1_right, u2_left), myCop)
    C_val_4 <- pCopula(cbind(u1_right, u2_right), myCop)

  }

  if (copula == "gaussian"){

    eta <- para
    myCop <- normalCopula(param=eta, dim=2)

    C_val_1 <- pCopula(cbind(u1_left, u2_left), myCop)
    C_val_2 <- pCopula(cbind(u1_left, u2_right), myCop)
    C_val_3 <- pCopula(cbind(u1_right, u2_left), myCop)
    C_val_4 <- pCopula(cbind(u1_right, u2_right), myCop)

  }

  if (copula == "clayton"){

    eta <- para

    C_val_1 <-(u1_left^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_left^(-eta)+u2_right^(-eta)-1)^(-1/eta)
    C_val_3 <-(u1_right^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_4 <-(u1_right^(-eta)+u2_right^(-eta)-1)^(-1/eta)

  }

  if (copula == "gumbel"){

    eta <- para

    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_left))^eta + (-log(u2_right))^eta)^(1/eta))
    C_val_3 <- exp(-((-log(u1_right))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_4 <- exp(-((-log(u1_right))^eta + (-log(u2_right))^eta)^(1/eta))

  }

  if (copula == "frank"){

    eta <- para

    C_val_1 <- 1/(-eta) * log(1 + (exp(-eta*u1_left)-1)*(exp(-eta*u2_left)-1)/(exp(-eta)-1))
    C_val_2 <- 1/(-eta) * log(1 + (exp(-eta*u1_left)-1)*(exp(-eta*u2_right)-1)/(exp(-eta)-1))
    C_val_3 <- 1/(-eta) * log(1 + (exp(-eta*u1_right)-1)*(exp(-eta*u2_left)-1)/(exp(-eta)-1))
    C_val_4 <- 1/(-eta) * log(1 + (exp(-eta*u1_right)-1)*(exp(-eta*u2_right)-1)/(exp(-eta)-1))

  }

  if (copula == "joe"){

    eta <- para

    C_val_1 <- 1 - ((1-u1_left)^eta + (1-u2_left)^eta - ((1-u1_left)^eta)*((1-u2_left)^eta) )^(1/eta)
    C_val_2 <- 1 - ((1-u1_left)^eta + (1-u2_right)^eta - ((1-u1_left)^eta)*((1-u2_right)^eta) )^(1/eta)
    C_val_3 <- 1 - ((1-u1_right)^eta + (1-u2_left)^eta - ((1-u1_right)^eta)*((1-u2_left)^eta) )^(1/eta)
    C_val_4 <- 1 - ((1-u1_right)^eta + (1-u2_right)^eta - ((1-u1_right)^eta)*((1-u2_right)^eta) )^(1/eta)

  }

  term1 <- log(C_val_1 - C_val_2 - C_val_3 + C_val_4)
  term1 <- ifelse((data1[,"status"] == 1) & (data2[,"status"] == 1), term1, 0)

  term2 <- log(C_val_1 - C_val_3)
  term2 <- ifelse((data1[,"status"] == 1) & (data2[,"status"] == 0), term2, 0)

  term3 <- log(C_val_1 - C_val_2)
  term3 <- ifelse((data1[,"status"] == 0) & (data2[,"status"] == 1), term3, 0)

  term4 <- log(C_val_1)
  term4 <- ifelse((data1[,"status"] == 0) & (data2[,"status"] == 0), term4, 0)

  logL<-(-1)*sum( term1 + term2 + term3 + term4 )
  return(logL)

}


gen_copula_IC <- function(data, fit) {

  n <- nrow(data)
  param <- fit$param
  copula = fit$copula

    u1_left <- fit$u1_left
    u1_right <- fit$u1_right
    u2_left <- fit$u2_left
    u2_right <- fit$u2_right

    S1 <- data.frame(x = c(u1_left$x, u1_right$x), y = c(u1_left$y, u1_right$y))
    S2 <- data.frame(x = c(u2_left$x, u2_right$x), y = c(u2_left$y, u2_right$y))

    if (copula != "copula2") {
      eta <- param
    }

    if (copula != "copula2" & copula != "gaussian" & copula != "custom" & copula != "plackett") {
      cl <- archmCopula(family = tolower(copula), param = eta, dim = 2)
      Cop <- rCopula(n, cl)
      u<-Cop[,1]
      v<-Cop[,2]
    }

    if (copula == "gaussian") {
      cl <- normalCopula(param = eta, dim = 2)
      Cop <- rCopula(n, cl)
      u<-Cop[,1]
      v<-Cop[,2]
    }

    if (copula == "plackett") {
      cl <- plackettCopula(param = eta)
      Cop <- rCopula(n, cl)
      u<-Cop[,1]
      v<-Cop[,2]
    }


    # if (copula == "custom") {
    #   Cop <- custom_gen(eta = eta, n = n)
    #   u<-Cop[,1]
    #   v<-Cop[,2]
    # }

    dat <- data
    dat1 <- dat[,c("Left.L","Right.L","status.L")]
    dat2 <- dat[,c("Left.R","Right.R","status.R")]
    colnames(dat1) = colnames(dat2) = c("Left","Right","status")

    dat1<-data.frame(id=seq(1:n),enum=rep(1,n))
    dat2<-data.frame(id=seq(1:n),enum=rep(2,n))

    t1 <- sapply(1:n, function(x) surv_2_time(p=u[x], S=S1))
    t2 <- sapply(1:n, function(x) surv_2_time(p=v[x], S=S2))


  seeds <- sample(1:10^3,1)
  mu_estimate <- 1/mean(mean(dat$Right.L[dat$status.L==1] - dat$Left.L[dat$status.L==1]),
                      mean(dat$Right.R[dat$status.R==1] - dat$Left.R[dat$status.R==1]))
  K <- optim(10, cen_exponential_IC, n=n, mu=mu_estimate, p2 = 1, t1=t1, t2=t2, rate=1-mean(c(dat$status.L,dat$status.R)),
             method = "Brent", lower = 0, upper = 50, seeds=seeds)
  K <- round(K$par)


  dat1<-cbind(dat1,data.frame(time=t1))
  dat2<-cbind(dat2,data.frame(time=t2))

  dat<-rbind(dat1,dat2)
  dat<-dat[order(dat$id),]

  dat = event_to_Copula_interval(dat,K=K,mu=mu_estimate,p=1)

    colnames(dat) <- c("id","ind","event_time","Left","Right","status","event")
    dat$EYE <- "L"
    dat$EYE[dat$ind==2] <- "R"
    dat <- reshape(dat[,c("id","EYE","Left","Right","status")], idvar = "id", timevar = "EYE", direction = "wide")

  return(dat)
}


event_to_Copula_interval <- function(data,K,mu,p) {
  data <- data.frame(data)
  enum = 1
  dat1 <- subset(data,enum==1)
  dat2 <- subset(data,enum==2)

  for (i in 1:nrow(dat1)) {

    a = cumsum(rexp(K,rate=1/mu))
    u = runif(K)
    ind = ifelse(u<=p,1,0)
    A = a[ind==1]

    if(dat1$time[i]>max(A)) {
      dat1$Left[i]  <- max(A)
      dat1$Right[i] <- Inf
      dat1$status[i] = 0
      dat1$event[i] = 0
    } else if (dat1$time[i]<=min(A)) {
      dat1$Left[i] <- 0
      dat1$Right[i] <- min(A)
      dat1$status[i] = 1
      dat1$event[i] = 2
    } else {
      dat1$Left[i]  <- max(A[dat1$time[i]>A])
      dat1$Right[i] <- min(A[dat1$time[i]<=A])
      dat1$status[i] = 1
      dat1$event[i] = 3
    }

    if(dat2$time[i]>max(A)) {
      dat2$Left[i]  <- max(A)
      dat2$Right[i] <- Inf
      dat2$status[i] = 0
      dat2$event[i] = 0
    } else if (dat2$time[i]<=min(A)) {
      dat2$Left[i] <- 0
      dat2$Right[i] <- min(A)
      dat2$status[i] = 1
      dat2$event[i] = 2
    } else {
      dat2$Left[i]  <- max(A[dat2$time[i]>A])
      dat2$Right[i] <- min(A[dat2$time[i]<=A])
      dat2$status[i] = 1
      dat2$event[i] = 3
    }
  }

  dat = rbind(dat1,dat2)
  dat = dat[order(dat$id),]
  return(dat)
}

cen_exponential_IC <- function(p, n, mu=0.6, p2 = 1, t1, t2, rate, seeds){

  set.seed(seeds)

  status1 <- seq(1:n)
  status2 <- seq(1:n)

  for (i in 1:n) {

    a = cumsum(rexp(p,rate=1/mu))
    u = runif(p)
    ind = ifelse(u<=p2,1,0)
    A = a[ind==1]

    if(t1[i]>max(A)) {
      status1[i] = 0
    } else if (t1[i]<=min(A)) {
      status1[i] = 1
    } else {
      status1[i] = 1
    }

    if(t2[i]>max(A)) {
      status2[i] = 0
    } else if (t2[i]<=min(A)) {
      status2[i] = 1
    } else {
      status2[i] = 1
    }

  }

  status <- 1- mean(c(status1,status2))
  l <- (status - rate)^2
  return(l)

}



# Vine
## function for setting lower bounds for parameters of pair-copulas in optimization procedure
lower.bounds <- function(fam){
  switch(fam,
         '1' = -1,
         '2' = -1,
         '3' = 0.01,
         '4' = 1.01,
         '5' = -Inf,
         '6' = 0,
         '13' = 0,
         '14' = 1.01,
         '16' = 0,
         '23' = -Inf,
         '24' = -Inf,
         '26' = -Inf,
         '33' = -Inf,
         '34' = -Inf,
         '36' = -Inf)
}
lower.bounds <- Vectorize(lower.bounds)

## function for setting upper bounds for parameters of pair-copulas in optimization procedure
upper.bounds <- function(fam){
  switch(fam,
         '1' = 1,
         '2' = 1,
         '3' = Inf,
         '4' = Inf,
         '5' = Inf,
         '6' = Inf,
         '13' = Inf,
         '14' = Inf,
         '16' = Inf,
         '23' = 0,
         '24' = -1,
         '26' = 0,
         '33' = 0,
         '34' = -1,
         '36' = 0)
}
upper.bounds <- Vectorize(upper.bounds)

## function for setting starting values for parameters of pair-copulas in optimization procedure
init.par <- function(fam){
  switch(fam,
         '1' = 0,
         '2' = 0,
         '3' = 2,
         '4' = 2,
         '5' = 1,
         '6' = 1,
         '13' = 2,
         '14' = 2,
         '16' = 1,
         '23' = -2,
         '24' = -2,
         '26' = -1,
         '33' = -2,
         '34' = -2,
         '36' = -1)
}
init.par <- Vectorize(init.par)


integrated.density <- function(data, DVM){

  temp <- RVineLogLik(data = data, RVM = DVM, separate = TRUE)

  return(BiCopHfunc2(u1 = temp$V$direct[2, 1, ],
                     u2 = temp$V$indirect[2, 2, ],
                     family = DVM$family[2, 1],
                     par = DVM$par[2, 1],
                     par2 = DVM$par2[2, 1])
         * apply(X = temp$V$value,
                 FUN = function(M) exp(sum(M[-1, -1])),
                 MARGIN = 3)
  )

}


# llk for sequential estimation
loglik2stage_dim2 <- function(par, fam, gap1, gap2, status2){

  par2 <- par
  fam2 <- fam

  ### cluster size 2 ###
  data2 <- cbind(gap1, gap2)

  ll2 <- 0
  for (i in 1:length(gap2)){

    ll2.temp <- ifelse(status2[i] == 1, log(BiCopPDF(u1 = data2[i, 1], u2 = data2[i, 2], family = fam2, par = par2)),
                       log(BiCopHfunc1(u1 = data2[i, 1], u2 = data2[i, 2], family = fam2, par = par2)))

    if(is.finite(ll2.temp)){
      ll2 <- ll2 + ll2.temp
    }

  }

  return(ll2)

}



# llk in full loglik2stage
loglik2stage_dim3 <- function(par, fam, gap1, gap2, gap3, status3){

  par3 <- par

  RVM3 <- RVineMatrix(Matrix = matrix(c(3, 1, 2,
                                        0, 2, 1,
                                        0, 0, 1), 3, 3),
                      family = matrix(c(0, fam[3], fam[2],
                                        0, 0, fam[1],
                                        0, 0, 0), 3, 3),
                      par = matrix(c(0, par3[3], par3[2],
                                     0, 0, par3[1],
                                     0, 0, 0), 3, 3))

  ### cluster size 3 ###
  data3 <- cbind(gap1, gap2, gap3)

  ll3 <- 0
  for (i in 1:length(gap3)){

    ll3.temp <- ifelse(status3[i] == 1, log(RVinePDF(data3[i,], RVM = RVM3)),
                       log(integrated.density(data3[i, ], DVM = RVM3)))

    if(is.finite(ll3.temp)){
      ll3 <- ll3 + ll3.temp
    }

  }

  return(ll3)

}



# llk in full loglik2stage
loglik2stage_dim4 <- function(par, fam, gap1, gap2, gap3, gap4, status4){

  par4 <- par

  RVM4 <- RVineMatrix(Matrix = matrix(c(4, 1, 2, 3,
                                        0, 3, 1, 2,
                                        0, 0, 2, 1,
                                        0, 0, 0, 1), 4, 4),
                      family = matrix(c(0, fam[6], fam[5], fam[3],
                                        0,      0, fam[4], fam[2],
                                        0,      0,      0, fam[1],
                                        0,      0,      0,      0), 4, 4),
                      par = matrix(c(0, par4[6], par4[5], par4[3],
                                     0,       0, par4[4], par4[2],
                                     0,       0,       0, par4[1],
                                     0,       0,       0,      0), 4, 4))


  ### cluster size 4 ###
  ll4 <- 0

  data4 <- cbind(gap1, gap2, gap3, gap4)
  for (i in 1:length(gap4)){

    ll4.temp <- ifelse(status4[i] == 1, log(RVinePDF(data4[i,], RVM = RVM4)),
                       log(integrated.density(data4[i,], DVM = RVM4)))

    if(is.finite(ll4.temp)){
      ll4 <- ll4 + ll4.temp
    }

  }

  return(ll4)

}



# full llk for global estimation
loglik2stage <- function(par, data, fam){

  ll       <- 0
  par.temp <- par

  #
  # sum of loglikelihood contributions of clusters of size 4
  #
  gap1   <- data$gap1[1:length(data$gap4)]
  gap2   <- data$gap2[1:length(data$gap4)]
  gap3   <- data$gap3[1:length(data$gap4)]
  gap4   <- data$gap4
  status <- data$status4
  fam    <- fam
  par    <- par.temp
  ll     <- loglik2stage_dim4(par = par, fam = fam, gap1 = gap1, gap2 = gap2, gap3 = gap3, gap4 = gap4, status4 = status)

  #
  # adding loglikelihood contributions of clusters of size 3
  #
  gap1   <- data$gap1[(length(data$gap4)+1):length(data$gap3)]
  gap2   <- data$gap2[(length(data$gap4)+1):length(data$gap3)]
  gap3   <- data$gap3[(length(data$gap4)+1):length(data$gap3)]
  status <- data$status3[(length(data$gap4)+1):length(data$gap3)]
  fam    <- c(fam[1:2], fam[4])
  par    <- c(par.temp[1:2], par.temp[4])
  ll     <- ll + loglik2stage_dim3(par = par, fam = fam, gap1 = gap1, gap2 = gap2, gap3 = gap3, status3 = status)

  #
  # adding loglikelihood contributions of clusters of size 2
  #
  gap1   <- data$gap1[(length(data$gap3)+1):length(data$gap2)]
  gap2   <- data$gap2[(length(data$gap3)+1):length(data$gap2)]
  status <- data$status2[(length(data$gap3)+1):length(data$gap2)]
  fam    <- fam[1]
  par    <- par.temp[1]
  ll     <- ll + loglik2stage_dim2(par = par, fam = fam, gap1 = gap1, gap2 = gap2, status2 = status)

  return(ll)

}


# sequential estimation of eta: faster than joint simultaneous estimation
MLE_2stage_seq <- function(data, paras.init){

  res <- rep(NA, 6)
  data.temp <- list()

  #
  # estimation of all pair-copulas in first tree level
  #
  data.temp[[1]] <- cbind(data$gap1[1:length(data$gap2)], data$gap2, data$status2)
  data.temp[[2]] <- cbind(data$gap2[1:length(data$gap3)], data$gap3, data$status3)
  data.temp[[3]] <- cbind(data$gap3[1:length(data$gap4)], data$gap4, data$status4)
  for(i.cop in 1:(paras.init$dim_cop-1)){

    fam2 <- paras.init$RVM$family[4,(4-i.cop)]
    res[i.cop] <- optim(par = init.par(fam2),
                        fn = loglik2stage_dim2,
                        gap1 = data.temp[[i.cop]][, 1],
                        gap2 = data.temp[[i.cop]][, 2],
                        status2 = data.temp[[i.cop]][, 3],
                        fam = fam2,
                        # method = "L-BFGS-B",
                        upper = upper.bounds(fam2),
                        lower = lower.bounds(fam2),
                        control = list(fnscale = -1))$par # loglik is positive)

  }

  #
  # generation of pseudo copula data in second tree level
  #
  data.temp[[1]] <- cbind(BiCopHfunc2(u1 = data$gap1[1:length(data$gap3)],
                                      u2 = data$gap2[1:length(data$gap3)],
                                      family = paras.init$RVM$family[4, 3], par = res[1]),
                          BiCopHfunc1(u1 = data$gap2[1:length(data$gap3)],
                                      u2 = data$gap3,
                                      family = paras.init$RVM$family[4, 2], par = res[2]),
                          data$status3)
  data.temp[[2]] <- cbind(BiCopHfunc2(u1 = data$gap2[1:length(data$gap4)],
                                      u2 = data$gap3[1:length(data$gap4)],
                                      family = paras.init$RVM$family[4, 2], par = res[2]),
                          BiCopHfunc1(u1 = data$gap3[1:length(data$gap4)],
                                      u2 = data$gap4,
                                      family = paras.init$RVM$family[4, 1], par = res[3]),
                          data$status4)

  #
  # estimation of all pair-copulas in second tree level
  #
  for(i.cop in 1:2){

    fam2 <- paras.init$RVM$family[3,(4-1-i.cop)]
    res[paras.init$dim_cop-1+i.cop] <- optim(par = init.par(fam2),
                                             fn = loglik2stage_dim2,
                                             gap1 = data.temp[[i.cop]][, 1],
                                             gap2 = data.temp[[i.cop]][, 2],
                                             status2 = data.temp[[i.cop]][, 3],
                                             fam = fam2,
                                             # method = "L-BFGS-B",
                                             upper = upper.bounds(fam2),
                                             lower = lower.bounds(fam2),
                                             control = list(fnscale = -1))$par

  }

  #
  # generation of pseudo copula data in third tree level
  #
  data.temp_old <- cbind(data.temp[[1]][1:length(data$gap4), 1],
                         data.temp[[1]][1:length(data$gap4), 2],
                         data.temp[[2]][1:length(data$gap4), 1],
                         data.temp[[2]][1:length(data$gap4), 2])
  data.temp     <- cbind(BiCopHfunc2(u1 = data.temp_old[, 1],
                                     u2 = data.temp_old[, 2],
                                     family = paras.init$RVM$family[3, 2],
                                     par = res[4]),
                         BiCopHfunc1(u1 = data.temp_old[, 3],
                                     u2 = data.temp_old[, 4],
                                     family = paras.init$RVM$family[3, 1], par = res[5]),
                         data$status4)

  #
  # estimation of pair-copulas in third tree level
  #
  fam2 <- paras.init$RVM$family[2, 1]
  res[6] <- optim(par = init.par(fam2),
                  fn = loglik2stage_dim2,
                  gap1 = data.temp[, 1],
                  gap2 = data.temp[, 2],
                  status2 = data.temp[, 3],
                  fam = fam2,
                  # method = "L-BFGS-B",
                  upper = upper.bounds(fam2),
                  lower = lower.bounds(fam2),
                  control = list(fnscale = -1))$par

  return(res)

}
