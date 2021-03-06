# rc, cox margin
rc_copula_log_lik_cox <-function(p, x1, x2, indata1, indata2, copula)
{

  if (copula != "Copula2") {
    eta <- p[1]
    beta <- p[2:length(p)] # coefficients
  }

  if (copula == "Copula2") {
    alpha <- p[1]
    kappa <- p[2]
    beta <- p[3:length(p)] # coefficients
  }

  t1<-indata1[, "obs_time"]
  t2<-indata2[, "obs_time"]
  c1<-indata1[, "status"]
  c2<-indata2[, "status"]

  t <- c(t1, t2)
  c <- c(c1, c2)
  x1 <- matrix(x1, nrow = length(t1))
  x2 <- matrix(x2, nrow = length(t2))
  x12 <- rbind(x1, x2)
  tab <- data.frame(table(t[c == 1]))
  y <- as.numeric(levels(tab[, 1]))[tab[, 1]] # sorted unique event times
  d <- tab[, 2] # number of events at each unique time (in case of tied events)

  Lambda <- 0
  Lambda <- 0

  # Create hazard function and calculate observed cumulative hazard and instanenuous hazard
  h0 <- d/sapply(y,
                 function(x) {
                   sum(exp(as.matrix(x12[t >= x, ]) %*% beta))
                 })
  lambda <- rowSums(sapply(y, function(x) t==x) %*% h0 * c)
  # # new_c: ones who are censored but have the same observation times as event times
  # new_c <- c
  # new_c[which(t %in% y)] <- 1
  # lambda <- rowSums(sapply(y, function(x) t==x) %*% h0 * new_c)
  Lambda <- as.numeric(sapply(y, function(x) t>=x) %*% h0)

  lambda1 <- lambda[1:length(t1)]
  lambda2 <- lambda[-(1:length(t1))]
  Lambda1 <- Lambda[1:length(t1)]
  Lambda2 <- Lambda[-(1:length(t1))]

  # tab1 <- data.frame(table(t1[c1 == 1]))
  # y1 <- as.numeric(levels(tab1[, 1]))[tab1[, 1]] # sorted unique event times
  # d1 <- tab1[, 2] # number of events at each unique time (in case of tied events)
  #
  # tab2 <- data.frame(table(t2[c2 == 1]))
  # y2 <- as.numeric(levels(tab2[, 1]))[tab2[, 1]]
  # d2 <- tab2[, 2]
  #
  # Lambda1 <- 0
  # Lambda2 <- 0
  # lambda1 <- 0
  # lambda2 <- 0
  #
  # # Create hazard function and calculate observed cumulative hazard and instanenuous hazard
  # h01<-d1/sapply(y1,
  #                function(x) {
  #                  sum(exp(as.matrix(x1[t1 >= x, ]) %*% beta))
  #                })
  # lambda1 <- rowSums(sapply(y1, function(x) t1==x) %*% h01 * c1)
  # Lambda1 <- as.numeric(sapply(y1, function(x) t1>=x) %*% h01)
  #
  # h02 <- d2/sapply(y2,
  #                  function(x) {
  #                    sum(exp(as.matrix(x1[t2 >= x, ]) %*% beta))
  #                  })
  # lambda2 <- rowSums(sapply(y2, function(x) t2 == x) %*% h02 * c2)
  # Lambda2 <- as.numeric(sapply(y2, function(x) t2 >= x) %*% h02)

  # survival and density
  u1 <- exp(-Lambda1 * exp(x1 %*% beta))
  u2 <- exp(-Lambda2 * exp(x2 %*% beta))
  u1_t1 <- -u1 * lambda1 * exp(x1 %*% beta)
  u2_t2 <- -u2 * lambda2 * exp(x2 %*% beta)

  ### copula models ######
  if (copula == "Copula2") {

    # Copula2 Copula function for joint distribution probability
    C_val <- (1 + ((u1^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    c_u1_val <- C_val^(1+1/kappa) * (1 + ((u2^(-1/kappa)-1)/(u1^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u1^(-1/kappa-1) # wrt to u1
    c_u2_val <- C_val^(1+1/kappa) * (1 + ((u1^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2
    c_val <- (1+1/kappa) * C_val^(1/kappa) * c_u2_val * (1 + ((u2^(-1/kappa)-1)/(u1^(-1/kappa)-1))^(1/alpha) )^(alpha-1) * u1^(-1-1/kappa) +
      (-1+1/alpha) * (1/kappa) * C_val^(1+1/kappa) * (u1*u2)^(-1-1/kappa) * (1 + ((u2^(-1/kappa)-1)/(u1^(-1/kappa)-1))^(1/alpha) )^(alpha-2) * (u2^(-1/kappa)-1)^(-1+1/alpha) * (u1^(-1/kappa)-1)^(-1/alpha)

  }



  if (copula == "Joe") {

    # Joe joint and density functions
    C_val <- 1 - ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta)
    c_u1_val <- ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u1)^(eta-1) - (1-u1)^(eta-1)*(1-u2)^eta)
    c_u2_val <- ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1)^eta)
    c_val <- ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 2) * (1/eta - 1) * ((1-u1)^(eta-1) - (1-u1)^(eta-1)*(1-u2)^eta) *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1)^eta) * (-eta)+
      ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 1) * (eta * ((1-u1)^(eta-1)) * ((1-u2)^(eta-1)) )

  }


  if (copula == "Gumbel") {

    c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-
      gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
    c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

  }


  if (copula == "AMH") {

    # AHM joint and desity functions
    c_val <- ((eta^2)*(1-u1)*(1-u2) + eta*(u1*u2+u1+u2-2) + 1)/(1-eta*(1-u1)*(1-u2))^3
    c_u1_val<- u2*(1-eta*(1-u2))/(1-eta*(1-u1)*(1-u2))^2
    c_u2_val<- u1*(1-eta*(1-u1))/(1-eta*(1-u1)*(1-u2))^2
    C_val<- u1*u2/(1-eta*(1-u1)*(1-u2))

  }


  if (copula == "Frank") {

    # Frank Copula function for joint distribution probability
    c_val <- (-1) * eta*exp((-1) * eta*u1)*exp((-1) * eta*u2)/((exp((-1) * eta)-1)+(exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)) -
      (-1) * eta*exp((-1) * eta*u1)*exp((-1) * eta*u2)*(exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/((exp((-1) * eta)-1)+(exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1))^2
    c_u1_val <- (exp((-1) * eta*u2)-1)*exp((-1) * eta*u1)/((1 + (exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    c_u2_val <- (exp((-1) * eta*u1)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    C_val <- (1/(-1) * eta) * log(1 + (exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    # Clayton Copula function for joint distribution probability
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)

  }


  C_val <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), C_val, 1)
  term1 <- log(abs(C_val))

  term2 <- c_u1_val * (-u1_t1)
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- (c_u2_val) * (-u2_t2)
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
  term3 <- log(abs(term3))

  term4 <- c_val * u1_t1 * u2_t2
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
  term4 <- log(abs(term4))

  logL <- (-1) * sum(term1 + term2 + term3 + term4)
  return(logL)
}

# rc, parametric margins
rc_copula_log_lik <-function(p, x1, x2,indata1,indata2, quantiles = NULL, copula, m.dist)
{

  if (m.dist == "Weibull") {

    lambda <- p[1]
    k <- p[2]

    if (copula != "Copula2") {

      beta<- p[3:(length(p)-1)] # coefficients
      eta <- p[length(p)]

    }

    if (copula == "Copula2") {

      beta<- p[3:(length(p)-2)] # coefficients
      alpha <- p[(length(p)-1)]
      kappa <- p[length(p)]

    }


    t1<-indata1[,"obs_time"]
    t2<-indata2[,"obs_time"]

    u1<-exp(-(t1/lambda)^k*exp(x1%*%beta))
    u2<-exp(-(t2/lambda)^k*exp(x2%*%beta))
    u1_t1<- -u1*k*(1/lambda)*(t1/lambda)^(k-1)*exp(x1%*%beta)
    u2_t2<- -u2*k*(1/lambda)*(t2/lambda)^(k-1)*exp(x2%*%beta)

  }


  if (m.dist == "Loglogistic") {

    lambda <- p[1]
    k <- p[2]

    if (copula != "Copula2") {

      beta<- p[3:(length(p)-1)] # coefficients
      eta <- p[length(p)]

    }

    if (copula == "Copula2") {

      beta<- p[3:(length(p)-2)] # coefficients
      alpha <- p[(length(p)-1)]
      kappa <- p[length(p)]

    }


    t1<-indata1[,"obs_time"]
    t2<-indata2[,"obs_time"]

    # loglogistic surv and density
    u1 <- (1+((t1/lambda)^k)*exp(x1%*%beta))^(-1)
    u2 <- (1+((t2/lambda)^k)*exp(x2%*%beta))^(-1)
    u1_t1 <- -1*((u1)^(2)) * (k/lambda) * ((t1/lambda)^(k-1)) * exp(x1%*%beta)
    u2_t2 <- -1*((u2)^(2)) * (k/lambda) * ((t2/lambda)^(k-1)) * exp(x2%*%beta)

  }


  if (m.dist == "Gompertz") {

    a <- p[1]
    b <- p[2]

    if (copula != "Copula2") {

      beta<- p[3:(length(p)-1)] # coefficients
      eta <- p[length(p)]

    }

    if (copula == "Copula2") {

      beta<- p[3:(length(p)-2)] # coefficients
      alpha <- p[(length(p)-1)]
      kappa <- p[length(p)]

    }

    t1<-indata1[,"obs_time"]
    t2<-indata2[,"obs_time"]
    # x1<-as.matrix(indata1[,var_list])
    # x2<-as.matrix(indata2[,var_list])

    u1<-exp(b/a*(1-exp(a*t1))*exp(x1%*%beta))
    u2<-exp(b/a*(1-exp(a*t2))*exp(x2%*%beta))
    u1_t1<- -u1*b*exp(a*t1)*exp(x1%*%beta)
    u2_t2<- -u2*b*exp(a*t2)*exp(x2%*%beta)

  }


  if (copula == "Copula2") {

    # Copula2 Copula function for joint distribution probability
    C_val <- (1 + ((u1^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    c_u1_val <- C_val^(1+1/kappa) * (1 + ((u2^(-1/kappa)-1)/(u1^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u1^(-1/kappa-1) # wrt to u1
    c_u2_val <- C_val^(1+1/kappa) * (1 + ((u1^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2
    c_val <- (1+1/kappa) * C_val^(1/kappa) * c_u2_val * (1 + ((u2^(-1/kappa)-1)/(u1^(-1/kappa)-1))^(1/alpha) )^(alpha-1) * u1^(-1-1/kappa) +
      (-1+1/alpha) * (1/kappa) * C_val^(1+1/kappa) * (u1*u2)^(-1-1/kappa) * (1 + ((u2^(-1/kappa)-1)/(u1^(-1/kappa)-1))^(1/alpha) )^(alpha-2) * (u2^(-1/kappa)-1)^(-1+1/alpha) * (u1^(-1/kappa)-1)^(-1/alpha)

  }



  if (copula == "Joe") {

    # Joe joint and density functions
    C_val <- 1 - ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta)
    c_u1_val <- ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u1)^(eta-1) - (1-u1)^(eta-1)*(1-u2)^eta)
    c_u2_val <- ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1)^eta)
    c_val <- ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 2) * (1/eta - 1) * ((1-u1)^(eta-1) - (1-u1)^(eta-1)*(1-u2)^eta) *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1)^eta) * (-eta)+
      ((1-u1)^eta + (1-u2)^eta - ((1-u1)^eta)*((1-u2)^eta) )^(1/eta - 1) * (eta * ((1-u1)^(eta-1)) * ((1-u2)^(eta-1)) )

  }


  if (copula == "Gumbel") {

    c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-
      gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
    c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

  }


  if (copula == "AMH") {

    # AHM joint and desity functions
    c_val <- ((eta^2)*(1-u1)*(1-u2) + eta*(u1*u2+u1+u2-2) + 1)/(1-eta*(1-u1)*(1-u2))^3
    c_u1_val<- u2*(1-eta*(1-u2))/(1-eta*(1-u1)*(1-u2))^2
    c_u2_val<- u1*(1-eta*(1-u1))/(1-eta*(1-u1)*(1-u2))^2
    C_val<- u1*u2/(1-eta*(1-u1)*(1-u2))

  }


  if (copula == "Frank") {

    # Frank Copula function for joint distribution probability
    c_val <- (-1) * eta*exp((-1) * eta*u1)*exp((-1) * eta*u2)/((exp((-1) * eta)-1)+(exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)) -
      (-1) * eta*exp((-1) * eta*u1)*exp((-1) * eta*u2)*(exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/((exp((-1) * eta)-1)+(exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1))^2
    c_u1_val <- (exp((-1) * eta*u2)-1)*exp((-1) * eta*u1)/((1 + (exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    c_u2_val <- (exp((-1) * eta*u1)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    C_val <- (1/(-1) * eta) * log(1 + (exp((-1) * eta*u1)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    # Clayton Copula function for joint distribution probability
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)

  }

  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), C_val, 1)
  term1 <- log(abs(term1))

  term2 <- c_u1_val * (-u1_t1)
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
  term2 <- log(abs(term2))
  # term2 <- log(c_u1_val)+log(-u1_t1)
  # term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 0)

  term3 <- (c_u2_val) * (-u2_t2)
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
  term3 <- log(abs(term3))
  # term3 <- log(c_u2_val)+log(-u2_t2)
  # term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 0)

  term4 <- c_val * u1_t1 * u2_t2
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
  term4[term4 < 0] <- 1
  term4 <- log(abs(term4))
  # term4 <- log(c_val)+log(-u1_t1)+log(-u2_t2)
  # term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 0)

  logL<-(-1)*sum( term1 + term2 + term3 + term4 )
  return(logL)
}


# ic, sieve
ic_copula_log_lik_sieve <- function(para, p, m, x1, x2, bl1, br1, bl2, br2, indata1, indata2, r, copula) {


  if (copula != "Copula2") {

    beta<-para[1:p] #1:3, p=3
    phi<-para[(p+1):(p+1+m)] #4:7, m=3
    eta<-para[p+1+m+1] #8

  }

  if (copula == "Copula2") {

    beta<-para[1:p]
    phi<-para[(p+1):(p+1+m)] #4:7, m=3
    alpha<-para[p+1+m+1] #8
    kappa<-para[p+1+m+1+1] #9

  }


  ep<-cumsum(exp(phi))

  gLl1<-G(exp(x1%*%beta)*(bl1%*%ep),r) #left eye, left end
  gLr1<-G(exp(x1%*%beta)*(br1%*%ep),r) #left eye, right end
  gLl2<-G(exp(x2%*%beta)*(bl2%*%ep),r) #right eye, left end
  gLr2<-G(exp(x2%*%beta)*(br2%*%ep),r)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2_left = exp(-gLl2)
  u2_right = exp(-gLr2)


  ### copula models ###

  if (copula == "Copula2") {

    # copula2 Copula function for joint distribution probability
    C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2_left^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_2 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2_right^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_3 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2_left^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_4 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2_right^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)

  }



  if (copula == "Joe") {

    # Joe Copula function for joint distribution probability
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2_left)^(eta) - (1-u1_left)^(eta)*(1-u2_left)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_left)^(eta) + (1-u2_right)^(eta) - (1-u1_left)^(eta)*(1-u2_right)^(eta) )^(1/eta)
    C_val_3 <- 1 - ( (1-u1_right)^(eta) + (1-u2_left)^(eta) - (1-u1_right)^(eta)*(1-u2_left)^(eta) )^(1/eta)
    C_val_4 <- 1 - ( (1-u1_right)^(eta) + (1-u2_right)^(eta) - (1-u1_right)^(eta)*(1-u2_right)^(eta) )^(1/eta)

  }


  if (copula == "Gumbel") {

    # Gumbel Copula function for joint distribution probability
    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_left))^eta + (-log(u2_right))^eta)^(1/eta))
    C_val_3 <- exp(-((-log(u1_right))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_4 <- exp(-((-log(u1_right))^eta + (-log(u2_right))^eta)^(1/eta))

  }


  if (copula == "AMH") {

    # AMH Copula function for joint distribution probability
    C_val_1 <- u1_left * u2_left /(1 - eta * (1-u1_left) * (1-u2_left))
    C_val_2 <- u1_left * u2_right /(1 - eta * (1-u1_left) * (1-u2_right))
    C_val_3 <- u1_right * u2_left /(1 - eta * (1-u1_right) * (1-u2_left))
    C_val_4 <- u1_right * u2_right /(1 - eta * (1-u1_right) * (1-u2_right))

  }


  if (copula == "Frank") {

    # Frank Copula function for joint distribution probability
    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2_left)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2_right)-1)/(exp((-1) * eta)-1))
    C_val_3 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2_left)-1)/(exp((-1) * eta)-1))
    C_val_4 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2_right)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    # Copula function for joint distribution probability
    C_val_1 <-(u1_left^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_left^(-eta)+u2_right^(-eta)-1)^(-1/eta)
    C_val_3 <-(u1_right^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_4 <-(u1_right^(-eta)+u2_right^(-eta)-1)^(-1/eta)

  }


  # Use Copula functions to write each block of likelihood function
  term1 <- log(abs(C_val_1 - C_val_2 - C_val_3 + C_val_4))
  term1 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term1, 0)

  term2 <- log(abs(C_val_1 - C_val_3))
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 0)

  term3 <- log(abs(C_val_1 - C_val_2))
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 0)

  term4 <- log(abs(C_val_1))
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term4, 0)

  logL<-(-1)*sum( term1 + term2 + term3 + term4 )

  return(logL)

}


# ic, parametric margins
ic_copula_log_lik_param <- function(para, p, x1, x2, t1_left, t1_right, t2_left, t2_right, indata1, indata2, copula, m.dist) {


  ### Weibull ###
  if (m.dist == "Weibull") {

    lambda<-para[1]
    k<-para[2]

    if (copula != "Copula2") {

      beta<-para[3:(length(para)-1)]
      eta <- para[length(para)]

    }

    if (copula == "Copula2") {

      beta<-para[3:(p+2)]
      alpha <- para[length(para)-1]
      kappa <- para[length(para)]
    }

    # Survival functions with weibull margins under PH assumption
    u1_left<-exp(-(t1_left/lambda)^k*exp(x1%*%beta))
    u1_right<-exp(-(t1_right/lambda)^k*exp(x1%*%beta))
    u2_left<-exp(-(t2_left/lambda)^k*exp(x2%*%beta))
    u2_right<-exp(-(t2_right/lambda)^k*exp(x2%*%beta))

  }



  ### Loglog ###

  if (m.dist == "Loglogistic") {

    # loglog parameters
    lambda<-para[1]
    k<-para[2]

    if (copula != "Copula2") {

      beta<-para[3:(p+2)]
      eta <- para[length(para)]

    }

    if (copula == "Copula2") {

      beta<-para[3:(p+2)]
      alpha <- para[length(para)-1]
      kappa <- para[length(para)]

    }


    # Survival functions for loglogistic distribution under PO assumption
    u1_left<- (1+((t1_left/lambda)^k)*exp(x1%*%beta))^(-1)
    u1_right<- (1+((t1_right/lambda)^k)*exp(x1%*%beta))^(-1)
    u2_left<- (1+((t2_left/lambda)^k)*exp(x2%*%beta))^(-1)
    u2_right<- (1+((t2_right/lambda)^k)*exp(x2%*%beta))^(-1)


  }



  ### Gompertz ###


  if (m.dist == "Gompertz") {

    a<-para[1]
    b<-para[2]

    if (copula != "Copula2") {

      beta<-para[3:(length(para)-1)]
      eta <- para[length(para)]

    }

    if (copula == "Copula2") {

      beta<-para[3:(p+2)]
      alpha <- para[length(para)-1]
      kappa <- para[length(para)]


    }

    # Survival functions with gompertz margins under PH assumption
    u1_left<-exp(b/a*(1-exp(a*t1_left))*exp(x1%*%beta))
    u1_right<-exp(b/a*(1-exp(a*t1_right))*exp(x1%*%beta))
    u2_left<-exp(b/a*(1-exp(a*t2_left))*exp(x2%*%beta))
    u2_right<-exp(b/a*(1-exp(a*t2_right))*exp(x2%*%beta))

  }



  ### copula models ###

  if (copula == "Copula2") {

    # copula2 Copula function for joint distribution probability
    C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2_left^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_2 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2_right^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_3 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2_left^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_4 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2_right^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)

  }



  if (copula == "Joe") {

    # Joe Copula function for joint distribution probability
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2_left)^(eta) - (1-u1_left)^(eta)*(1-u2_left)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_left)^(eta) + (1-u2_right)^(eta) - (1-u1_left)^(eta)*(1-u2_right)^(eta) )^(1/eta)
    C_val_3 <- 1 - ( (1-u1_right)^(eta) + (1-u2_left)^(eta) - (1-u1_right)^(eta)*(1-u2_left)^(eta) )^(1/eta)
    C_val_4 <- 1 - ( (1-u1_right)^(eta) + (1-u2_right)^(eta) - (1-u1_right)^(eta)*(1-u2_right)^(eta) )^(1/eta)

  }


  if (copula == "Gumbel") {

    # Gumbel Copula function for joint distribution probability
    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_left))^eta + (-log(u2_right))^eta)^(1/eta))
    C_val_3 <- exp(-((-log(u1_right))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_4 <- exp(-((-log(u1_right))^eta + (-log(u2_right))^eta)^(1/eta))

  }


  if (copula == "AMH") {

    # AMH Copula function for joint distribution probability
    C_val_1 <- u1_left * u2_left /(1 - eta * (1-u1_left) * (1-u2_left))
    C_val_2 <- u1_left * u2_right /(1 - eta * (1-u1_left) * (1-u2_right))
    C_val_3 <- u1_right * u2_left /(1 - eta * (1-u1_right) * (1-u2_left))
    C_val_4 <- u1_right * u2_right /(1 - eta * (1-u1_right) * (1-u2_right))

  }


  if (copula == "Frank") {

    # Frank Copula function for joint distribution probability
    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2_left)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2_right)-1)/(exp((-1) * eta)-1))
    C_val_3 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2_left)-1)/(exp((-1) * eta)-1))
    C_val_4 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2_right)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    # Copula function for joint distribution probability
    C_val_1 <-(u1_left^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_left^(-eta)+u2_right^(-eta)-1)^(-1/eta)
    C_val_3 <-(u1_right^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_4 <-(u1_right^(-eta)+u2_right^(-eta)-1)^(-1/eta)

  }


  # Use Copula functions to write each block of likelihood function
  term1 <- log(abs(C_val_1 - C_val_2 - C_val_3 + C_val_4))
  term1 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term1, 0)

  term2 <- log(abs(C_val_1 - C_val_3))
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 0)

  term3 <- log(abs(C_val_1 - C_val_2))
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 0)

  term4 <- log(abs(C_val_1))
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term4, 0)

  logL<-(-1)*sum( term1 + term2 + term3 + term4 )
  return(logL)

}


### marginal sieve model
### also used in step 1a when fitting a copula-based sieve model
estimate_sieve_step1a <- function(para, p, m, x1, x2, bl1, br1, bl2, br2, indata1, indata2,r) { #total 7 parameters = 3beta + 4polynomial

  beta<-para[1:p] #p=3
  phi<-para[(p+1):(p+1+m)] #4:7,m=3

  ep<-cumsum(exp(phi))

  gLl1<-G(exp(x1%*%beta)*(bl1%*%ep),r) #left eye, left end
  gLr1<-G(exp(x1%*%beta)*(br1%*%ep),r)
  gLl2<-G(exp(x2%*%beta)*(bl2%*%ep),r) #right eye, left end
  gLr2<-G(exp(x2%*%beta)*(br2%*%ep),r)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2_left = exp(-gLl2)
  u2_right = exp(-gLr2)

  # joint distribution probability under independence
  C_val_1 <- u1_left*u2_left
  C_val_2 <- u1_left*u2_right
  C_val_3 <- u1_right*u2_left
  C_val_4 <- u1_right*u2_right

  term1 <- log(abs(C_val_1 - C_val_2 - C_val_3 + C_val_4))
  term1 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term1, 0)

  term2 <- log(abs(C_val_1 - C_val_3))
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 0)

  term3 <- log(abs(C_val_1 - C_val_2))
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 0)

  term4 <- log(abs(C_val_1))
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term4, 0)


  logL <- (-1)*sum( term1 + term2 + term3 + term4 )


  return(logL)
}


