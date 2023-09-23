# Berstein polynomials derivatives
bern_derivative <- function(j,m,l,u,t){

  if (j != 0) {
    # b = m*bern(j-1, m-1, l, u, t) - m*bern(j, m-1, l, u, t)
    b = (m/(u - l)) * (bern(j-1, m-1, l, u, t)-bern(j, m-1, l, u, t))
  } else if (j == 0) {
    b = -1 * m * (u - t)^(m - 1)/((u - l)^m)
  }
  return(b)
}


# ic, semi-competing, sieve
data_process_scmprisk_ic_sp_LT_A1A2 <- function(data, var_list, l1, u1, m1, l2, u2, m2) {

  # replace Inf by constant u
  data$Right[data$status==0] <- u1

  indata1 <- data[,c("Left","Right","status","A1",var_list)]
  indata2 <- data[,c("timeD","statusD","A2",var_list)]
  t1_left <- indata1[,"Left"]
  t1_right <- indata1[,"Right"]
  t2 <- indata2[,"timeD"]
  A1 <- indata1[,"A1"]
  A2 <- indata2[,"A2"]

  dim_m <- dim(as.matrix(indata1[,var_list]))
  n <- dim_m[1]

  # x1, x2, new var_list
  tmp1 <- var_transform(indata1, var_list)
  tmp2 <- var_transform(indata2, var_list)

  x1 <- tmp1$x
  x2 <- tmp2$x
  var_list <- tmp1$var_list

  # matix
  x1 <- as.matrix(x1,dim_m[1])
  x2 <- as.matrix(x2,dim_m[1])
  p <- dim(x1)[2]

  # BP
  bl1 <-matrix(0,nrow = n,ncol = m1+1) # AD, left end
  br1 <-matrix(0,nrow = n,ncol = m1+1) # AD, right end
  b2 <- matrix(0,nrow = n,ncol = m2+1) # terminal
  b1_A <- matrix(0,nrow = n,ncol = m1+1) # AD, left truncation
  b2_A <- matrix(0,nrow = n,ncol = m2+1) # terminal, left truncation

  for (i in 0:m1) {
    bl1[,(i+1)] <- bern(i,m1,l1,u1,t1_left)
    br1[,(i+1)] <- bern(i,m1,l1,u1,t1_right)
    b1_A[,(i+1)] <- bern(i,m1,l1,u1,A1)
  }

  for (i in 0:m2) {
    b2[,(i+1)] <- bern(i,m2,l2,u2,t2)
    b2_A[,(i+1)] <- bern(i,m2,l2,u2,A2)
  }

  # BP derivatives
  bl1_d <-matrix(0,nrow = n,ncol = m1+1) # AD, left end
  br1_d <-matrix(0,nrow = n,ncol = m1+1) # AD, right end
  b2_d  <- matrix(0,nrow = n,ncol = m2+1) # terminal

  for (i in 0:m1) {
    bl1_d[,(i+1)] <- bern_derivative(i,m1,l1,u1,t1_left)
    br1_d[,(i+1)] <- bern_derivative(i,m1,l1,u1,t1_right)
  }

  for (i in 0:m2) {
    b2_d[,(i+1)] <- bern_derivative(i,m2,l2,u2,t2)
  }

  return(list(indata1=indata1, indata2=indata2, t1_left=t1_left, t1_right=t1_right, t2=t2,
              n=n, p=p, x1=x1, x2=x2, var_list = var_list, bl1=bl1, br1=br1, b2=b2, bl1_d=bl1_d, br1_d=br1_d, b2_d=b2_d,
              A1=A1, A2=A2, b1_A=b1_A, b2_A=b2_A))

}

### marginal sieve model under right censoring and left truncation
### also used in step 1a when fitting a semi-parametric transformation model in LTRC data
# omit LT for now
estimate_sieve_step1a_scmprisk_rc_LT <- function(para, p2, m2, x2, b2, b2_d, b2_A, indata2, r2) { #total 7 parameters = 3beta + 4polynomial

  beta2<-para[1:p2] #p=3
  phi2<-para[(p2+1):(p2+1+m2)] #4:7,m=3

  ep2<-cumsum(exp(phi2))

  # survival probability
  gL2 <- G(exp(x2%*%beta2)*(b2%*%ep2),r2)
  u2 = exp(-gL2)

  # derivative of survival (= negative density)
  u2_t2 <- u2 * (-1) * pG(exp(x2%*%beta2)*(b2%*%ep2),r2) * (b2_d %*% ep2) * exp(x2 %*% beta2)

  # Left truncation
  gL2_A <- G(exp(x2%*%beta2)*(b2_A%*%ep2),r2)
  u2_A = exp(-gL2_A)

  # likelihood
  term1 <- ifelse((indata2[,"statusD"] == 0), u2, 1)
  term1 <- log(abs(term1))

  term2 <- ifelse((indata2[,"statusD"] == 1), -1 * u2_t2, 1)
  term2 <- log(abs(term2))

  # left truncation
  term3 <- ifelse((indata2[, "A"] == 0), 1, u2_A) # A = 0 means S(A)=1, term3 = log1 = 0
  term3 <- log(abs(term3))

  # # left truncation
  # term3 <- ifelse((indata2[, "A"] == 0), 1, u2_A) # A = 0 means S(A)=1, term3 = log1 = 0
  # term3 <- -1*log(abs(term3))
  #
  # term1 <- ifelse((indata2[,"statusD"] == 0), u2, 1)
  # term1 <- log(abs(term1))
  #
  # term2 <- ifelse((indata2[,"statusD"] == 0), -1 * u2_t2, 1)
  # term2 <- log(abs(term2))



  # logL <- (-1)*sum( term1 + term2 - term3)
  logL <- (-1)*sum( term1 + term2) # simulations


  return(logL)
}


# consider LT for control methods
estimate_sieve_step1a_scmprisk_rc_LT_2 <- function(para, p2, m2, x2, b2, b2_d, b2_A, indata2, r2) { #total 7 parameters = 3beta + 4polynomial

  beta2<-para[1:p2] #p=3
  phi2<-para[(p2+1):(p2+1+m2)] #4:7,m=3

  ep2<-cumsum(exp(phi2))

  # survival probability
  gL2 <- G(exp(x2%*%beta2)*(b2%*%ep2),r2)
  u2 = exp(-gL2)

  # derivative of survival (= negative density)
  u2_t2 <- u2 * (-1) * pG(exp(x2%*%beta2)*(b2%*%ep2),r2) * (b2_d %*% ep2) * exp(x2 %*% beta2)

  # Left truncation
  gL2_A <- G(exp(x2%*%beta2)*(b2_A%*%ep2),r2)
  u2_A = exp(-gL2_A)

  # likelihood
  term1 <- ifelse((indata2[,"statusD"] == 0), u2, 1)
  term1 <- log(abs(term1))

  term2 <- ifelse((indata2[,"statusD"] == 1), -1 * u2_t2, 1)
  term2 <- log(abs(term2))

  # left truncation
  term3 <- ifelse((indata2[, "A"] == 0), 1, u2_A) # A = 0 means S(A)=1, term3 = log1 = 0
  term3 <- log(abs(term3))

  # # left truncation
  # term3 <- ifelse((indata2[, "A"] == 0), 1, u2_A) # A = 0 means S(A)=1, term3 = log1 = 0
  # term3 <- -1*log(abs(term3))
  #
  # term1 <- ifelse((indata2[,"statusD"] == 0), u2, 1)
  # term1 <- log(abs(term1))
  #
  # term2 <- ifelse((indata2[,"statusD"] == 0), -1 * u2_t2, 1)
  # term2 <- log(abs(term2))



  logL <- (-1)*sum( term1 + term2 - term3, log=T)
  # logL <- (-1)*sum( term1 + term2) # simulations


  return(logL)
}


# marginal, ic, semi-parametric margins, first step omitting left truncation
estimate_sieve_step1a_scmprisk_ic_LT_1 <- function(para, p1, m1, x1, bl1, br1, b1_A, indata1, r1) { #total 7 parameters = 3beta + 4polynomial

  beta1<-para[1:p1] #p=3
  phi1<-para[(p1+1):(p1+1+m1)] #4:7,m=3

  ep1<-cumsum(exp(phi1))

  gLl1<-G(exp(x1%*%beta1)*(bl1%*%ep1),r1) # AD, left end
  gLr1<-G(exp(x1%*%beta1)*(br1%*%ep1),r1)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)

  # Left truncation
  gL1_A <- G(exp(x1%*%beta1)*(b1_A%*%ep1),r1)
  u1_A = exp(-gL1_A)

  # likelihood
  term1 <- log(abs(u1_left - u1_right))
  term1 <- ifelse((indata1[,"status"] == 1), term1, 0)

  term2 <- log(abs(u1_left))
  term2 <- ifelse((indata1[,"status"] == 0), term2, 0)

  # # left truncation
  # term3 <- ifelse((indata1[, "A"] == 0), 1, u1_A) # A = 0 means S(A)=1, term3 = log1 = 0
  # term3 <- log(abs(term3))

  # logL <- (-1)*sum(term1 + term2 - term3)
  logL <- (-1)*sum(term1 + term2)


  return(logL)
}


# marginal, ic, semi-parametric margins, left truncation
estimate_sieve_step1a_scmprisk_ic_LT_2 <- function(para, p1, m1, x1, bl1, br1, b1_A, indata1, r1) { #total 7 parameters = 3beta + 4polynomial

  beta1<-para[1:p1] #p=3
  phi1<-para[(p1+1):(p1+1+m1)] #4:7,m=3

  ep1<-cumsum(exp(phi1))

  gLl1<-G(exp(x1%*%beta1)*(bl1%*%ep1),r1) # AD, left end
  gLr1<-G(exp(x1%*%beta1)*(br1%*%ep1),r1)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)

  # Left truncation
  gL1_A <- G(exp(x1%*%beta1)*(b1_A%*%ep1),r1)
  u1_A = exp(-gL1_A)

  # likelihood
  term1 <- log(abs(u1_left - u1_right))
  term1 <- ifelse((indata1[,"status"] == 1), term1, 0)

  term2 <- log(abs(u1_left))
  term2 <- ifelse((indata1[,"status"] == 0), term2, 0)

  # left truncation
  term3 <- ifelse((indata1[, "A"] == 0), 1, u1_A) # A = 0 means S(A)=1, term3 = log1 = 0
  term3 <- log(abs(term3))

  logL <- (-1)*sum(term1 + term2 - term3)

  return(logL)
}


# general LT, sieve, omitting LT
ic_scmprisk_copula_log_lik_sieve_pseudo_LT_1_copula2 <- function(p, fitted, x1, x2, t1_left, t1_right, t2, indata1, indata2, bl1, br1, b2, b2_d, b1_A, b2_A, m1, m2, r1, r2, p1, p2, quantiles = NULL, copula)
{

  if (copula != "Copula2") {

    eta <- (p[1]) # anti-log
    phi1 <- p[2:(1+1+m1)]
    beta1 <- p[(1+1+m1+1):length(p)]
    ep1 <- cumsum(exp(phi1))


    phi2 <- fitted[1:(1+m2)]
    beta2 <- fitted[(1+m2+1):length(fitted)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  if (copula == "Copula2") {

    alpha <- exp(p[1])/(1+exp(p[1])) # anti-log
    kappa <- exp(p[2]) # anti-log
    phi1 <- p[3:(2+1+m1)]
    beta1 <- p[(2+1+m1+1):length(p)]
    ep1 <- cumsum(exp(phi1))

    phi2 <- fitted[1:(1+m2)]
    beta2 <- fitted[(1+m2+1):length(fitted)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  # survival probability
  gLl1<-G(exp(x1%*%beta1)*(bl1%*%ep1),r1) # AD, left end
  gLr1<-G(exp(x1%*%beta1)*(br1%*%ep1),r1)
  gL2 <- G(exp(x2%*%beta2)*(b2%*%ep2),r2)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2 = exp(-gL2)

  # derivative of survival (= negative density)
  u2_t2 <- u2 * (-1) * pG(exp(x2%*%beta2)*(b2%*%ep2),r2) * (b2_d %*% ep2) * exp(x2 %*% beta2)

  # Left truncation
  gL1_A <- G(exp(x1%*%beta1)*(b1_A%*%ep1),r1)
  u1_A = exp(-gL1_A)
  gL2_A <- G(exp(x2%*%beta2)*(b2_A%*%ep2),r2)
  u2_A = exp(-gL2_A)

  if (copula == "Copula2") {

    # Copula2 Copula function for joint distribution probability
    C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_2 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    c_u2_val_1 <- C_val_1^(1+1/kappa) * (1 + ((u1_left^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2
    c_u2_val_2 <- C_val_2^(1+1/kappa) * (1 + ((u1_right^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2

    C_val_A <- (1 + ((u1_A^(-1/kappa)-1)^(1/alpha) + (u2_A^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)

  }



  if (copula == "Joe") {

    # Joe joint and density functions
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2)^(eta) - (1-u1_left)^(eta)*(1-u2)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_right)^(eta) + (1-u2)^(eta) - (1-u1_right)^(eta)*(1-u2)^(eta) )^(1/eta)
    c_u2_val_1 <- ((1-u1_left)^eta + (1-u2)^eta - ((1-u1_left)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_left)^eta)
    c_u2_val_2 <- ((1-u1_right)^eta + (1-u2)^eta - ((1-u1_right)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_right)^eta)

    C_val_A <- 1 - ( (1-u1_A)^(eta) + (1-u2_A)^(eta) - (1-u1_A)^(eta)*(1-u2_A)^(eta) )^(1/eta)

  }


  if (copula == "Gumbel") {

    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_right))^eta + (-log(u2))^eta)^(1/eta))
    c_u2_val_1 <- gh_F(u1_left,u2,eta)*((-log(u1_left))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    c_u2_val_2 <- gh_F(u1_right,u2,eta)*((-log(u1_right))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2

    C_val_A <- exp(-((-log(u1_A))^eta + (-log(u2_A))^eta)^(1/eta))

  }


  if (copula == "AMH") {

    # AHM joint and desity functions
    C_val_1 <- u1_left * u2 /(1 - eta * (1-u1_left) * (1-u2))
    C_val_2 <- u1_right * u2 /(1 - eta * (1-u1_right) * (1-u2))
    c_u2_val_1 <- u1_left*(1-eta*(1-u1_left))/(1-eta*(1-u1_left)*(1-u2))^2
    c_u2_val_2 <- u1_right*(1-eta*(1-u1_right))/(1-eta*(1-u1_right)*(1-u2))^2

    C_val_A <- u1_A * u2_A /(1 - eta * (1-u1_A) * (1-u2_A))

  }


  if (copula == "Frank") {

    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    c_u2_val_1 <- (exp((-1) * eta*u1_left)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    c_u2_val_2 <- (exp((-1) * eta*u1_right)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))

    C_val_A <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_A)-1)*(exp((-1) * eta*u2_A)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    C_val_1 <-(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta)
    c_u2_val_1 <- u2^(-eta-1)*(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val_2 <- u2^(-eta-1)*(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta-1)

    C_val_A <-(u1_A^(-eta)+u2_A^(-eta)-1)^(-1/eta)

  }


  # Use Copula functions to write each block of likelihood function
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 0), C_val_1, 1)
  term1 <- log(abs(term1))

  term2 <- C_val_1 - C_val_2
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- (c_u2_val_1 - c_u2_val_2) * (-u2_t2)
  term3 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 1), term3, 1)
  term3[term3 < 0] <- 1
  term3 <- log(abs(term3))

  term4 <- (c_u2_val_1) * (-u2_t2)
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 1), term4, 1)
  term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  # # left truncation
  # term5 <- ifelse((indata1[, "A1"] == 0), u2_A, C_val_A) # A1 = 0 means S(A1, A2) = S(0, A2) = S2(A2)
  # term5 <- ifelse((indata1[, "A1"] == 0) & (indata2[, "A2"] == 0), 1, term5) # A1 = A2 = 0 means S(0, 0) = 1
  # term5 <- log(abs(term5))


  # logL<-(-1)*sum( term1 + term2 + term3 + term4 - term5)

  logL<-(-1)*sum( term1 + term2 + term3 + term4)
  return(logL)
}


# general LT, sieve
ic_scmprisk_copula_log_lik_sieve_pseudo_LT_copula2 <- function(p, fitted, x1, x2, t1_left, t1_right, t2, indata1, indata2, bl1, br1, b2, b2_d, b1_A, b2_A, m1, m2, r1, r2, p1, p2, quantiles = NULL, copula)
{

  if (copula != "Copula2") {

    eta <- (p[1]) # anti-log
    phi1 <- p[2:(1+1+m1)]
    beta1 <- p[(1+1+m1+1):length(p)]
    ep1 <- cumsum(exp(phi1))


    phi2 <- fitted[1:(1+m2)]
    beta2 <- fitted[(1+m2+1):length(fitted)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  if (copula == "Copula2") {

    alpha <- exp(p[1])/(1+exp(p[1])) # anti-log
    kappa <- exp(p[2]) # anti-log
    phi1 <- p[3:(2+1+m1)]
    beta1 <- p[(2+1+m1+1):length(p)]
    ep1 <- cumsum(exp(phi1))

    phi2 <- fitted[1:(1+m2)]
    beta2 <- fitted[(1+m2+1):length(fitted)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  # survival probability
  gLl1<-G(exp(x1%*%beta1)*(bl1%*%ep1),r1) # AD, left end
  gLr1<-G(exp(x1%*%beta1)*(br1%*%ep1),r1)
  gL2 <- G(exp(x2%*%beta2)*(b2%*%ep2),r2)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2 = exp(-gL2)

  # derivative of survival (= negative density)
  u2_t2 <- u2 * (-1) * pG(exp(x2%*%beta2)*(b2%*%ep2),r2) * (b2_d %*% ep2) * exp(x2 %*% beta2)

  # Left truncation
  gL1_A <- G(exp(x1%*%beta1)*(b1_A%*%ep1),r1)
  u1_A = exp(-gL1_A)
  gL2_A <- G(exp(x2%*%beta2)*(b2_A%*%ep2),r2)
  u2_A = exp(-gL2_A)

  if (copula == "Copula2") {

    # Copula2 Copula function for joint distribution probability
    C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_2 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    c_u2_val_1 <- C_val_1^(1+1/kappa) * (1 + ((u1_left^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2
    c_u2_val_2 <- C_val_2^(1+1/kappa) * (1 + ((u1_right^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2

    C_val_A <- (1 + ((u1_A^(-1/kappa)-1)^(1/alpha) + (u2_A^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)

  }



  if (copula == "Joe") {

    # Joe joint and density functions
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2)^(eta) - (1-u1_left)^(eta)*(1-u2)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_right)^(eta) + (1-u2)^(eta) - (1-u1_right)^(eta)*(1-u2)^(eta) )^(1/eta)
    c_u2_val_1 <- ((1-u1_left)^eta + (1-u2)^eta - ((1-u1_left)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_left)^eta)
    c_u2_val_2 <- ((1-u1_right)^eta + (1-u2)^eta - ((1-u1_right)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_right)^eta)

    C_val_A <- 1 - ( (1-u1_A)^(eta) + (1-u2_A)^(eta) - (1-u1_A)^(eta)*(1-u2_A)^(eta) )^(1/eta)

  }


  if (copula == "Gumbel") {

    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_right))^eta + (-log(u2))^eta)^(1/eta))
    c_u2_val_1 <- gh_F(u1_left,u2,eta)*((-log(u1_left))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    c_u2_val_2 <- gh_F(u1_right,u2,eta)*((-log(u1_right))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2

    C_val_A <- exp(-((-log(u1_A))^eta + (-log(u2_A))^eta)^(1/eta))

  }


  if (copula == "AMH") {

    # AHM joint and desity functions
    C_val_1 <- u1_left * u2 /(1 - eta * (1-u1_left) * (1-u2))
    C_val_2 <- u1_right * u2 /(1 - eta * (1-u1_right) * (1-u2))
    c_u2_val_1 <- u1_left*(1-eta*(1-u1_left))/(1-eta*(1-u1_left)*(1-u2))^2
    c_u2_val_2 <- u1_right*(1-eta*(1-u1_right))/(1-eta*(1-u1_right)*(1-u2))^2

    C_val_A <- u1_A * u2_A /(1 - eta * (1-u1_A) * (1-u2_A))

  }


  if (copula == "Frank") {

    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    c_u2_val_1 <- (exp((-1) * eta*u1_left)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    c_u2_val_2 <- (exp((-1) * eta*u1_right)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))

    C_val_A <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_A)-1)*(exp((-1) * eta*u2_A)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    C_val_1 <-(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta)
    c_u2_val_1 <- u2^(-eta-1)*(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val_2 <- u2^(-eta-1)*(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta-1)

    C_val_A <-(u1_A^(-eta)+u2_A^(-eta)-1)^(-1/eta)

  }


  # Use Copula functions to write each block of likelihood function
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 0), C_val_1, 1)
  term1 <- log(abs(term1))

  term2 <- C_val_1 - C_val_2
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- (c_u2_val_1 - c_u2_val_2) * (-u2_t2)
  term3 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 1), term3, 1)
  term3[term3 < 0] <- 1
  term3 <- log(abs(term3))

  term4 <- (c_u2_val_1) * (-u2_t2)
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 1), term4, 1)
  term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  # left truncation
  term5 <- ifelse((indata1[, "A1"] == 0), u2_A, C_val_A) # A1 = 0 means S(A1, A2) = S(0, A2) = S2(A2)
  term5 <- ifelse((indata1[, "A1"] == 0) & (indata2[, "A2"] == 0), 1, term5) # A1 = A2 = 0 means S(0, 0) = 1
  term5 <- log(abs(term5))


  logL<-(-1)*sum( term1 + term2 + term3 + term4 - term5)
  return(logL)
}

# scmprisk, ic, semi-parametric margins, left truncation, general, omitting LT
ic_scmprisk_copula_log_lik_sieve_LT_1_copula2 <-function(p, x1, x2, t1_left, t1_right, t2, indata1, indata2, bl1, br1, b2, b2_d, b1_A, b2_A, m1, m2, r1, r2, p1, p2, quantiles = NULL, copula)
{

  if (copula != "Copula2") {

    eta <- (p[1])
    phi1 <- p[2:(1+1+m1)]
    beta1 <- p[(1+1+m1+1):(1+1+m1+p1)]
    ep1 <- cumsum(exp(phi1))

    phi2 <- p[(1+1+m1+p1+1):(1+1+m1+p1+1+m2)]
    beta2 <- p[(1+1+m1+p1+1+m2+1):length(p)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  if (copula == "Copula2") {

    alpha <- exp(p[1])/(1+exp(p[1])) # anti-log
    kappa <- exp(p[2]) # anti-log
    phi1 <- p[3:(1+1+1+m1)]
    beta1 <- p[(1+1+1+m1+1):(1+1+1+m1+p1)]
    ep1 <- cumsum(exp(phi1))

    phi2 <- p[(1+1+1+m1+p1+1):(1+1+1+m1+p1+1+m2)]
    beta2 <- p[(1+1+1+m1+p1+1+m2+1):length(p)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  # survival probability
  gLl1<-G(exp(x1%*%beta1)*(bl1%*%ep1),r1) # AD, left end
  gLr1<-G(exp(x1%*%beta1)*(br1%*%ep1),r1)
  gL2 <- G(exp(x2%*%beta2)*(b2%*%ep2),r2)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2 = exp(-gL2)

  # derivative of survival (= negative density)
  u2_t2 <- u2 * (-1) * pG(exp(x2%*%beta2)*(b2%*%ep2),r2) * (b2_d %*% ep2) * exp(x2 %*% beta2)

  # Left truncation
  gL1_A <- G(exp(x1%*%beta1)*(b1_A%*%ep1),r1)
  u1_A = exp(-gL1_A)
  gL2_A <- G(exp(x2%*%beta2)*(b2_A%*%ep2),r2)
  u2_A = exp(-gL2_A)

  if (copula == "Copula2") {

    # Copula2 Copula function for joint distribution probability
    C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_2 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    c_u2_val_1 <- C_val_1^(1+1/kappa) * (1 + ((u1_left^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2
    c_u2_val_2 <- C_val_2^(1+1/kappa) * (1 + ((u1_right^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2

    C_val_A <- (1 + ((u1_A^(-1/kappa)-1)^(1/alpha) + (u2_A^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)

  }



  if (copula == "Joe") {

    # Joe joint and density functions
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2)^(eta) - (1-u1_left)^(eta)*(1-u2)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_right)^(eta) + (1-u2)^(eta) - (1-u1_right)^(eta)*(1-u2)^(eta) )^(1/eta)
    c_u2_val_1 <- ((1-u1_left)^eta + (1-u2)^eta - ((1-u1_left)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_left)^eta)
    c_u2_val_2 <- ((1-u1_right)^eta + (1-u2)^eta - ((1-u1_right)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_right)^eta)

    C_val_A <- 1 - ( (1-u1_A)^(eta) + (1-u2_A)^(eta) - (1-u1_A)^(eta)*(1-u2_A)^(eta) )^(1/eta)

  }


  if (copula == "Gumbel") {

    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_right))^eta + (-log(u2))^eta)^(1/eta))
    c_u2_val_1 <- gh_F(u1_left,u2,eta)*((-log(u1_left))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    c_u2_val_2 <- gh_F(u1_right,u2,eta)*((-log(u1_right))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2

    C_val_A <- exp(-((-log(u1_A))^eta + (-log(u2_A))^eta)^(1/eta))

  }


  if (copula == "AMH") {

    # AHM joint and desity functions
    C_val_1 <- u1_left * u2 /(1 - eta * (1-u1_left) * (1-u2))
    C_val_2 <- u1_right * u2 /(1 - eta * (1-u1_right) * (1-u2))
    c_u2_val_1 <- u1_left*(1-eta*(1-u1_left))/(1-eta*(1-u1_left)*(1-u2))^2
    c_u2_val_2 <- u1_right*(1-eta*(1-u1_right))/(1-eta*(1-u1_right)*(1-u2))^2

    C_val_A <- u1_A * u2_A /(1 - eta * (1-u1_A) * (1-u2_A))

  }


  if (copula == "Frank") {

    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    c_u2_val_1 <- (exp((-1) * eta*u1_left)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    c_u2_val_2 <- (exp((-1) * eta*u1_right)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))

    C_val_A <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_A)-1)*(exp((-1) * eta*u2_A)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    C_val_1 <-(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta)
    c_u2_val_1 <- u2^(-eta-1)*(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val_2 <- u2^(-eta-1)*(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta-1)

    C_val_A <-(u1_A^(-eta)+u2_A^(-eta)-1)^(-1/eta)

  }

  # Use Copula functions to write each block of likelihood function
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 0), C_val_1, 1)
  term1 <- log(abs(term1))

  term2 <- C_val_1 - C_val_2
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- (c_u2_val_1 - c_u2_val_2) * (-u2_t2)
  term3 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 1), term3, 1)
  term3[term3 < 0] <- 1
  term3 <- log(abs(term3))

  term4 <- (c_u2_val_1) * (-u2_t2)
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 1), term4, 1)
  term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  # # left truncation
  # term5 <- ifelse((indata1[, "A1"] == 0), u2_A, C_val_A)
  # term5 <- ifelse((indata1[, "A1"] == 0) & (indata2[, "A2"] == 0), 1, term5) # A1 = A2 = 0 means S(0, 0) = 1
  # term5 <- log(abs(term5))

  # logL<-(-1)*sum( term1 + term2 + term3 + term4 - term5)

  logL<-(-1)*sum( term1 + term2 + term3 + term4)
  return(logL)
}

# scmprisk, ic, semi-parametric margins, left truncation, general, for Copula2 only
ic_scmprisk_copula_log_lik_sieve_LT_copula2 <-function(p, x1, x2, t1_left, t1_right, t2, indata1, indata2, bl1, br1, b2, b2_d, b1_A, b2_A, m1, m2, r1, r2, p1, p2, quantiles = NULL, copula)
{

  if (copula != "Copula2") {

    eta <- exp(p[1])
    phi1 <- p[2:(1+1+m1)]
    beta1 <- p[(1+1+m1+1):(1+1+m1+p1)]
    ep1 <- cumsum(exp(phi1))

    phi2 <- p[(1+1+m1+p1+1):(1+1+m1+p1+1+m2)]
    beta2 <- p[(1+1+m1+p1+1+m2+1):length(p)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  if (copula == "Copula2") {

    alpha <- exp(p[1])/(1+exp(p[1]))
    kappa <- exp(p[2])
    phi1 <- p[3:(1+1+1+m1)]
    beta1 <- p[(1+1+1+m1+1):(1+1+1+m1+p1)]
    ep1 <- cumsum(exp(phi1))

    phi2 <- p[(1+1+1+m1+p1+1):(1+1+1+m1+p1+1+m2)]
    beta2 <- p[(1+1+1+m1+p1+1+m2+1):length(p)] # coefficients
    ep2 <- cumsum(exp(phi2))

  }

  # survival probability
  gLl1<-G(exp(x1%*%beta1)*(bl1%*%ep1),r1) # AD, left end
  gLr1<-G(exp(x1%*%beta1)*(br1%*%ep1),r1)
  gL2 <- G(exp(x2%*%beta2)*(b2%*%ep2),r2)

  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2 = exp(-gL2)

  # derivative of survival (= negative density)
  u2_t2 <- u2 * (-1) * pG(exp(x2%*%beta2)*(b2%*%ep2),r2) * (b2_d %*% ep2) * exp(x2 %*% beta2)

  # Left truncation
  gL1_A <- G(exp(x1%*%beta1)*(b1_A%*%ep1),r1)
  u1_A = exp(-gL1_A)
  gL2_A <- G(exp(x2%*%beta2)*(b2_A%*%ep2),r2)
  u2_A = exp(-gL2_A)

  if (copula == "Copula2") {

    # Copula2 Copula function for joint distribution probability
    C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    C_val_2 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
    c_u2_val_1 <- C_val_1^(1+1/kappa) * (1 + ((u1_left^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2
    c_u2_val_2 <- C_val_2^(1+1/kappa) * (1 + ((u1_right^(-1/kappa)-1)/(u2^(-1/kappa)-1))^(1/alpha) )^(alpha-1)  * u2^(-1/kappa-1) # wrt to u2

    C_val_A <- (1 + ((u1_A^(-1/kappa)-1)^(1/alpha) + (u2_A^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)

  }



  if (copula == "Joe") {

    # Joe joint and density functions
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2)^(eta) - (1-u1_left)^(eta)*(1-u2)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_right)^(eta) + (1-u2)^(eta) - (1-u1_right)^(eta)*(1-u2)^(eta) )^(1/eta)
    c_u2_val_1 <- ((1-u1_left)^eta + (1-u2)^eta - ((1-u1_left)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_left)^eta)
    c_u2_val_2 <- ((1-u1_right)^eta + (1-u2)^eta - ((1-u1_right)^eta)*((1-u2)^eta) )^(1/eta - 1)  *  ((1-u2)^(eta-1) - (1-u2)^(eta-1)*(1-u1_right)^eta)

    C_val_A <- 1 - ( (1-u1_A)^(eta) + (1-u2_A)^(eta) - (1-u1_A)^(eta)*(1-u2_A)^(eta) )^(1/eta)

  }


  if (copula == "Gumbel") {

    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_right))^eta + (-log(u2))^eta)^(1/eta))
    c_u2_val_1 <- gh_F(u1_left,u2,eta)*((-log(u1_left))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    c_u2_val_2 <- gh_F(u1_right,u2,eta)*((-log(u1_right))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2

    C_val_A <- exp(-((-log(u1_A))^eta + (-log(u2_A))^eta)^(1/eta))

  }


  if (copula == "AMH") {

    # AHM joint and desity functions
    C_val_1 <- u1_left * u2 /(1 - eta * (1-u1_left) * (1-u2))
    C_val_2 <- u1_right * u2 /(1 - eta * (1-u1_right) * (1-u2))
    c_u2_val_1 <- u1_left*(1-eta*(1-u1_left))/(1-eta*(1-u1_left)*(1-u2))^2
    c_u2_val_2 <- u1_right*(1-eta*(1-u1_right))/(1-eta*(1-u1_right)*(1-u2))^2

    C_val_A <- u1_A * u2_A /(1 - eta * (1-u1_A) * (1-u2_A))

  }


  if (copula == "Frank") {

    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))
    c_u2_val_1 <- (exp((-1) * eta*u1_left)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))
    c_u2_val_2 <- (exp((-1) * eta*u1_right)-1)*exp((-1) * eta*u2)/((1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2)-1)/(exp((-1) * eta)-1))*(exp((-1) * eta)-1))

    C_val_A <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_A)-1)*(exp((-1) * eta*u2_A)-1)/(exp((-1) * eta)-1))

  }


  if (copula == "Clayton") {

    C_val_1 <-(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta)
    c_u2_val_1 <- u2^(-eta-1)*(u1_left^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val_2 <- u2^(-eta-1)*(u1_right^(-eta)+u2^(-eta)-1)^(-1/eta-1)

    C_val_A <-(u1_A^(-eta)+u2_A^(-eta)-1)^(-1/eta)

  }

  # Use Copula functions to write each block of likelihood function
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 0), C_val_1, 1)
  term1 <- log(abs(term1))

  term2 <- C_val_1 - C_val_2
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- (c_u2_val_1 - c_u2_val_2) * (-u2_t2)
  term3 <- ifelse((indata1[,"status"] == 1) & (indata2[,"statusD"] == 1), term3, 1)
  term3[term3 < 0] <- 1
  term3 <- log(abs(term3))

  term4 <- (c_u2_val_1) * (-u2_t2)
  term4 <- ifelse((indata1[,"status"] == 0) & (indata2[,"statusD"] == 1), term4, 1)
  term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  # left truncation
  term5 <- ifelse((indata1[, "A1"] == 0), u2_A, C_val_A)
  term5 <- ifelse((indata1[, "A1"] == 0) & (indata2[, "A2"] == 0), 1, term5) # A1 = A2 = 0 means S(0, 0) = 1
  term5 <- log(abs(term5))

  logL<-(-1)*sum( term1 + term2 + term3 + term4 - term5)
  return(logL)
}



