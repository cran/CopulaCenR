# penalty derivatives

# penalty 1
dLASSO1 <- function(z,a) a
dSCAD1 <- function(z,alph,a) a*((z<=a)+ifelse(alph*a-z>0,alph*a-z,0)/(a*(alph-1))*(z>a))
dMCP1 <- function(z,tau2,a) (z<=tau2*a)*(a-z/tau2)
dSICA1 <- function(z,tau3,a) a*tau3*(tau3+1)/(z+tau3)^2
dSELO1 <- function(z,gam2,a) a*gam2/(log(2)*(2*z+gam2)*(z+gam2))

# penalty 2
dLASSO2 <- function(z,b) b
dSCAD2 <- function(z,alph,b) b*((z<=b)+ifelse(alph*b-z>0,alph*b-z,0)/(b*(alph-1))*(z>b))
dMCP2 <- function(z,tau2,b) (z<=tau2*b)*(b-z/tau2)
dSICA2 <- function(z,tau3,b) b*tau3*(tau3+1)/(z+tau3)^2
dSELO2 <- function(z,gam2,b) b*gam2/(log(2)*(2*z+gam2)*(z+gam2))

# rc, semi-competing, sieve
data_process_scmprisk_rc_sp <- function(indata, var_list, l, u, m)
{

  indata <- indata[,c('obs_time','status',var_list)]

  t <- indata[,'obs_time']

  dim_m <- dim(as.matrix(indata[,var_list]))
  n <- dim_m[1]

  # matix
  x <- as.matrix(indata[,var_list], nrow = n)
  p <- dim(x)[2]

  # BP Bernstein?????????????????????
  b <- matrix(0,nrow = n,ncol = m+1) # terminal

  for (i in 0:m) {
    b[,(i+1)] <- bern(i,m,l,u,t)
  }


  # BP derivatives???Bernstein????????????t???????????????t????????????????????????
  b_d  <- matrix(0,nrow = n,ncol = m+1) # terminal

  for (i in 0:m) {
    b_d[,(i+1)] <- bern_derivative(i,m,l,u,t)
  }


  return(list(indata=indata, t=t,
              n=n, p=p, x=x, var_list = var_list,
              b=b, b_d=b_d))
}


# # Berstein polynomials: j from 0-m, m for degree, l/u for range of time, t for specific times
# bern <- function(j,m,l,u,t){
#
#   if (j > m) {
#     b = 0
#   } else {
#     b = (factorial(m)/(factorial(j)*factorial(m-j)))*(((t-l)/(u-l))^j)*((1-(t-l)/(u-l))^(m-j))
#     # Berstein polynomials?????????
#   }
#   return(b)
# }


# # Berstein polynomials derivatives
# bern_derivative <- function(j,m,l,u,t){
#
#   if (j != 0) {
#     # b = m*bern(j-1, m-1, l, u, t) - m*bern(j, m-1, l, u, t)
#
#     b = (m/(u - l)) * (bern(j-1, m-1, l, u, t)-bern(j, m-1, l, u, t))
#   } else if (j == 0) {
#     b = -1 * m * (u - t)^(m - 1)/((u - l)^m)
#   }
#   return(b)
# }



# rc, scmprisk, cox margin 2

scmprisk_log_lik_baseline2 <- function(par,fitted,x2,indata2,b2,b2_d,m2, quantiles = NULL)
{

  phi2 <- par
  ep2 <- cumsum(exp(phi2))

  beta2 <- fitted[1:length(fitted)]

  t2<-indata2[,"obs_time"]

  # baseline
  Lambda2 <- b2%*%ep2
  lambda2 <- b2_d%*%ep2

  # survival and density
  u2 <- exp((-1)*Lambda2*exp(x2%*%beta2))
  f2 <- u2 * exp(x2%*%beta2) * lambda2


  term1 <- ifelse(indata2[,'status']==1, f2, 1)
  term1 <- log(abs(term1))

  term2 <- ifelse(indata2[,'status']==0, u2, 1)
  term2 <- log(abs(term2))

  logL <- (-1)*sum(term1 + term2)
  return(logL)

}


# rc, scmprisk, cox margin 1

scmprisk_log_lik_baseline1 <- function(par,fitted,x1,indata1,b1,b1_d,m1, quantiles = NULL)
{

  phi1 <- par
  ep1 <- cumsum(exp(phi1))

  beta1 <- fitted[1:length(fitted)]

  t1<-indata1[,"obs_time"]

  # baseline
  Lambda1 <- b1%*%ep1
  lambda1 <- b1_d%*%ep1

  # survival and density
  u1 <- exp((-1)*Lambda1*exp(x1%*%beta1))
  f1 <- u1 * exp(x1%*%beta1) * lambda1


  term1 <- ifelse(indata1[,'status']==1, f1, 1)
  term1 <- log(abs(term1))

  term2 <- ifelse(indata1[,'status']==0, u1, 1)
  term2 <- log(abs(term2))

  logL <- (-1)*sum(term1 + term2)
  return(logL)

}


# rc, scmprisk, copula, baseline
rc_scmprisk_copula_loglik_sieve_baseline <- function(par,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula)
{

  phi1 <- par[1:(m1+1)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- par[(m1+2):(m2+m1+2)]
  ep2 <- cumsum(exp(phi2))

  eta <- exp(fitted[1])
  beta1 <- fitted[2:(1+dim(x1)[2])]
  beta2 <- fitted[(2+dim(x1)[2]):length(fitted)]

  t1<-indata1[,"obs_time"]
  t2<-indata2[,"obs_time"]

  # survival and density
  u1 <- exp((-1)*(b1%*%ep1)*exp(x1%*%beta1))
  u2 <- exp((-1)*(b2%*%ep2)*exp(x2%*%beta2))
  f1 <- u1 * exp(x1%*%beta1) * (b1_d%*%ep1)
  f2 <- u2 * exp(x2%*%beta2) * (b2_d%*%ep2)

  # copula function
  if (copula == "Gumbel") {

    c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
    c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

  }


  if (copula == "Clayton") {

    # Clayton Copula function for joint distribution probability
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta) # C(u,v)

  }


  C_val <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), C_val, 1)
  term1 <- log(abs(C_val))

  term2 <- c_u1_val * f1
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- c_u2_val * f2
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
  term3 <- log(abs(term3))

  term4 <- c_val * f1 * f2
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
  # term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  if (all(is.finite(term1), is.finite(term2), is.finite(term3), is.finite(term4))) {
    logL <- (-1) * sum(term1 + term2 + term3 + term4)
  } else {
    logL <- 0
  }

  return(logL)

}



# rc_scmprisk_copula_loglik_sieve_baseline2 <- function(par,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula)
# {
#
#   phi2 <- par
#   ep2 <- cumsum(exp(phi2))
#
#   eta <- exp(fitted[1])
#   phi1 <- fitted[2:(m1+2)]
#   ep1 <- cumsum(exp(phi1))
#   beta1 <- fitted[(m1+3):(m1+2+dim(x1)[2])]
#   beta2 <- fitted[(m1+3+dim(x1)[2]):length(fitted)]
#
#   t1<-indata1[,"obs_time"]
#   t2<-indata2[,"obs_time"]
#
#   # survival and density
#   u1 <- exp((-1)*(b1%*%ep1)*exp(x1%*%beta1))
#   u2 <- exp((-1)*(b2%*%ep2)*exp(x2%*%beta2))
#   f1 <- u1 * exp(x1%*%beta1) * (b1_d%*%ep1)
#   f2 <- u2 * exp(x2%*%beta2) * (b2_d%*%ep2)
#
#   # copula function
#   if (copula == "Gumbel") {
#
#     c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
#     c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
#     c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
#     C_val<-gh_F(u1,u2,eta)
#
#   }
#
#
#   if (copula == "Clayton") {
#
#     # Clayton Copula function for joint distribution probability
#     c_val<-clt_f(u1,u2,eta)
#     c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta) # C(u,v)
#
#   }
#
#
#   C_val <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), C_val, 1)
#   term1 <- log(abs(C_val))
#
#   term2 <- c_u1_val * f1
#   term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
#   term2 <- log(abs(term2))
#
#   term3 <- c_u2_val * f2
#   term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
#   term3 <- log(abs(term3))
#
#   term4 <- c_val * f1 * f2
#   term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
#   # term4[term4 < 0] <- 1
#   term4 <- log(abs(term4))
#
#   logL<-(-1)*sum( term1 + term2 + term3 + term4 )
#   return(logL)
#
# }



# rc_scmprisk_copula_loglik_sieve_baseline1 <- function(par,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula)
# {
#
#   phi1 <- par
#   ep1 <- cumsum(exp(phi1))
#
#   eta <- exp(fitted[1])
#   phi2 <- fitted[2:(m2+2)]
#   ep2 <- cumsum(exp(phi2))
#
#   beta1 <- fitted[(m2+3):(m2+2+dim(x1)[2])]
#   beta2 <- fitted[(m2+3+dim(x1)[2]):length(fitted)]
#
#   t1<-indata1[,"obs_time"]
#   t2<-indata2[,"obs_time"]
#
#   # survival and density
#   u1 <- exp((-1)*(b1%*%ep1)*exp(x1%*%beta1))
#   u2 <- exp((-1)*(b2%*%ep2)*exp(x2%*%beta2))
#   f1 <- u1 * exp(x1%*%beta1) * (b1_d%*%ep1)
#   f2 <- u2 * exp(x2%*%beta2) * (b2_d%*%ep2)
#
#   # copula function
#   if (copula == "Gumbel") {
#
#     c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
#     c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
#     c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
#     C_val<-gh_F(u1,u2,eta)
#
#   }
#
#
#   if (copula == "Clayton") {
#
#     # Clayton Copula function for joint distribution probability
#     c_val<-clt_f(u1,u2,eta)
#     c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta) # C(u,v)
#
#   }
#
#
#   C_val <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), C_val, 1)
#   term1 <- log(abs(C_val))
#
#   term2 <- c_u1_val * f1
#   term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
#   term2 <- log(abs(term2))
#
#   term3 <- c_u2_val * f2
#   term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
#   term3 <- log(abs(term3))
#
#   term4 <- c_val * f1 * f2
#   term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
#   # term4[term4 < 0] <- 1
#   term4 <- log(abs(term4))
#
#   logL<-(-1)*sum( term1 + term2 + term3 + term4 )
#   return(logL)
#
# }

# d_logL_beta1_group <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,kk)
# {
#   beta1 <- p
#
#   eta <- exp(fitted[1])
#   phi1 <- fitted[2:(m1+2)]
#   ep1 <- cumsum(exp(phi1))
#   phi2 <- fitted[(m1+3):(m1+3+m2)]
#   ep2 <- cumsum(exp(phi2))
#   beta2 <- fitted[(m1+4+m2):(m1+3+m2+dim(x2)[2])]
#
#   t1<-indata1[,"obs_time"]
#   t2<-indata2[,"obs_time"]
#
#   # baseline
#   Lambda1 <- b1%*%ep1
#   Lambda2 <- b2%*%ep2
#   lambda1 <- b1_d%*%ep1
#   lambda2 <- b2_d%*%ep2
#
#   # survival and density
#   u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
#   u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
#   f1 <- u1 * lambda1 * exp(x1 %*% beta1)
#   f2 <- u2 * lambda2 * exp(x2 %*% beta2)
#
#   x1_kk <- array(data = t(x1[,kk]), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#   # ?????????
#   # ?????????????????????????????????1*1?????????????????????????????????(1,1,n)???????????????(k,1,n)???????????????????????????????????????(k,1,n)
#   u1_beta1_part1 <- array(data = (-1) * Lambda1 * u1 * exp(x1%*%beta1), dim = c(1,1,n))
#   u1_beta1 <- array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%u1_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#
#   f1_beta1_part1 <- array(data = exp(x1%*%beta1) * lambda1, dim = c(1,1,n))
#   f1_beta1_part2 <- array(data = u1 * exp(x1%*%beta1) * lambda1, dim = c(1,1,n))
#   f1_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f1_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%f1_beta1_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   if (copula == "Clayton") {
#
#     C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
#     C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
#
#     C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
#     c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
#     C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))
#
#     C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
#     c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
#     C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))
#
#   }
#
#   if (copula == "Gumbel") {
#     # ??????????????????????????????????????????
#     A <- (-log(u1))^(eta-1)/u1
#     B <- (-log(u2))^(eta-1)/u2
#     D <- (-log(u1))^eta + (-log(u2))^eta
#     E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
#     A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
#     A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))
#
#     c_val <- E*(D^(1/eta)+eta-1)
#     C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
#     C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
#     C_val<-gh_F(u1,u2,eta)
#
#     E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)
#
#     C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
#     c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)
#
#     E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)
#
#     c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
#     C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
#     # --------------------------------------
#     B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
#     B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
#     E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)
#
#     C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
#     c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)
#
#     E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)
#
#     c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
#     C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2
#
#   }
#
#   # ?????????
#   f_val_beta1_part1 <- array(data = f2*c_1*f1, dim = c(1,1,n))
#   f_val_beta1_part2 <- array(data = f2*c_val, dim = c(1,1,n))
#   f_val_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f_val_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1[,,x], function(z){z%*%f_val_beta1_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   C_val_beta1_part1 <- array(data = C_1, dim = c(1,1,n))
#   C_val_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_val_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   C_2_beta1_part1 <- array(data = c_val, dim = c(1,1,n))
#   C_2_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_2_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_1_beta1_part1 <- array(data = C_1_1, dim = c(1,1,n))
#   C_1_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_1_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   term <- array(data = c(0), dim = c(length(kk),1,n))
#
#   section1 <- array(data = 1/C_val, dim = c(1,1,n))
#   term1 <- array(data = sapply(1:n, function(x){sapply(C_val_beta1[,,x], function(z){z%*%section1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#   index_00 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 0))  # ???which?????????????????????????????????????????????
#   term[,,index_00] <- term1[,,index_00]
#
#   section2_1 <- array(data = 1/C_1, dim = c(1,1,n))
#   section2_2 <- array(data = 1/f1, dim = c(1,1,n))
#   term2 <- array(data = sapply(1:n, function(x){sapply(C_1_beta1[,,x], function(z){z%*%section2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1[,,x], function(z){z%*%section2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) # ??????(k,1,n)
#   index_10 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 0))
#   term[,,index_10] <- term2[,,index_10]
#
#   section3 <- array(data = 1/C_2, dim = c(1,1,n))
#   term3 <- array(data = sapply(1:n, function(x){sapply(C_2_beta1[,,x], function(z){z%*%section3[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#   index_01 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 1))
#   term[,,index_01] <- term3[,,index_01]
#
#   section4 <- array(data = 1/(c_val*f1*f2), dim = c(1,1,n))
#   term4 <- array(data = sapply(1:n, function(x){sapply(f_val_beta1[,,x], function(z){z%*%section4[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#   index_11 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 1))
#   term[,,index_11] <- term4[,,index_11]
#
#   term[which(is.na(term))] <- 0
#   logL_beta1_first <- apply(term, c(1,2), sum) # ?????????????????????
#
#   return(logL_beta1_first)
#
# }


# dd_logL_beta1_group <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,kk)
# {
#
#   beta1 <- p
#
#   eta <- exp(fitted[1])
#   phi1 <- fitted[2:(m1+2)]
#   ep1 <- cumsum(exp(phi1))
#   phi2 <- fitted[(m1+3):(m1+3+m2)]
#   ep2 <- cumsum(exp(phi2))
#   beta2 <- fitted[(m1+4+m2):(m1+3+m2+dim(x2)[2])]
#
#   t1<-indata1[,"obs_time"]
#   t2<-indata2[,"obs_time"]
#
#   # baseline
#   Lambda1 <- b1%*%ep1
#   Lambda2 <- b2%*%ep2
#   lambda1 <- b1_d%*%ep1
#   lambda2 <- b2_d%*%ep2
#
#   # survival and density
#   u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
#   u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
#   f1 <- u1 * lambda1 * exp(x1 %*% beta1)
#   f2 <- u2 * lambda2 * exp(x2 %*% beta2)
#
#   x1_kk <- array(data = t(x1[,kk]), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#   # ?????????
#   u1_beta1_part1 <- array(data = (-1) * Lambda1 * u1 * exp(x1%*%beta1), dim = c(1,1,n))
#   u1_beta1 <- array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%u1_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#
#   f1_beta1_part1 <- array(data = exp(x1%*%beta1) * lambda1, dim = c(1,1,n))
#   f1_beta1_part2 <- array(data = u1 * exp(x1%*%beta1) * lambda1, dim = c(1,1,n))
#   f1_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f1_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%f1_beta1_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   # ?????????
#   u1_beta1_2_part1_1 <- array(data = (-1) * Lambda1 * exp(x1%*%beta1), dim = c(1,1,n))
#   u1_beta1_2_part2_1 <- array(data = (-1) * Lambda1 * u1 * exp(x1%*%beta1), dim = c(1,1,n))
#   u1_beta1_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%u1_beta1_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   u1_beta1_2_part2_2 <- array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%u1_beta1_2_part2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   u1_beta1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%u1_beta1_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%u1_beta1_2_part2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n))
#
#
#   f1_beta1_2_part1 <- array(data = lambda1 * exp(x1%*%beta1), dim = c(1,1,n))
#   f1_beta1_2_part2_1 <- array(data = 2 * lambda1 * exp(x1%*%beta1), dim = c(1,1,n))
#   f1_beta1_2_part3_1 <- array(data = lambda1 * u1 * exp(x1%*%beta1), dim = c(1,1,n))
#   f1_beta1_2_part2_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f1_beta1_2_part2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f1_beta1_2_part3_2 <- array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%f1_beta1_2_part3_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f1_beta1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1_2[,,x], function(z){z%*%f1_beta1_2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%f1_beta1_2_part2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(x1_kk[,,x], function(z){z%*%f1_beta1_2_part3_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n))
#
#   if (copula == "Clayton") {
#
#     C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
#     C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
#
#     C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
#     c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
#     C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))
#
#     C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
#     c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
#     C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))
#
#   }
#
#   if (copula == "Gumbel") {
#     # ??????????????????????????????????????????
#     A <- (-log(u1))^(eta-1)/u1
#     B <- (-log(u2))^(eta-1)/u2
#     D <- (-log(u1))^eta + (-log(u2))^eta
#     E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
#     A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
#     A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))
#
#     c_val <- E*(D^(1/eta)+eta-1)
#     C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
#     C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
#     C_val<-gh_F(u1,u2,eta)
#
#     E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)
#
#     C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
#     c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)
#
#     E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)
#
#     c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
#     C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
#     # --------------------------------------
#     B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
#     B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
#     E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)
#
#     C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
#     c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)
#
#     E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)
#
#     c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
#     C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2
#
#   }
#
#   # ?????????
#   f_val_beta1_part1 <- array(data = f2*c_1*f1, dim = c(1,1,n))
#   f_val_beta1_part2 <- array(data = f2*c_val, dim = c(1,1,n))
#   f_val_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f_val_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1[,,x], function(z){z%*%f_val_beta1_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   C_val_beta1_part1 <- array(data = C_1, dim = c(1,1,n))
#   C_val_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_val_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   C_2_beta1_part1 <- array(data = c_val, dim = c(1,1,n))
#   C_2_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_2_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_1_beta1_part1 <- array(data = C_1_1, dim = c(1,1,n))
#   C_1_beta1 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_1_beta1_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   # ?????????
#   f_val_beta1_2_part1_1 <- array(data = f2 * c_1_1 * f1, dim = c(1,1,n))
#   f_val_beta1_2_part2 <- array(data = f2 * c_1 * f1, dim = c(1,1,n))
#   f_val_beta1_2_part3_1 <- array(data = 2 * f2 * c_1, dim = c(1,1,n))
#   f_val_beta1_2_part4 <- array(data = f2 * c_val, dim = c(1,1,n))
#   f_val_beta1_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f_val_beta1_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f_val_beta1_2_part3_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f_val_beta1_2_part3_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f_val_beta1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%f_val_beta1_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u1_beta1_2[,,x], function(z){z%*%f_val_beta1_2_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1[,,x], function(z){z%*%f_val_beta1_2_part3_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1_2[,,x], function(z){z%*%f_val_beta1_2_part4[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n))
#
#
#   C_val_beta1_2_part1_1<- array(data = C_1_1, dim = c(1,1,n))
#   C_val_beta1_2_part2 <- array(data = C_1, dim = c(1,1,n))
#   C_val_beta1_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_val_beta1_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   C_val_beta1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_val_beta1_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u1_beta1_2[,,x], function(z){z%*%C_val_beta1_2_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n))
#
#
#   C_1_beta1_2_part1_1 <- array(data = C_1_1_1, dim = c(1,1,n))
#   C_1_beta1_2_part2 <- array(data = C_1_1, dim = c(1,1,n))
#   C_1_beta1_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_1_beta1_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   C_1_beta1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_1_beta1_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u1_beta1_2[,,x], function(z){z%*%C_1_beta1_2_part2[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#
#   C_2_beta1_2_part1_1 <- array(data = c_1, dim = c(1,1,n))
#   C_2_beta1_2_part2 <- array(data = c_val, dim = c(1,1,n))
#   C_2_beta1_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_2_beta1_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   C_2_beta1_2 <- array(data = sapply(1:n, function(x){sapply(u1_beta1[,,x], function(z){z%*%C_2_beta1_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u1_beta1_2[,,x], function(z){z%*%C_2_beta1_2_part2[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#
#   term <- array(data = c(0), dim = c(length(kk),length(kk),n))
#
#   section1_1 <- array(data = -1/C_val^2, dim = c(1,1,n))
#   section1_3 <- array(data = 1/C_val, dim = c(1,1,n))
#   section1_2 <- array(data = sapply(1:n, function(x){sapply(C_val_beta1[,,x], function(z){z%*%section1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   term1 <- array(data = sapply(1:n, function(x){sapply(C_val_beta1[,,x], function(z){z%*%section1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(C_val_beta1_2[,,x], function(z){z%*%section1_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_00 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 0))
#   term[,,index_00] <- term1[,,index_00]
#
#   section2_1 <- array(data = -1/C_1^2, dim = c(1,1,n))
#   section2_2 <- array(data = sapply(1:n, function(x){sapply(C_1_beta1[,,x], function(z){z%*%section2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section2_3 <- array(data = 1/C_1, dim = c(1,1,n))
#   section2_4 <- array(data = -1/f1^2, dim = c(1,1,n))
#   section2_5 <- array(data = sapply(1:n, function(x){sapply(f1_beta1[,,x], function(z){z%*%section2_4[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section2_6 <- array(data = 1/f1, dim = c(1,1,n))
#   term2 <- array(data = sapply(1:n, function(x){sapply(C_1_beta1[,,x], function(z){z%*%section2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(C_1_beta1_2[,,x], function(z){z%*%section2_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1[,,x], function(z){z%*%section2_5[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f1_beta1_2[,,x], function(z){z%*%section2_6[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_10 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 0))
#   term[,,index_10] <- term2[,,index_10]
#
#   section3_1 <- array(data = -1/C_2^2, dim = c(1,1,n))
#   section3_2 <- array(data = sapply(1:n, function(x){sapply(C_2_beta1[,,x], function(z){z%*%section3_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section3_3 <- array(data = 1/C_2, dim = c(1,1,n))
#   term3 <- array(data = sapply(1:n, function(x){sapply(C_2_beta1[,,x], function(z){z%*%section3_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(C_2_beta1_2[,,x], function(z){z%*%section3_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_01 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 1))
#   term[,,index_01] <- term3[,,index_01]
#
#   section4_1 <- array(data = -1/(c_val*f1*f2)^2, dim = c(1,1,n))
#   section4_2 <- array(data = sapply(1:n, function(x){sapply(f_val_beta1[,,x], function(z){z%*%section4_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section4_3 <- array(data = 1/(c_val*f1*f2), dim = c(1,1,n))
#   term4 <- array(data = sapply(1:n, function(x){sapply(f_val_beta1[,,x], function(z){z%*%section4_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f_val_beta1_2[,,x], function(z){z%*%section4_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_11 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 1))
#   term[,,index_11] <- term4[,,index_11]
#
#   term[which(is.na(term))] <- 0
#   logL_beta1_second <- apply(term, c(1,2), sum)
#
#   return(logL_beta1_second)
#
# }

# rc, scmprisk, copula, first order derivative on beta1
d_logL_beta1 <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,dd)
{
  beta1 <- p

  eta <- exp(fitted[1])
  phi1 <- fitted[2:(m1+2)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- fitted[(m1+3):(m1+3+m2)]
  ep2 <- cumsum(exp(phi2))
  beta2 <- fitted[(m1+4+m2):(m1+3+m2+dim(x2)[2])]

  t1<-indata1[,"obs_time"]
  t2<-indata2[,"obs_time"]

  # baseline
  Lambda1 <- b1%*%ep1
  Lambda2 <- b2%*%ep2
  lambda1 <- b1_d%*%ep1
  lambda2 <- b2_d%*%ep2

  # survival and density
  u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
  u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
  f1 <- u1 * lambda1 * exp(x1 %*% beta1)
  f2 <- u2 * lambda2 * exp(x2 %*% beta2)

  x1_dd <- x1[,dd]
  # ?????????
  u1_beta1 <- function(Lambda1,x1,beta1){
    # u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
    result <- (-1) * Lambda1 * u1 * exp(x1%*%beta1) * x1_dd
    return(result)
  }

  f1_beta1 <- function(Lambda1,lambda1,x1,beta1){
    # u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
    result <- lambda1 * (u1_beta1(Lambda1,x1,beta1)*exp(x1%*%beta1) + u1*exp(x1%*%beta1)*x1_dd)
    return(result)
  }


  if (copula == "Clayton") {

    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)

    C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
    c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
    C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))

    C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
    c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
    C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))

  }

  if (copula == "Gumbel") {
    # ??????????????????????????????????????????
    A <- (-log(u1))^(eta-1)/u1
    B <- (-log(u2))^(eta-1)/u2
    D <- (-log(u1))^eta + (-log(u2))^eta
    E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
    A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
    A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))

    c_val <- E*(D^(1/eta)+eta-1)
    C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

    E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)

    C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
    c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)

    E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)

    c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
    C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
    # --------------------------------------
    B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
    B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
    E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)

    C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
    c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)

    E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)

    c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
    C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2

  }

  # first order derivative
  f_val_beta1 <- function(Lambda1,lambda1,x1,beta1){
    result <- f2*(c_1*u1_beta1(Lambda1,x1,beta1)*f1+c_val*f1_beta1(Lambda1,lambda1,x1,beta1))
    return(result)
  }

  C_val_beta1 <- function(Lambda1,x1,beta1){
    result <- C_1*u1_beta1(Lambda1,x1,beta1)
    return(result)
  }

  C_1_beta1 <- function(Lambda1,x1,beta1){
    result <- C_1_1*u1_beta1(Lambda1,x1,beta1)
    return(result)
  }

  C_2_beta1 <- function(Lambda1,x1,beta1){
    result <- c_val*u1_beta1(Lambda1,x1,beta1)
    return(result)
  }


  term1 <- (1/C_val)*C_val_beta1(Lambda1,x1,beta1)
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term1, 0)
  term1[which(is.na(term1))] <- 0

  term2 <- (1/C_1)*C_1_beta1(Lambda1,x1,beta1)+(1/f1)*f1_beta1(Lambda1,lambda1,x1,beta1)
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 0)
  term2[which(is.na(term2))] <- 0

  term3 <- (1/C_2)*C_2_beta1(Lambda1,x1,beta1)
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 0)
  term3[which(is.na(term3))] <- 0

  term4 <- (1/(c_val*f1*f2))*f_val_beta1(Lambda1,lambda1,x1,beta1)
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 0)
  term4[which(is.na(term4))] <- 0

  logL_beta1_first <- sum( term1 + term2 + term3 + term4 )

  return(logL_beta1_first)

}

# rc, scmprisk, copula, second order derivative on beta1
dd_logL_beta1 <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,dd)
{

  beta1 <- p

  eta <- exp(fitted[1])
  phi1 <- fitted[2:(m1+2)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- fitted[(m1+3):(m1+3+m2)]
  ep2 <- cumsum(exp(phi2))
  beta2 <- fitted[(m1+4+m2):(m1+3+m2+dim(x2)[2])]

  t1<-indata1[,"obs_time"]
  t2<-indata2[,"obs_time"]

  # baseline
  Lambda1 <- b1%*%ep1
  Lambda2 <- b2%*%ep2
  lambda1 <- b1_d%*%ep1
  lambda2 <- b2_d%*%ep2

  # survival and density
  u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
  u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
  f1 <- u1 * lambda1 * exp(x1 %*% beta1)
  f2 <- u2 * lambda2 * exp(x2 %*% beta2)

  x1_dd <- x1[,dd]
  # first order derivative
  u1_beta1 <- function(Lambda1,x1,beta1){
    # u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
    result <- (-1) * Lambda1 * u1 * exp(x1%*%beta1) * x1_dd
    return(result)
  }

  f1_beta1 <- function(Lambda1,lambda1,x1,beta1){
    # u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
    result <- lambda1 * (u1_beta1(Lambda1,x1,beta1)*exp(x1%*%beta1) + u1*exp(x1%*%beta1)*x1_dd)
    return(result)
  }

  # second order derivative
  u1_beta1_2 <- function(Lambda1,x1,beta1){
    # u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
    result <- (-1)*Lambda1*x1_dd*(u1_beta1(Lambda1,x1,beta1)*exp(x1%*%beta1) + u1*exp(x1%*%beta1)*x1_dd)
    return(result)
  }

  f1_beta1_2 <- function(Lambda1,lambda1,x1,beta1){
    # u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
    result <- lambda1*(u1_beta1_2(Lambda1,x1,beta1)*exp(x1%*%beta1)+2*u1_beta1(Lambda1,x1,beta1)*exp(x1%*%beta1)*x1_dd+u1*exp(x1%*%beta1)*(x1_dd^2))
    return(result)
  }

  if (copula == "Clayton") {

    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)

    C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
    c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
    C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))

    C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
    c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
    C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))

  }

  if (copula == "Gumbel") {
    # ??????????????????????????????????????????
    A <- (-log(u1))^(eta-1)/u1
    B <- (-log(u2))^(eta-1)/u2
    D <- (-log(u1))^eta + (-log(u2))^eta
    E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
    A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
    A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))

    c_val <- E*(D^(1/eta)+eta-1)
    C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

    E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)

    C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
    c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)

    E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)

    c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
    C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
    # --------------------------------------
    B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
    B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
    E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)

    C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
    c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)

    E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)

    c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
    C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2

  }

  # first order derivative
  f_val_beta1 <- function(Lambda1,lambda1,x1,beta1){
    result <- f2*(c_1*u1_beta1(Lambda1,x1,beta1)*f1+c_val*f1_beta1(Lambda1,lambda1,x1,beta1))
    return(result)
  }

  C_val_beta1 <- function(Lambda1,x1,beta1){
    result <- C_1*u1_beta1(Lambda1,x1,beta1)
    return(result)
  }

  C_1_beta1 <- function(Lambda1,x1,beta1){
    result <- C_1_1*u1_beta1(Lambda1,x1,beta1)
    return(result)
  }

  C_2_beta1 <- function(Lambda1,x1,beta1){
    result <- c_val*u1_beta1(Lambda1,x1,beta1)
    return(result)
  }

  # second order derivative
  f_val_beta1_2 <- function(Lambda1,lambda1,x1,beta1){
    result <- f2*(c_1_1*(u1_beta1(Lambda1,x1,beta1)^2)*f1+c_1*u1_beta1_2(Lambda1,x1,beta1)*f1+2*c_1*u1_beta1(Lambda1,x1,beta1)*f1_beta1(Lambda1,lambda1,x1,beta1)+c_val*f1_beta1_2(Lambda1,lambda1,x1,beta1))
    return(result)
  }

  C_val_beta1_2 <- function(Lambda1,x1,beta1){
    result <- C_1_1*(u1_beta1(Lambda1,x1,beta1)^2)+C_1*u1_beta1_2(Lambda1,x1,beta1)
    return(result)
  }

  C_1_beta1_2 <- function(Lambda1,x1,beta1){
    result <- C_1_1_1*(u1_beta1(Lambda1,x1,beta1)^2)+C_1_1*u1_beta1_2(Lambda1,x1,beta1)
    return(result)
  }

  C_2_beta1_2 <- function(Lambda1,x1,beta1){
    result <- c_1*(u1_beta1(Lambda1,x1,beta1)^2)+c_val*u1_beta1_2(Lambda1,x1,beta1)
    return(result)
  }


  term1 <- (-1/C_val^2)*C_val_beta1(Lambda1,x1,beta1)^2+(1/C_val)*C_val_beta1_2(Lambda1,x1,beta1)
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term1, 0)
  term1[which(is.na(term1))] <- 0

  term2 <- (-1/C_1^2)*(C_1_beta1(Lambda1,x1,beta1)^2)+(1/C_1)*C_1_beta1_2(Lambda1,x1,beta1)+(-1/f1^2)*(f1_beta1(Lambda1,lambda1,x1,beta1)^2)+(1/f1)*f1_beta1_2(Lambda1,lambda1,x1,beta1)
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 0)
  term2[which(is.na(term2))] <- 0

  term3 <- (-1/C_2^2)*(C_2_beta1(Lambda1,x1,beta1)^2)+(1/C_2)*C_2_beta1_2(Lambda1,x1,beta1)
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 0)
  term3[which(is.na(term3))] <- 0

  term4 <- (-1/(c_val*f1*f2)^2)*(f_val_beta1(Lambda1,lambda1,x1,beta1)^2)+(c_val*f1*f2)^(-1)*f_val_beta1_2(Lambda1,lambda1,x1,beta1)
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 0)
  term4[which(is.na(term4))] <- 0

  logL_beta1_second <- sum( term1 + term2 + term3 + term4 )

  return(logL_beta1_second)

}


# d_logL_beta2_group <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,kk)
# {
#   beta2 <- p
#
#   eta <- exp(fitted[1])
#   phi1 <- fitted[2:(m1+2)]
#   ep1 <- cumsum(exp(phi1))
#   phi2 <- fitted[(m1+3):(m1+3+m2)]
#   ep2 <- cumsum(exp(phi2))
#   beta1 <- fitted[(m1+4+m2):(m1+3+m2+dim(x1)[2])]
#
#   t1<-indata1[,"obs_time"]
#   t2<-indata2[,"obs_time"]
#
#   # baseline
#   Lambda1 <- b1%*%ep1
#   Lambda2 <- b2%*%ep2
#   lambda1 <- b1_d%*%ep1
#   lambda2 <- b2_d%*%ep2
#
#   # survival and density
#   u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
#   u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
#   f1 <- u1 * lambda1 * exp(x1 %*% beta1)
#   f2 <- u2 * lambda2 * exp(x2 %*% beta2)
#
#   x2_kk <- array(data = t(x2[,kk]), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#   # ?????????
#   # ?????????????????????????????????1*1?????????????????????????????????(1,1,n)???????????????(k,1,n)???????????????????????????????????????(k,1,n)
#   u2_beta2_part1 <- array(data = (-1) * Lambda2 * u2 * exp(x2%*%beta2), dim = c(1,1,n))
#   u2_beta2 <- array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%u2_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#
#   f2_beta2_part1 <- array(data = exp(x2%*%beta2) * lambda2, dim = c(1,1,n))
#   f2_beta2_part2 <- array(data = u2 * exp(x2%*%beta2) * lambda2, dim = c(1,1,n))
#   f2_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f2_beta2_part1[,,x]})}, simplify = 'array') + sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%f2_beta2_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   if (copula == "Clayton") {
#
#     C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
#     C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
#
#     C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
#     c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
#     C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))
#
#     C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
#     c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
#     C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))
#
#   }
#
#   if (copula == "Gumbel") {
#     # ??????????????????????????????????????????
#     A <- (-log(u1))^(eta-1)/u1
#     B <- (-log(u2))^(eta-1)/u2
#     D <- (-log(u1))^eta + (-log(u2))^eta
#     E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
#     A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
#     A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))
#
#     c_val <- E*(D^(1/eta)+eta-1)
#     C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
#     C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
#     C_val<-gh_F(u1,u2,eta)
#
#     E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)
#
#     C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
#     c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)
#
#     E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)
#
#     c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
#     C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
#     # --------------------------------------
#     B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
#     B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
#     E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)
#
#     C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
#     c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)
#
#     E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)
#
#     c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
#     C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2
#
#   }
#
#   # ?????????
#   f_val_beta2_part1 <- array(data = f1*c_2*f2, dim = c(1,1,n))
#   f_val_beta2_part2 <- array(data = f1*c_val, dim = c(1,1,n))
#   f_val_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f_val_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2[,,x], function(z){z%*%f_val_beta2_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_val_beta2_part1 <- array(data = C_2, dim = c(1,1,n))
#   C_val_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_val_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_1_beta2_part1 <- array(data = c_val, dim = c(1,1,n))
#   C_1_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_1_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_2_beta2_part1 <- array(data = C_2_2, dim = c(1,1,n))
#   C_2_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_2_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   term <- array(data = c(0), dim = c(length(kk),1,n))
#
#   section1 <- array(data = 1/C_val, dim = c(1,1,n))
#   term1 <- array(data = sapply(1:n, function(x){sapply(C_val_beta2[,,x], function(z){z%*%section1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#   index_00 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 0))  # ???which?????????????????????????????????????????????
#   term[,,index_00] <- term1[,,index_00]
#
#
#   section2_1 <- array(data = 1/C_2, dim = c(1,1,n))
#   section2_2 <- array(data = 1/f2, dim = c(1,1,n))
#   term2 <- array(data = sapply(1:n, function(x){sapply(C_2_beta2[,,x], function(z){z%*%section2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2[,,x], function(z){z%*%section2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) # ??????(k,1,n)
#   index_01 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 1))
#   term[,,index_01] <- term2[,,index_01]
#
#   section3 <- array(data = 1/C_1, dim = c(1,1,n))
#   term3 <- array(data = sapply(1:n, function(x){sapply(C_1_beta2[,,x], function(z){z%*%section3[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#   index_10 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 0))
#   term[,,index_10] <- term3[,,index_10]
#
#   section4 <- array(data = 1/(c_val*f1*f2), dim = c(1,1,n))
#   term4 <- array(data = sapply(1:n, function(x){sapply(f_val_beta2[,,x], function(z){z%*%section4[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#   index_11 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 1))
#   term[,,index_11] <- term4[,,index_11]
#
#   term[which(is.na(term))] <- 0
#   logL_beta2_first <- apply(term, c(1,2), sum) # ?????????????????????
#
#   return(logL_beta2_first)
#
# }


# dd_logL_beta2_group <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,kk)
# {
#
#   beta2 <- p
#
#   eta <- exp(fitted[1])
#   phi1 <- fitted[2:(m1+2)]
#   ep1 <- cumsum(exp(phi1))
#   phi2 <- fitted[(m1+3):(m1+3+m2)]
#   ep2 <- cumsum(exp(phi2))
#   beta1 <- fitted[(m1+4+m2):(m1+3+m2+dim(x1)[2])]
#
#   t1<-indata1[,"obs_time"]
#   t2<-indata2[,"obs_time"]
#
#   # baseline
#   Lambda1 <- b1%*%ep1
#   Lambda2 <- b2%*%ep2
#   lambda1 <- b1_d%*%ep1
#   lambda2 <- b2_d%*%ep2
#
#   # survival and density
#   u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
#   u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
#   f1 <- u1 * lambda1 * exp(x1 %*% beta1)
#   f2 <- u2 * lambda2 * exp(x2 %*% beta2)
#
#   x2_kk <- array(data = t(x2[,kk]), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#   # ?????????
#   # ?????????????????????????????????1*1?????????????????????????????????(1,1,n)???????????????(k,1,n)???????????????????????????????????????(k,1,n)
#   u2_beta2_part1 <- array(data = (-1) * Lambda2 * u2 * exp(x2%*%beta2), dim = c(1,1,n))
#   u2_beta2 <- array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%u2_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) # ??????(k,1,n)
#
#
#   f2_beta2_part1 <- array(data = exp(x2%*%beta2) * lambda2, dim = c(1,1,n))
#   f2_beta2_part2 <- array(data = u2 * exp(x2%*%beta2) * lambda2, dim = c(1,1,n))
#   f2_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f2_beta2_part1[,,x]})}, simplify = 'array') + sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%f2_beta2_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#   # ?????????
#   u2_beta2_2_part1_1 <- array(data = (-1) * Lambda2 * exp(x2%*%beta2), dim = c(1,1,n))
#   u2_beta2_2_part2_1 <- array(data = (-1) * Lambda2 * u2 * exp(x2%*%beta2), dim = c(1,1,n))
#   u2_beta2_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%u2_beta2_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   u2_beta2_2_part2_2 <- array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%u2_beta2_2_part2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   u2_beta2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%u2_beta2_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%u2_beta2_2_part2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n))
#
#
#   f2_beta2_2_part1 <- array(data = lambda2 * exp(x2%*%beta2), dim = c(1,1,n))
#   f2_beta2_2_part2_1 <- array(data = 2 * lambda2 * exp(x2%*%beta2), dim = c(1,1,n))
#   f2_beta2_2_part3_1 <- array(data = lambda2 * u2 * exp(x2%*%beta2), dim = c(1,1,n))
#   f2_beta2_2_part2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f2_beta2_2_part2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f2_beta2_2_part3_2 <- array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%f2_beta2_2_part3_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f2_beta2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2_2[,,x], function(z){z%*%f2_beta2_2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%f2_beta2_2_part2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(x2_kk[,,x], function(z){z%*%f2_beta2_2_part3_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n))
#
#   if (copula == "Clayton") {
#
#     C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
#     C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
#     c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
#
#     C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
#     c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
#     C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))
#
#     C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
#     c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
#     c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
#     C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))
#
#   }
#
#   if (copula == "Gumbel") {
#     # ??????????????????????????????????????????
#     A <- (-log(u1))^(eta-1)/u1
#     B <- (-log(u2))^(eta-1)/u2
#     D <- (-log(u1))^eta + (-log(u2))^eta
#     E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
#     A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
#     A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))
#
#     c_val <- E*(D^(1/eta)+eta-1)
#     C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
#     C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
#     C_val<-gh_F(u1,u2,eta)
#
#     E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)
#
#     C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
#     c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)
#
#     E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)
#
#     c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
#     C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
#     # --------------------------------------
#     B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
#     B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
#     E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)
#
#     C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
#     c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)
#
#     E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)
#
#     c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
#     C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2
#
#   }
#
#   # ?????????
#   f_val_beta2_part1 <- array(data = f1*c_2*f2, dim = c(1,1,n))
#   f_val_beta2_part2 <- array(data = f1*c_val, dim = c(1,1,n))
#   f_val_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f_val_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2[,,x], function(z){z%*%f_val_beta2_part2[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_val_beta2_part1 <- array(data = C_2, dim = c(1,1,n))
#   C_val_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_val_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_1_beta2_part1 <- array(data = c_val, dim = c(1,1,n))
#   C_1_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_1_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   C_2_beta2_part1 <- array(data = C_2_2, dim = c(1,1,n))
#   C_2_beta2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_2_beta2_part1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))  # ??????(k,1,n)
#
#
#   # ?????????
#   f_val_beta2_2_part1_1 <- array(data = f1 * c_2_2 * f2, dim = c(1,1,n))
#   f_val_beta2_2_part2 <- array(data = f1 * c_2 * f2, dim = c(1,1,n))
#   f_val_beta2_2_part3_1 <- array(data = 2 * f1 * c_2, dim = c(1,1,n))
#   f_val_beta2_2_part4 <- array(data = f1 * c_val, dim = c(1,1,n))
#   f_val_beta2_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f_val_beta2_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f_val_beta2_2_part3_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f_val_beta2_2_part3_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   f_val_beta2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%f_val_beta2_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u2_beta2_2[,,x], function(z){z%*%f_val_beta2_2_part2[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2[,,x], function(z){z%*%f_val_beta2_2_part3_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2_2[,,x], function(z){z%*%f_val_beta2_2_part4[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#
#   C_val_beta2_2_part1_1 <- array(data = C_2_2, dim = c(1,1,n))
#   C_val_beta2_2_part2 <- array(data = C_2, dim = c(1,1,n))
#   C_val_beta2_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_val_beta2_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   C_val_beta2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_val_beta2_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u2_beta2_2[,,x], function(z){z%*%C_val_beta2_2_part2[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#
#   C_2_beta2_2_part1_1 <- array(data = C_2_2_2, dim = c(1,1,n))
#   C_2_beta2_2_part2 <- array(data = C_2_2, dim = c(1,1,n))
#   C_2_beta2_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_2_beta2_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   C_2_beta2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_2_beta2_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u2_beta2_2[,,x], function(z){z%*%C_2_beta2_2_part2[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#
#   C_1_beta2_2_part1_1 <- array(data = c_2, dim = c(1,1,n))
#   C_1_beta2_2_part2 <- array(data = c_val, dim = c(1,1,n))
#   C_1_beta2_2_part1_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_1_beta2_2_part1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   C_1_beta2_2 <- array(data = sapply(1:n, function(x){sapply(u2_beta2[,,x], function(z){z%*%C_1_beta2_2_part1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(u2_beta2_2[,,x], function(z){z%*%C_1_beta2_2_part2[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   term <- array(data = c(0), dim = c(length(kk),length(kk),n))
#
#   section1_1 <- array(data = -1/C_val^2, dim = c(1,1,n))
#   section1_3 <- array(data = 1/C_val, dim = c(1,1,n))
#   section1_2 <- array(data = sapply(1:n, function(x){sapply(C_val_beta2[,,x], function(z){z%*%section1_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   term1 <- array(data = sapply(1:n, function(x){sapply(C_val_beta2[,,x], function(z){z%*%section1_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(C_val_beta2_2[,,x], function(z){z%*%section1_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_00 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 0))
#   term[,,index_00] <- term1[,,index_00]
#
#   section2_1 <- array(data = -1/C_2^2, dim = c(1,1,n))
#   section2_2 <- array(data = sapply(1:n, function(x){sapply(C_2_beta2[,,x], function(z){z%*%section2_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section2_3 <- array(data = 1/C_2, dim = c(1,1,n))
#   section2_4 <- array(data = -1/f2^2, dim = c(1,1,n))
#   section2_5 <- array(data = sapply(1:n, function(x){sapply(f2_beta2[,,x], function(z){z%*%section2_4[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section2_6 <- array(data = 1/f2, dim = c(1,1,n))
#   term2 <- array(data = sapply(1:n, function(x){sapply(C_2_beta2[,,x], function(z){z%*%section2_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(C_2_beta2_2[,,x], function(z){z%*%section2_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2[,,x], function(z){z%*%section2_5[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f2_beta2_2[,,x], function(z){z%*%section2_6[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_01 <- which((indata1[,"status"] == 0) & (indata2[,"status"] == 1))
#   term[,,index_01] <- term2[,,index_01]
#
#   section3_1 <- array(data = -1/C_1^2, dim = c(1,1,n))
#   section3_2 <- array(data = sapply(1:n, function(x){sapply(C_1_beta2[,,x], function(z){z%*%section3_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section3_3 <- array(data = 1/C_1, dim = c(1,1,n))
#   term3 <- array(data = sapply(1:n, function(x){sapply(C_1_beta2[,,x], function(z){z%*%section3_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(C_1_beta2_2[,,x], function(z){z%*%section3_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_10 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 0))
#   term[,,index_10] <- term3[,,index_10]
#
#   section4_1 <- array(data = -1/(c_val*f1*f2)^2, dim = c(1,1,n))
#   section4_2 <- array(data = sapply(1:n, function(x){sapply(f_val_beta2[,,x], function(z){z%*%section4_1[,,x]})}, simplify = 'array'), dim = c(length(kk),1,n))
#   section4_3 <- array(data = 1/(c_val*f1*f2), dim = c(1,1,n))
#   term4 <- array(data = sapply(1:n, function(x){sapply(f_val_beta2[,,x], function(z){z%*%section4_2[,,x]})}, simplify = 'array'), dim = c(length(kk),length(kk),n)) + array(data = sapply(1:n, function(x){sapply(f_val_beta2_2[,,x], function(z){z%*%section4_3[,,x]})}, simplify = 'array'),dim = c(length(kk),length(kk),n))
#
#   index_11 <- which((indata1[,"status"] == 1) & (indata2[,"status"] == 1))
#   term[,,index_11] <- term4[,,index_11]
#
#   term[which(is.na(term))] <- 0
#   logL_beta2_second <- apply(term, c(1,2), sum)
#
#   return(logL_beta2_second)
#
# }

# rc, scmprisk, copula, first order derivative on beta2
d_logL_beta2 <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,dd)
{

  beta2 <- p

  eta <- exp(fitted[1])
  phi1 <- fitted[2:(m1+2)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- fitted[(m1+3):(m1+3+m2)]
  ep2 <- cumsum(exp(phi2))
  beta1 <- fitted[(m1+4+m2):(m1+3+m2+dim(x1)[2])]

  t1<-indata1[,"obs_time"]
  t2<-indata2[,"obs_time"]

  # baseline
  Lambda1 <- b1%*%ep1
  Lambda2 <- b2%*%ep2
  lambda1 <- b1_d%*%ep1
  lambda2 <- b2_d%*%ep2

  # survival and density
  u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
  u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
  f1 <- u1 * lambda1 * exp(x1 %*% beta1)
  f2 <- u2 * lambda2 * exp(x2 %*% beta2)

  x2_dd <- x2[,dd]

  u2_beta2 <- function(Lambda2,x2,beta2){
    # u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
    result <- (-1) * Lambda2 * u2 * exp(x2%*%beta2) * x2_dd
    return(result)
  }

  f2_beta2 <- function(Lambda2,lambda2,x2,beta2){
    # u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
    result <- lambda2 * (u2_beta2(Lambda2,x2,beta2)*exp(x2%*%beta2) + u2*exp(x2%*%beta2)*x2_dd)
    return(result)
  }


  if (copula == "Clayton") {

    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)

    C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
    c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
    C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))

    C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
    c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
    C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))

  }

  if (copula == "Gumbel") {
    # ??????????????????????????????????????????
    A <- (-log(u1))^(eta-1)/u1
    B <- (-log(u2))^(eta-1)/u2
    D <- (-log(u1))^eta + (-log(u2))^eta
    E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
    A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
    A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))

    c_val <- E*(D^(1/eta)+eta-1)
    C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

    E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)

    C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
    c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)

    E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)

    c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
    C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
    # --------------------------------------
    B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
    B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
    E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)

    C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
    c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)

    E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)

    c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
    C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2

  }

  # first order derivative
  f_val_beta2 <- function(Lambda2,lambda2,x2,beta2){
    result <- f1*(c_2*u2_beta2(Lambda2,x2,beta2)*f2+c_val*f2_beta2(Lambda2,lambda2,x2,beta2))
    return(result)
  }

  C_val_beta2 <- function(Lambda2,x2,beta2){
    result <- C_2*u2_beta2(Lambda2,x2,beta2)
    return(result)
  }

  C_1_beta2 <- function(Lambda2,x2,beta2){
    result <- c_val*u2_beta2(Lambda2,x2,beta2)
    return(result)
  }

  C_2_beta2 <- function(Lambda2,x2,beta2){
    result <- C_2_2*u2_beta2(Lambda2,x2,beta2)
    return(result)
  }


  term1 <- (1/C_val)*C_val_beta2(Lambda2,x2,beta2)
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term1, 0)
  term1[which(is.na(term1))] <- 0

  term2 <- (1/C_2)*C_2_beta2(Lambda2,x2,beta2)+(1/f2)*f2_beta2(Lambda2,lambda2,x2,beta2)
  term2 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term2, 0)
  term2[which(is.na(term2))] <- 0

  term3 <- (1/C_1)*C_1_beta2(Lambda2,x2,beta2)
  term3 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term3, 0)
  term3[which(is.na(term3))] <- 0

  term4 <- (1/(c_val*f1*f2))*f_val_beta2(Lambda2,lambda2,x2,beta2)
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 0)
  term4[which(is.na(term4))] <- 0

  logL_beta2_first <- sum( term1 + term2 + term3 + term4 )

  return(logL_beta2_first)

}


# rc, scmprisk, copula, second order derivative on beta2
dd_logL_beta2 <- function(p,fitted,x1,x2,indata1,indata2,b1,b2,b1_d,b2_d,m1,m2, quantiles = NULL, copula,dd)
{

  beta2 <- p

  eta <- exp(fitted[1])
  phi1 <- fitted[2:(m1+2)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- fitted[(m1+3):(m1+3+m2)]
  ep2 <- cumsum(exp(phi2))
  beta1 <- fitted[(m1+4+m2):(m1+3+m2+dim(x1)[2])]

  t1<-indata1[,"obs_time"]
  t2<-indata2[,"obs_time"]

  # baseline
  Lambda1 <- b1%*%ep1
  Lambda2 <- b2%*%ep2
  lambda1 <- b1_d%*%ep1
  lambda2 <- b2_d%*%ep2

  # survival and density
  u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
  u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
  f1 <- u1 * lambda1 * exp(x1 %*% beta1)
  f2 <- u2 * lambda2 * exp(x2 %*% beta2)

  x2_dd <- x2[,dd]
  # first order derivative
  u2_beta2 <- function(Lambda2,x2,beta2){
    # u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
    result <- (-1) * Lambda2 * u2 * exp(x2%*%beta2) * x2_dd
    return(result)
  }

  f2_beta2 <- function(Lambda2,lambda2,x2,beta2){
    # u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
    result <- lambda2 * (u2_beta2(Lambda2,x2,beta2)*exp(x2%*%beta2) + u2*exp(x2%*%beta2)*x2_dd)
    return(result)
  }

  # second order derivative
  u2_beta2_2 <- function(Lambda2,x2,beta2){
    # u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
    result <- (-1)*Lambda2*x2_dd*(u2_beta2(Lambda2,x2,beta2)*exp(x2%*%beta2) + u2*exp(x2%*%beta2)*x2_dd)
    return(result)
  }

  f2_beta2_2 <- function(Lambda2,lambda2,x2,beta2){
    # u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
    result <- lambda2*(u2_beta2_2(Lambda2,x2,beta2)*exp(x2%*%beta2)+2*u2_beta2(Lambda2,x2,beta2)*exp(x2%*%beta2)*x2_dd+u2*exp(x2%*%beta2)*(x2_dd^2))
    return(result)
  }

  if (copula == "Clayton") {

    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta)
    C_1<-u1^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_2<-u2^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_val<-(1+eta) * (u1*u2)^(-eta-1) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)

    C_1_1<-(1+eta) * (u1^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u1^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_1<-c_val*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)
    c_1_1<-c_1*(((2*eta+1)*u1^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u1)+c_val*(eta*(2*eta+1)*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u1^(-2))
    C_1_1_1<-(1+eta)*C_1_1*(u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u1))+(1+eta)*C_1*(eta*u1^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u1^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u1^2))

    C_2_2<-(1+eta) * (u2^(-eta-2)) * (u1^(-eta)+u2^(-eta)-1)^(-1/eta-1) * (u2^(-eta)/(u1^(-eta)+u2^(-eta)-1)-1)
    c_2<-c_val*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)
    c_2_2<-c_2*(((2*eta+1)*u2^(-eta-1))/(u1^(-eta)+u2^(-eta)-1)-(1+eta)/u2)+c_val*(eta*(2*eta+1)*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*(2*eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1+eta)*u2^(-2))
    C_2_2_2<-(1+eta)*C_2_2*(u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1)-(1/u2))+(1+eta)*C_2*(eta*u2^(-2*eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-2)-(eta+1)*u2^(-eta-2)*(u1^(-eta)+u2^(-eta)-1)^(-1)+(1/u2^2))

  }

  if (copula == "Gumbel") {
    # ??????????????????????????????????????????
    A <- (-log(u1))^(eta-1)/u1
    B <- (-log(u2))^(eta-1)/u2
    D <- (-log(u1))^eta + (-log(u2))^eta
    E <- gh_F(u1,u2,eta) * A*B*D^(1/eta-2)
    A_u1 <- (-1)*u1^(-2)*((-log(u1))^(eta-2))*(eta-1-log(u1))
    A_u1_u1 <- 2*u1^(-3)*((-log(u1))^(eta-2))*(eta-1-log(u1))-(eta-2)*u1^(-3)*((-log(u1))^(eta-3))*(eta-1-log(u1))-u1^(-3)*((-log(u1))^(eta-2))

    c_val <- E*(D^(1/eta)+eta-1)
    C_1<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    C_2<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

    E_u1 <- C_1*A*B*D^(1/eta-2)+C_val*A_u1*B*D^(1/eta-2)+(2*eta-1)*C_val*A^2*B*D^(1/eta-3)

    C_1_1 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*A^2+(eta-1)*A^2*(1/D)+A_u1)
    c_1 <- E_u1*(D^(1/eta)+eta-1)-E*A*D^(1/eta-1)

    E_u1_u1 <- C_1_1*A*B*D^(1/eta-2)+2*C_1*A_u1*B*D^(1/eta-2)+2*(2*eta-1)*C_1*A^2*B*D^(1/eta-3)+C_val*A_u1_u1*B*D^(1/eta-2)+3*(2*eta-1)*C_val*A_u1*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A^3*B*D^(1/eta-4)

    c_1_1 <- E_u1_u1*(D^(1/eta)+eta-1)-2*E_u1*A*D^(1/eta-1)-E*A_u1*D^(1/eta-1)+(1-eta)*E*A^2*D^(1/eta-2)
    C_1_1_1 <- C_1*A^2*D^(2/eta-2)+2*(eta-1)*C_val*A^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*A*A_u1+(eta-1)*C_1*A^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*A^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*A*A_u1+C_1*D^(1/eta-1)*A_u1+C_val*D^(1/eta-1)*A_u1_u1
    # --------------------------------------
    B_u2 <- (-1)*u2^(-2)*((-log(u2))^(eta-2))*(eta-1-log(u2))
    B_u2_u2 <- 2*u2^(-3)*((-log(u2))^(eta-2))*(eta-1-log(u2))-(eta-2)*u2^(-3)*((-log(u2))^(eta-3))*(eta-1-log(u2))-u2^(-3)*((-log(u2))^(eta-2))
    E_u2 <- C_2*A*B*D^(1/eta-2)+C_val*A*B_u2*D^(1/eta-2)+(2*eta-1)*C_val*A*B^2*D^(1/eta-3)

    C_2_2 <- C_val*D^(1/eta-1)*(D^(1/eta-1)*B^2+(eta-1)*B^2*(1/D)+B_u2)
    c_2 <- E_u2*(D^(1/eta)+eta-1)-E*B*D^(1/eta-1)

    E_u2_u2 <- C_2_2*A*B*D^(1/eta-2)+2*C_2*A*B_u2*D^(1/eta-2)+2*(2*eta-1)*C_2*A*B^2*D^(1/eta-3)+C_val*A*B_u2_u2*D^(1/eta-2)+3*(2*eta-1)*C_val*B_u2*A*B*D^(1/eta-3)+(2*eta-1)*(3*eta-1)*C_val*A*B^3*D^(1/eta-4)

    c_2_2 <- E_u2_u2*(D^(1/eta)+eta-1)-2*E_u2*B*D^(1/eta-1)-E*B_u2*D^(1/eta-1)+(1-eta)*E*B^2*D^(1/eta-2)
    C_2_2_2 <- C_2*B^2*D^(2/eta-2)+2*(eta-1)*C_val*B^3*D^(2/eta-3)+2*C_val*D^(2/eta-2)*B*B_u2+(eta-1)*C_2*B^2*D^(1/eta-2)+(eta-1)*(2*eta-1)*C_val*B^3*D^(1/eta-3)+3*(eta-1)*C_val*D^(1/eta-2)*B*B_u2+C_2*D^(1/eta-1)*B_u2+C_val*D^(1/eta-1)*B_u2_u2

  }

  # first order derivative
  f_val_beta2 <- function(Lambda2,lambda2,x2,beta2){
    result <- f1*(c_2*u2_beta2(Lambda2,x2,beta2)*f2+c_val*f2_beta2(Lambda2,lambda2,x2,beta2))
    return(result)
  }

  C_val_beta2 <- function(Lambda2,x2,beta2){
    result <- C_2*u2_beta2(Lambda2,x2,beta2)
    return(result)
  }

  C_1_beta2 <- function(Lambda2,x2,beta2){
    result <- c_val*u2_beta2(Lambda2,x2,beta2)
    return(result)
  }

  C_2_beta2 <- function(Lambda2,x2,beta2){
    result <- C_2_2*u2_beta2(Lambda2,x2,beta2)
    return(result)
  }

  # second order derivative
  f_val_beta2_2 <- function(Lambda2,lambda2,x2,beta2){
    result <- f1*(c_2_2*(u2_beta2(Lambda2,x2,beta2)^2)*f2+c_2*u2_beta2_2(Lambda2,x2,beta2)*f2+2*c_2*u2_beta2(Lambda2,x2,beta2)*f2_beta2(Lambda2,lambda2,x2,beta2)+c_val*f2_beta2_2(Lambda2,lambda2,x2,beta2))
    return(result)
  }

  C_val_beta2_2 <- function(Lambda2,x2,beta2){
    result <- C_2_2*(u2_beta2(Lambda2,x2,beta2)^2)+C_2*u2_beta2_2(Lambda2,x2,beta2)
    return(result)
  }

  C_2_beta2_2 <- function(Lambda2,x2,beta2){
    result <- C_2_2_2*(u2_beta2(Lambda2,x2,beta2)^2)+C_2_2*u2_beta2_2(Lambda2,x2,beta2)
    return(result)
  }

  C_1_beta2_2 <- function(Lambda2,x2,beta2){
    result <- c_2*(u2_beta2(Lambda2,x2,beta2)^2)+c_val*u2_beta2_2(Lambda2,x2,beta2)
    return(result)
  }


  term1 <- (-1/C_val^2) * (C_val_beta2(Lambda2,x2,beta2)^2) + (1/C_val)*C_val_beta2_2(Lambda2,x2,beta2)
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term1, 0)
  term1[which(is.na(term1))] <- 0

  term2 <- (-1/C_2^2)*C_2_beta2(Lambda2,x2,beta2)^2+(1/C_2)*C_2_beta2_2(Lambda2,x2,beta2)+(-1/f2^2)*f2_beta2(Lambda2,lambda2,x2,beta2)^2+(1/f2)*f2_beta2_2(Lambda2,lambda2,x2,beta2)
  term2 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term2, 0)
  term2[which(is.na(term2))] <- 0

  term3 <- (-1/C_1^2) * (C_1_beta2(Lambda2,x2,beta2)^2) + (1/C_1)*C_1_beta2_2(Lambda2,x2,beta2)
  term3 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term3, 0)
  term3[which(is.na(term3))] <- 0

  term4 <- (-1/(c_val*f1*f2)^2) * (f_val_beta2(Lambda2,lambda2,x2,beta2)^2) +(c_val*f1*f2)^(-1)*f_val_beta2_2(Lambda2,lambda2,x2,beta2)
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 0)
  term4[which(is.na(term4))] <- 0

  logL_beta2_second <- sum( term1 + term2 + term3 + term4 )

  return(logL_beta2_second)

}


# rc, scmprisk, copula, eta
rc_scmprisk_copula_loglik_sieve_eta <- function(par, fitted, x1, x2,indata1, indata2, b1, b2, b1_d, b2_d, m1, m2, quantiles = NULL, copula)
{

  eta <- exp(par) # anti-log

  phi1 <- fitted[1:(m1+1)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- fitted[(m1+2):(m1+2+m2)]
  ep2 <- cumsum(exp(phi2))
  beta1 <- fitted[(m1+3+m2):(m1+2+m2+dim(x1)[2])]
  beta2 <- fitted[(m1+3+m2+dim(x1)[2]):length(fitted)]

  t1<-indata1[,"obs_time"]
  t2<-indata2[,"obs_time"]

  # survival and density
  u1 <- exp((-1)*(b1%*%ep1)*exp(x1%*%beta1))
  u2 <- exp((-1)*(b2%*%ep2)*exp(x2%*%beta2))
  f1 <- u1 * exp(x1%*%beta1) * (b1_d%*%ep1)
  f2 <- u2 * exp(x2%*%beta2) * (b2_d%*%ep2)

  # copula function
  if (copula == "Gumbel") {

    c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
    c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

  }


  if (copula == "Clayton") {

    # Clayton Copula function for joint distribution probability
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta) # C(u,v)

  }

  term1 <- C_val
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term1, 1)
  term1 <- log(abs(C_val))

  term2 <- c_u1_val * f1
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- c_u2_val * f2
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
  term3 <- log(abs(term3))

  term4 <- c_val * f1 * f2
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
  # term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  logL<-(-1)*sum( term1 + term2 + term3 + term4 )

  return(logL)

}


# rc, scmprisk, copula, log-likelihood
rc_scmprisk_copula_loglik <- function(beta, fitted, x1, x2,indata1, indata2, b1, b2, b1_d, b2_d, m1, m2, quantiles = NULL, copula)
{

  eta <- exp(fitted[1])
  phi1 <- fitted[2:(m1+2)]
  ep1 <- cumsum(exp(phi1))
  phi2 <- fitted[(m1+3):length(fitted)]
  ep2 <- cumsum(exp(phi2))
  beta1 <- beta[1:dim(x1)[2]]
  beta2 <- beta[(1+dim(x1)[2]):length(beta)]

  t1<-indata1[, "obs_time"]
  t2<-indata2[, "obs_time"]

  # baseline
  Lambda1 <- b1%*%ep1
  Lambda2 <- b2%*%ep2
  lambda1 <- b1_d%*%ep1
  lambda2 <- b2_d%*%ep2

  # survival and density
  u1 <- exp(-Lambda1 * exp(x1 %*% beta1))
  u2 <- exp(-Lambda2 * exp(x2 %*% beta2))
  f1 <- u1 * lambda1 * exp(x1 %*% beta1)
  f2 <- u2 * lambda2 * exp(x2 %*% beta2)

  if (copula == "Gumbel") {

    c_val<-((-log(u1))^(eta-1)*(-log(u2))^(eta-1)*gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(2/eta-2))/(u1*u2)-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-2)*(1/(u1*u2))*(1/eta-1)*eta*(-log(u1))^(eta-1)*(-log(u2))^(eta-1)
    c_u1_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u1))^(eta-1)/u1
    c_u2_val<-gh_F(u1,u2,eta)*((-log(u1))^eta+(-log(u2))^eta)^(1/eta-1)*(-log(u2))^(eta-1)/u2
    C_val<-gh_F(u1,u2,eta)

  }


  if (copula == "Clayton") {

    # Clayton Copula function for joint distribution probability
    c_val<-clt_f(u1,u2,eta)
    c_u1_val<-u1^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    c_u2_val<-u2^(-eta-1)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-1)
    C_val<-(u1^(-eta)+u2^(-eta)-1)^(-1/eta) # C(u,v)

  }

  term1 <- C_val
  term1 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 0), term1, 1)
  term1 <- log(abs(term1))

  term2 <- c_u1_val * f1
  term2 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 0), term2, 1)
  term2 <- log(abs(term2))

  term3 <- c_u2_val * f2
  term3 <- ifelse((indata1[,"status"] == 0) & (indata2[,"status"] == 1), term3, 1)
  term3 <- log(abs(term3))

  term4 <- c_val * f1 * f2
  term4 <- ifelse((indata1[,"status"] == 1) & (indata2[,"status"] == 1), term4, 1)
  # term4[term4 < 0] <- 1
  term4 <- log(abs(term4))

  logL <- sum( term1 + term2 + term3 + term4 )

  return(logL)

}
