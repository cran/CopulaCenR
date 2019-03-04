
###### functions for building likelihoods #######

# Berstein polynomials: j from 0-m, m for degree, l/u for range of time, t for specific times
bern <- function(j,m,l,u,t){
  b = (factorial(m)/(factorial(j)*factorial(m-j)))*(((t-l)/(u-l))^j)*((1-(t-l)/(u-l))^(m-j))
  return(b)
}

# Density function of Clayton copula
clt_f<-function(u1,u2,eta)
{
  c_val<-(1+eta)*(u1*u2)^(-1-eta)*(u1^(-eta)+u2^(-eta)-1)^(-1/eta-2)
  return(c_val)
}

# 1st order partial derivative of Gumbel copula
gh_F<-function(u1,u2,eta)
{
  c_val<-exp(-((-log(u1))^eta+(-log(u2))^eta)^(1/eta))
  return(c_val)
}

# transformation function
G <- function(para,r) {

  if (r <= 2) { # Box-Cox transformations; if r = 1, then PH model
    result = (1/r)*((1+para)^r-1)
  }

  else { # Logarithmic transfomations; if r = 3, then PO model
    result = (1/(r-2))*log(1+para*(r-2))
  }
  return(result)
}


# derivative of transformation function
pG <- function(para,r) {

  if (r == 1) { # PH model
    result = rep(1,length(para))
  }

  if (r == 3) { # PO model
    result = 1/(1+para)
  }
  return(result)
}



# a function to calculate B matrix in the sandwich variance estimator: A^(-1) B A^(-1)
# for robust marginal sieve model
# Based on Zhou (2017) JASA
B_mar <- function(para,m,p,n,x1,x2,bl1,br1,bl2,br2,indata1, indata2,r){

  beta = para[1:p]
  phi = para[(p+1):(p+m+1)]
  ep = cumsum(exp(phi))

  Ll1 = exp(x1%*%beta)*(bl1%*%ep)
  Lr1 = exp(x1%*%beta)*(br1%*%ep)
  Ll2 = exp(x2%*%beta)*(bl2%*%ep)
  Lr2 = exp(x2%*%beta)*(br2%*%ep)

  gLl1<-G(exp(x1%*%beta)*(bl1%*%ep),r) #left eye, left end
  gLr1<-G(exp(x1%*%beta)*(br1%*%ep),r) #left eye, right end
  gLl2<-G(exp(x2%*%beta)*(bl2%*%ep),r) #right eye, left end
  gLr2<-G(exp(x2%*%beta)*(br2%*%ep),r)

  pgLl1 = pG(exp(x1%*%beta)*(bl1%*%ep),r)
  pgLr1 = pG(exp(x1%*%beta)*(br1%*%ep),r)
  pgLl2 = pG(exp(x2%*%beta)*(bl2%*%ep),r)
  pgLr2 = pG(exp(x2%*%beta)*(br2%*%ep),r)

  sl1 = exp(-gLl1)
  sr1 = exp(-gLr1)
  sl2 = exp(-gLl2)
  sr2 = exp(-gLr2)

  term1 <- (sl1 - sr1)
  term1 <- ifelse(indata1[,"status"] == 1 , term1, 0)
  term2 <- (sl1)
  term2 <- ifelse(indata1[,"status"] == 0 , term2, 0)
  l1 = (term1 + term2)

  term3 <- (sl2 - sr2)
  term3 <- ifelse(indata2[,"status"] == 1 , term3, 0)
  term4 <- (sl2)
  term4 <- ifelse(indata2[,"status"] == 0 , term4, 0)
  l2 <- (term3 + term4)

  lb1 = matrix(0,n,p)
  lb2 = matrix(0,n,p)
  sl1b = matrix(0,n,p)
  sr1b = matrix(0,n,p)
  sl2b = matrix(0,n,p)
  sr2b = matrix(0,n,p)

  for (i in 1:p){

    sl1b[,i] = -(sl1^(0+1))*(pgLl1*Ll1)*x1[,i]
    sr1b[,i] = -(sr1^(0+1))*(pgLr1*Lr1)*x1[,i]
    sl2b[,i] = -(sl2^(0+1))*(pgLl2*Ll2)*x2[,i]
    sr2b[,i] = -(sr2^(0+1))*(pgLr2*Lr2)*x2[,i]

    tmp1 <- ifelse(indata1[,"status"] == 1 , (sl1b[,i]-sr1b[,i])/l1, 0)
    tmp2 <- ifelse(indata1[,"status"] == 0 , sl1b[,i]/l1, 0)
    lb1[,i] = tmp1 + tmp2

    tmp1 <- ifelse(indata2[,"status"] == 1 , (sl2b[,i]-sr2b[,i])/l2, 0)
    tmp2 <- ifelse(indata2[,"status"] == 0 , sl2b[,i]/l2, 0)
    lb2[,i] = tmp1 + tmp2

  }

  lb = lb1+lb2

  sl1p1 = matrix(0,n,m+1)
  sr1p1 = matrix(0,n,m+1)
  sl2p2 = matrix(0,n,m+1)
  sr2p2 = matrix(0,n,m+1)
  lp1 = matrix(0,n,m+1)
  lp2 = matrix(0,n,m+1)


  for (i in 1:(m+1)) {

    sl1p1[,i] = -(sl1^(0+1))*pgLl1*exp(x1%*%beta)*bl1[,i]*ep[i]
    sr1p1[,i] = -(sr1^(0+1))*pgLr1*exp(x1%*%beta)*br1[,i]*ep[i]
    sl2p2[,i] = -(sl2^(0+1))*pgLl2*exp(x2%*%beta)*bl2[,i]*ep[i]
    sr2p2[,i] = -(sr2^(0+1))*pgLr2*exp(x2%*%beta)*br2[,i]*ep[i]

    tmp1 <- ifelse(indata1[,"status"] == 1 , (sl1p1[,i]-sr1p1[,i])/l1, 0)
    tmp2 <- ifelse(indata1[,"status"] == 0 , sl1p1[,i]/l1, 0)
    lp1[,i] = tmp1 + tmp2

    tmp1 <- ifelse(indata2[,"status"] == 1 , (sl2p2[,i]-sr2p2[,i])/l2, 0)
    tmp2 <- ifelse(indata2[,"status"] == 0 , sl2p2[,i]/l2, 0)
    lp2[,i] = tmp1 + tmp2

  }

  lp = lp1+lp2


  score = cbind(lb,lp)
  B = matrix(0, (p+1+1*m),(p+1+1*m))
  for (i in 1:n) {
    B = B + (score[i,]) %*% t(score[i,])
  }

  return(B)

}




# a function to calculate B matrix in the sandwich variance estimator: A^(-1) B A^(-1)
# for generalized robust score test under the marginal sieve model
B_mar_null <- function(para,m,p,n,x1,x2,bl1,br1,bl2,br2,indata1, indata2,r){

  phi<-para[1:(m+1)] #4:7,m=3
  beta<-para[(m+1+1):(m+1+p)] #p=3
  ep = cumsum(exp(phi))

  Ll1 = exp(x1%*%beta)*(bl1%*%ep)
  Lr1 = exp(x1%*%beta)*(br1%*%ep)
  Ll2 = exp(x2%*%beta)*(bl2%*%ep)
  Lr2 = exp(x2%*%beta)*(br2%*%ep)

  gLl1<-G(exp(x1%*%beta)*(bl1%*%ep),r) #left eye, left end
  gLr1<-G(exp(x1%*%beta)*(br1%*%ep),r) #left eye, right end
  gLl2<-G(exp(x2%*%beta)*(bl2%*%ep),r) #right eye, left end
  gLr2<-G(exp(x2%*%beta)*(br2%*%ep),r)

  pgLl1 = pG(exp(x1%*%beta)*(bl1%*%ep),r)
  pgLr1 = pG(exp(x1%*%beta)*(br1%*%ep),r)
  pgLl2 = pG(exp(x2%*%beta)*(bl2%*%ep),r)
  pgLr2 = pG(exp(x2%*%beta)*(br2%*%ep),r)

  sl1 = exp(-gLl1)
  sr1 = exp(-gLr1)
  sl2 = exp(-gLl2)
  sr2 = exp(-gLr2)

  term1 <- (sl1 - sr1)
  term1 <- ifelse(indata1[,"status"] == 1 , term1, 0)
  term2 <- (sl1)
  term2 <- ifelse(indata1[,"status"] == 0 , term2, 0)
  l1 = (term1 + term2)

  term3 <- (sl2 - sr2)
  term3 <- ifelse(indata2[,"status"] == 1 , term3, 0)
  term4 <- (sl2)
  term4 <- ifelse(indata2[,"status"] == 0 , term4, 0)
  l2 <- (term3 + term4)

  lb1 = matrix(0,n,p)
  lb2 = matrix(0,n,p)
  sl1b = matrix(0,n,p)
  sr1b = matrix(0,n,p)
  sl2b = matrix(0,n,p)
  sr2b = matrix(0,n,p)

  for (i in 1:p){

    sl1b[,i] = -(sl1^(0+1))*(pgLl1*Ll1)*x1[,i]
    sr1b[,i] = -(sr1^(0+1))*(pgLr1*Lr1)*x1[,i]
    sl2b[,i] = -(sl2^(0+1))*(pgLl2*Ll2)*x2[,i]
    sr2b[,i] = -(sr2^(0+1))*(pgLr2*Lr2)*x2[,i]

    tmp1 <- ifelse(indata1[,"status"] == 1 , (sl1b[,i]-sr1b[,i])/l1, 0)
    tmp2 <- ifelse(indata1[,"status"] == 0 , sl1b[,i]/l1, 0)
    lb1[,i] = tmp1 + tmp2

    tmp1 <- ifelse(indata2[,"status"] == 1 , (sl2b[,i]-sr2b[,i])/l2, 0)
    tmp2 <- ifelse(indata2[,"status"] == 0 , sl2b[,i]/l2, 0)
    lb2[,i] = tmp1 + tmp2

  }

  lb = lb1+lb2

  sl1p1 = matrix(0,n,m+1)
  sr1p1 = matrix(0,n,m+1)
  sl2p2 = matrix(0,n,m+1)
  sr2p2 = matrix(0,n,m+1)
  lp1 = matrix(0,n,m+1)
  lp2 = matrix(0,n,m+1)


  for (i in 1:(m+1)) {

    sl1p1[,i] = -(sl1^(0+1))*pgLl1*exp(x1%*%beta)*bl1[,i]*ep[i]
    sr1p1[,i] = -(sr1^(0+1))*pgLr1*exp(x1%*%beta)*br1[,i]*ep[i]
    sl2p2[,i] = -(sl2^(0+1))*pgLl2*exp(x2%*%beta)*bl2[,i]*ep[i]
    sr2p2[,i] = -(sr2^(0+1))*pgLr2*exp(x2%*%beta)*br2[,i]*ep[i]

    tmp1 <- ifelse(indata1[,"status"] == 1 , (sl1p1[,i]-sr1p1[,i])/l1, 0)
    tmp2 <- ifelse(indata1[,"status"] == 0 , sl1p1[,i]/l1, 0)
    lp1[,i] = tmp1 + tmp2

    tmp1 <- ifelse(indata2[,"status"] == 1 , (sl2p2[,i]-sr2p2[,i])/l2, 0)
    tmp2 <- ifelse(indata2[,"status"] == 0 , sl2p2[,i]/l2, 0)
    lp2[,i] = tmp1 + tmp2

  }

  lp = lp1+lp2


  score = cbind(lp,lb)
  B = matrix(0, (p+1+1*m),(p+1+1*m))
  for (i in 1:n) {
    B = B + (score[i,]) %*% t(score[i,])
  }

  return(B)

}

