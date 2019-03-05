
######## functions to generate marginal, conditional and joint survival probabilities ###########
######## used for estimating and plotting predicted survival probabilities ######


###### internally used prediction function #######

internal_predict.CopulaCenR <- function(object, class = "joint", newdata, evalPoints = 50, evalTimes1 = NULL, evalTimes2 = NULL, cond_time = NULL, cond_margin = 2,...) { # for conditional, only need to assign for evalTimes1


  # first screen the inputs #
  if (class(object) != "CopulaCenR") stop('object must be a CopulaCenR class object')
  if (!class %in% c("joint","conditional","marginal")) stop('class must be one of joint, conditional and marginal')
  if (!is.data.frame(newdata) | !"id" %in% colnames(newdata) | !"ind" %in% colnames(newdata)) stop('newdata must be a data frame with columns id, ind and var_list')

  # sort by id and ind
  newdata <- newdata[order(newdata$id, newdata$ind), ]
  newdata1 = newdata[newdata$ind==1, object$var_list]
  newdata2 = newdata[newdata$ind==2, object$var_list]

  newdata1 = as.numeric(newdata1)
  newdata2 = as.numeric(newdata2)


  if (is.null(evalTimes1) & is.null(evalPoints)) {
    evalPoints = 50 # set default
  }

  if (is.null(evalTimes1) & !is.null(evalPoints)) {
    # calculate grid length
    # xmin_1 = min(c(object$indata1[,"Left"], object$indata1[,"Right"]))
    # xmax_1 = max(c(object$indata1[,"Left"], object$indata1[!is.infinite(object$indata1[,"Right"]),"Right"]))
    # xmin_2 = min(c(object$indata2[,"Left"], object$indata2[,"Right"]))
    # xmax_2 = max(c(object$indata2[,"Left"], object$indata2[!is.infinite(object$indata2[,"Right"]),"Right"]))
    xmin_1 = min(c(object$indata1$Left, object$indata1$Right, object$indata1$obs_time))
    xmax_1 = max(c(object$indata1$Left, object$indata1$Right[!is.infinite(object$indata1$Right)], object$indata1$obs_time))
    xmin_2 = min(c(object$indata2$Left, object$indata2$Right, object$indata2$obs_time))
    xmax_2 = max(c(object$indata2$Left, object$indata2$Right[!is.infinite(object$indata2$Right)], object$indata2$obs_time))
    grid.length1 = (xmax_1 - xmin_1)/evalPoints
    grid.length2 = (xmax_2 - xmin_2)/evalPoints

    if (class == "marginal") {
      output1 <- m_copula(object = object, grid.length = grid.length1, newdata = newdata1)
      output2 <- m_copula(object = object, grid.length = grid.length2, newdata = newdata2)
      output <- list(grid1 = output1$grid, m1 = output1$m, grid2 = output2$grid, m2 = output2$m,...)
    }

    if (class == "conditional") {
      output <- cond_copula(object = object, grid.length = grid.length1, newdata1 = newdata1, newdata2 = newdata2, cond_time = cond_time, cond_margin = cond_margin)
    }

    if (class == "joint") {
      output <- surv2_copula(object = object, grid.length1 = grid.length1, grid.length2 = grid.length2, newdata1 = newdata1, newdata2 = newdata2)
    }

  }


  if (!is.null(evalTimes1)) {

    grid.length1 = NULL
    grid.length2 = NULL

    if (class == "marginal") {
      output1 <- m_copula(object = object, grid.length = grid.length1, newdata = newdata1, evalTimes = evalTimes1)
      output2 <- m_copula(object = object, grid.length = grid.length2, newdata = newdata2, evalTimes = evalTimes2)
      output <- list(grid1 = output1$grid, m1 = output1$m, grid2 = output2$grid, m2 = output2$m,...)

    }

    if (class == "conditional") {
      output <- cond_copula(object = object, grid.length = grid.length1, newdata1 = newdata1, newdata2 = newdata2, cond_time = cond_time, cond_margin = cond_margin, evalTimes = evalTimes1)
    }

    if (class == "joint") {
      output <- surv2_copula(object = object, grid.length1 = grid.length1, grid.length2 = grid.length2, newdata1 = newdata1, newdata2 = newdata2, evalTimes1 = evalTimes1, evalTimes2 = evalTimes2)
    }

  }



  return(output)
}


####### wrapper functions ##########
m_copula <- function(object, grid.length, newdata, evalTimes = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # IC, transformation model
  if (is.numeric(object$m)){

    output <- m_sieve(object, grid.length, newdata, evalTimes = evalTimes)

  }

  # IC, parametric margins
  else if (!is.numeric(object$m) & !("obs_time" %in% colnames(object$indata1)) ) {

    output <- m_ic(object, grid.length, newdata, evalTimes = evalTimes)

  }

  # RC, parametric margins
  else if (!is.numeric(object$m) & ("obs_time" %in% colnames(object$indata1)) ) {

    output <- m_rc(object, grid.length, newdata, evalTimes = evalTimes)

  }

  return(output)
}




cond_copula <- function(object, grid.length, newdata1, newdata2, cond_time, cond_margin = 2, evalTimes = NULL){ # length of newdata1/newdata2 is the same as p = # of covariates in object



  # IC, transformation model
  if (is.numeric(object$m)){

    output <- cond_sieve(object, grid.length, newdata1, newdata2, cond_time, cond_margin = cond_margin, evalTimes = evalTimes)

  }

  # IC, parametric margins
  else if (!is.numeric(object$m) & !("obs_time" %in% colnames(object$indata1)) ) {

    output <- cond_ic(object, grid.length, newdata1, newdata2, cond_time, cond_margin = cond_margin, evalTimes = evalTimes)

  }

  # RC, parametric margins
  else if (!is.numeric(object$m) & ("obs_time" %in% colnames(object$indata1)) ) {

    output <- cond_rc(object, grid.length, newdata1, newdata2, cond_time, cond_margin = cond_margin, evalTimes = evalTimes)

  }

  return(output)
}



surv2_copula <- function(object, grid.length1, grid.length2, newdata1, newdata2, evalTimes1 = NULL, evalTimes2 = NULL){ # length of newdata1/newdata2 is the same as p = # of covariates in object



  # IC, transformation model
  if (is.numeric(object$m)){

    output <- surv2_sieve(object, grid.length1, grid.length2, newdata1, newdata2, evalTimes1 = evalTimes1, evalTimes2 = evalTimes2)

  }

  # IC, parametric margins
  else if (!is.numeric(object$m) & !("obs_time" %in% colnames(object$indata1)) ) {

    output <- surv2_ic(object, grid.length1, grid.length2, newdata1, newdata2, evalTimes1 = evalTimes1, evalTimes2 = evalTimes2)

  }

  # RC, parametric margins
  else if (!is.numeric(object$m) & ("obs_time" %in% colnames(object$indata1)) ) {

    output <- surv2_rc(object, grid.length1, grid.length2, newdata1, newdata2, evalTimes1 = evalTimes1, evalTimes2 = evalTimes2)

  }

  return(output)
}





######## utilising functions ######

### marginal survival, ic, parametric margins ###
m_ic <- function(object, grid.length, margin, evalTimes = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  l1 <- min(object$indata1[,"Left"], object$indata1[,"Right"])
  u1 <- max(object$indata1[,"Left"], object$indata1[,"Right"][is.finite(object$indata1[,"Right"])])
  l2 <-  min(object$indata2[,"Left"], object$indata2[,"Right"])
  u2 <- max(object$indata2[,"Left"], object$indata2[,"Right"][is.finite(object$indata2[,"Right"])])

  copula <- object$copula
  m.dist <- object$m.dist
  p <- dim(object$x1)[2]


  # create grids for two events; assume same grid for both events
  l <- min(l1, l2)
  u <- max(u1, u2)
  if (is.null(evalTimes)) {grid1 <- seq(l,u,grid.length)}
  if (!is.null(evalTimes)) {grid1 <- evalTimes}

  margin1 <- margin

  m_ic <-  as.numeric()

  baseline <- object$estimate[1:2]
  beta <- object$estimate[3:(2+p)]

  for (i in 1:length(grid1)){

    if (m.dist == "Weibull") {
      lambda <- baseline[1]
      k <- baseline[2]
      S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
    }

    if (m.dist == "Loglogistic") {
      lambda <- baseline[1]
      k <- baseline[2]
      S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
    }

    if (m.dist == "Gompertz") {
      a <- baseline[1]
      b <- baseline[2]
      S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
    }

    m_ic <- c(m_ic, S1)

  }

  output <- list(grid=grid1, m = m_ic)
  return(output)

}




### marginal survival, ic, transformation sieve margins ###
m_sieve <- function(object, grid.length, margin, evalTimes = NULL){ # m_margin = 2 means mitioing on 2nd eye



  # data
  l <- object$l
  u <- object$u
  copula <- object$copula
  r <- object$r
  m <- object$m
  p <- dim(object$x1)[2]

  # create grids and margins
  if (is.null(evalTimes)) {grid1 <- seq(l,u,grid.length)}
  if (!is.null(evalTimes)) {grid1 <- evalTimes}

  margin1 <- margin

  # generate bernstein polynomials based on grids
  BL = matrix(0,nrow = length(grid1),ncol = m+1)
  for (i in 0:m) {
    BL[,(i+1)] = bern(i,m,l,u,grid1) # used for 1st eye grid calculations
  }


  m_sieve <- as.numeric()

  for (i in 1:nrow(BL)){
    beta <- object$estimate[1:p]
    phi <- object$estimate[(p+1):(p+1+m)]
    ep<-cumsum(exp(phi))
    H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
    S1 = exp(-H1)
    m_sieve <- c(m_sieve, S1)
  }


  output <- list(grid=grid1, m = m_sieve)
  return(output)

}




### marginal survival, rc, parametric margins ###
m_rc <- function(object, grid.length, margin, evalTimes = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  l1 <- min(object$indata1[,"obs_time"])
  u1 <- max(object$indata1[,"obs_time"][is.finite(object$indata1[,"obs_time"])])
  l2 <-  min(object$indata2[,"obs_time"])
  u2 <- max(object$indata2[,"obs_time"][is.finite(object$indata2[,"obs_time"])])

  copula <- object$copula
  m.dist <- object$m.dist
  n.cons <- object$n.cons # NULL if not piecewise
  quantiles <- object$quantiles # NULL if not piecewise
  p <- dim(object$x1)[2]

  # create grids for two events; assume same grid for both events
  l <- min(l1, l2)
  u <- max(u1, u2)

  if (is.null(evalTimes)) {grid1 <- seq(l,u,grid.length)}
  if (!is.null(evalTimes)) {grid1 <- evalTimes}

  margin1 <- margin


  m_rc <-  as.numeric()

  baseline <- object$estimate[1:2]
  beta <- object$estimate[3:(2+p)]

  for (i in 1:length(grid1)){

    if (m.dist == "Weibull") {
      lambda <- baseline[1]
      k <- baseline[2]
      S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
    }

    if (m.dist == "Loglogistic") {
      lambda <- baseline[1]
      k <- baseline[2]
      S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
    }

    if (m.dist == "Gompertz") {
      a <- baseline[1]
      b <- baseline[2]
      S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
    }


    if (m.dist == "Piecewise") {

      Lambda1<-0 # left eye

      # different parameterization
      for (k in 1:n.cons)
      {
        Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
      }

      beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

      S1<-exp(-Lambda1*exp(margin1%*%beta))
    }


    m_rc <- c(m_rc, S1)

  }


  output <- list(grid=grid1, m = m_rc)
  return(output)

}



### conditional survival, ic, sieve margins ###
cond_sieve <- function(object, grid.length, margin1, margin2, cond_time, cond_margin = 2, evalTimes = NULL){ # cond_margin = 2 means conditioing on 2nd eye



  # data
  l <- object$l
  u <- object$u
  copula <- object$copula
  r <- object$r
  m <- object$m
  p <- dim(object$x1)[2]

  # create grids for two events; assume same grid for both events
  if (is.null(evalTimes)) {grid1 <- seq(0,(u-cond_time),grid.length) + cond_time}
  if (!is.null(evalTimes)) {grid1 <- evalTimes + cond_time}

    grid2 <- cond_time


  # generate bernstein polynomials based on grids
  BL = matrix(0,nrow = length(grid1),ncol = m+1)
  for (i in 0:m) {
    BL[,(i+1)] = bern(i,m,l,u,grid1) # used for 1st eye grid calculations
  }

  BR = matrix(0,nrow = length(grid2),ncol = m+1)
  for (i in 0:m) {
    BR[,(i+1)] = bern(i,m,l,u,grid2) # used for 2nd eye grid calculations
  }

  ### first, calculate the numerator of Pr(T1 >5+t1 and T2 < 5)
  surv2 <- matrix(NA,nrow=length(grid1),ncol=length(grid2))
  rownames(surv2) <- grid1
  colnames(surv2) <- grid2

  if (copula == "Copula2") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        alpha <- object$estimate[p+1+m+1]
        kappa <- object$estimate[p+1+m+1+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1 - (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)

      }
    }
  }


  if (copula == "Clayton") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1 - (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }


  if (copula == "Gumbel") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1 - exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }

  if (copula == "Frank") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1 - (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1 - S1*S2/(1-eta*(1-S1)(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1 - (1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta)) # Joe
      }
    }
  }




  ### Then, calculate the scalar value of joint probability of Pr(T1 > 5 and T2< 5)
  grid1_2 = cond_time
  grid2_2 = cond_time
  BL=BR=matrix(0,nrow = length(grid1_2),ncol = m+1)
  for (i in 0:m) {
    BL[,(i+1)] = bern(i,m,l,u,grid1_2)
    BR[,(i+1)] = bern(i,m,l,u,grid2_2)
  }
  denumerator <- matrix(NA,nrow=length(grid1_2),ncol=length(grid2_2))


  if (copula == "Copula2") {

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        alpha <- object$estimate[p+1+m+1]
        kappa <- object$estimate[p+1+m+1+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        denumerator[i,j] <- S1 - (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }


  if (copula == "Clayton") {

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        denumerator[i,j] <- S1 - (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }


  if (copula == "Gumbel") {

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        denumerator[i,j] <- S1 - exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }

  if (copula == "Frank") {

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        denumerator[i,j] <- S1 - (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        denumerator[i,j] <- S1 - S1*S2/(1-eta*(1-S1)(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        denumerator[i,j] <- S1 - (1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta)) # Joe
      }
    }
  }


  condition <- surv2%*%solve(denumerator)

  if (cond_margin == 2){
    output <- list(grid1=grid1, grid2=grid2, condition = condition) # vector grid is greater than cond_time, smaller than u
  }
  if (cond_margin == 1){
    output <- list(grid1=grid2, grid2=grid1, condition = condition) # vector grid is greater than cond_time, smaller than u
  }
  return(output)

}


### conditional survival, ic, parametric margins ###
cond_ic <- function(object, grid.length, margin1, margin2, cond_time, cond_margin = 2, evalTimes = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  u1 <- max(object$indata1[,"Left"], object$indata1[,"Right"][is.finite(object$indata1[,"Right"])])
  u2 <- max(object$indata2[,"Left"], object$indata2[,"Right"][is.finite(object$indata2[,"Right"])])

  copula <- object$copula
  m.dist <- object$m.dist
  p <- dim(object$x1)[2]

  # create grids for two events; assume same grid for both events
  if (cond_margin == 2){
    if (is.null(evalTimes)) {grid1 <- seq(0,(u1-cond_time),grid.length) + cond_time}
    if (!is.null(evalTimes)) {grid1 <- evalTimes + cond_time}
    grid2 <- cond_time
  }
  if (cond_margin == 1){
    if (is.null(evalTimes)) {grid1 <- seq(0,(u2-cond_time),grid.length) + cond_time}
    if (!is.null(evalTimes)) {grid1 <- evalTimes + cond_time}
    grid2 <- cond_time
  }

  ### first, calculate the numerator of Pr(T1 >5+t1 and T2 < 5)
  surv2 <- matrix(NA,nrow=length(grid1),ncol=length(grid2))
  rownames(surv2) <- grid1
  colnames(surv2) <- grid2

  if (copula == "Clayton") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))

        }

        surv2[i,j] <- S1 - (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }



  if (copula == "Gumbel") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }



  if (copula == "Frank") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }


  if (copula == "AMH") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - S1*S2/(1-eta*(1-S1)(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))

        }

        surv2[i,j] <- S1 - (1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta)) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    alpha <- object$estimate[(length(object$estimate)-1)]
    kappa <- object$estimate[(length(object$estimate))]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }


  ### Then, calculate the scalar value of joint probability of Pr(T1 > 5 and T2< 5)
  grid1_2 = cond_time
  grid2_2 = cond_time
  denumerator <- matrix(NA,nrow=length(grid1_2),ncol=length(grid2_2))


  if (copula == "Clayton") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))

        }

        denumerator[i,j] <- S1 - (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }



  if (copula == "Gumbel") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }



  if (copula == "Frank") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }


  if (copula == "AMH") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - S1*S2/(1-eta*(1-S1)(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))

        }

        denumerator[i,j] <- S1 - (1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta)) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    alpha <- object$estimate[(length(object$estimate)-1)]
    kappa <- object$estimate[(length(object$estimate))]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }


  condition <- surv2%*%solve(denumerator)

  if (cond_margin == 2){
    output <- list(grid1=grid1, grid2=grid2, condition = condition) # vector grid is greater than cond_time, smaller than u
  }
  if (cond_margin == 1){
    output <- list(grid1=grid2, grid2=grid1, condition = condition) # vector grid is greater than cond_time, smaller than u
  }
  return(output)
}



### conditional survival, rc, parametric margins ###
cond_rc <- function(object, grid.length, margin1, margin2, cond_time, cond_margin = 2, evalTimes = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  u1 <- max(object$indata1[,"obs_time"][is.finite(object$indata1[,"obs_time"])])
  u2 <- max(object$indata2[,"obs_time"][is.finite(object$indata2[,"obs_time"])])

  copula <- object$copula
  m.dist <- object$m.dist
  n.cons <- object$n.cons # NULL if not piecewise
  quantiles <- object$quantiles # NULL if not piecewise
  p <- dim(object$x1)[2]


  # create grids for two events; assume same grid for both events
  if (cond_margin == 2){
    if (is.null(evalTimes)) {grid1 <- seq(0,(u1-cond_time),grid.length) + cond_time}
    if (!is.null(evalTimes)) {grid1 <- evalTimes + cond_time}
    grid2 <- cond_time
  }
  if (cond_margin == 1){
    if (is.null(evalTimes)) {grid1 <- seq(0,(u2-cond_time),grid.length) + cond_time}
    if (!is.null(evalTimes)) {grid1 <- evalTimes + cond_time}
    grid2 <- cond_time
  }


  ### first, calculate the numerator of Pr(T1 >5+t1 and T2 < 5)
    surv2 <- matrix(NA,nrow=length(grid1),ncol=length(grid2))
    rownames(surv2) <- grid1
    colnames(surv2) <- grid2





  if (copula == "Clayton") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }


        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }


        surv2[i,j] <- S1 - (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }



  if (copula == "Gumbel") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }



  if (copula == "Frank") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - S1*S2/(1-eta*(1-S1)*(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))

        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1 - (1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta)) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    alpha <- object$estimate[(length(object$estimate)-1)]
    kappa <- object$estimate[(length(object$estimate))]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          alpha<-object$estimate[(n.cons+1)] # association
          kappa<-object$estimate[(n.cons+2)] # association
          beta<-object$estimate[(n.cons+3):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1- (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }


  ### Then, calculate the scalar value of joint probability of Pr(T1 > 5 and T2< 5)
  grid1_2 = cond_time
  grid2_2 = cond_time
  denumerator <- matrix(NA,nrow=length(grid1_2),ncol=length(grid2_2))


  if (copula == "Clayton") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }


        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1_2[i])), apply(cbind(rep(quantiles[k+1],length(grid1_2[i])),grid1_2[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2_2[j])), apply(cbind(rep(quantiles[k+1],length(grid2_2[j])),grid2_2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }


        denumerator[i,j] <- S1 - (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }



  if (copula == "Gumbel") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1_2[i])), apply(cbind(rep(quantiles[k+1],length(grid1_2[i])),grid1_2[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2_2[j])), apply(cbind(rep(quantiles[k+1],length(grid2_2[j])),grid2_2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }



  if (copula == "Frank") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1_2[i])), apply(cbind(rep(quantiles[k+1],length(grid1_2[i])),grid1_2[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2_2[j])), apply(cbind(rep(quantiles[k+1],length(grid2_2[j])),grid2_2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1_2[i])), apply(cbind(rep(quantiles[k+1],length(grid1_2[i])),grid1_2[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2_2[j])), apply(cbind(rep(quantiles[k+1],length(grid2_2[j])),grid2_2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - S1*S2/(1-eta*(1-S1)*(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))

        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1_2[i])), apply(cbind(rep(quantiles[k+1],length(grid1_2[i])),grid1_2[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2_2[j])), apply(cbind(rep(quantiles[k+1],length(grid2_2[j])),grid2_2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1 - (1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta)) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    alpha <- object$estimate[(length(object$estimate)-1)]
    kappa <- object$estimate[(length(object$estimate))]

    for (i in 1:nrow(denumerator)){
      for (j in 1:ncol(denumerator)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1_2[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2_2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1_2[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2_2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1_2[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2_2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1_2[i])), apply(cbind(rep(quantiles[k+1],length(grid1_2[i])),grid1_2[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2_2[j])), apply(cbind(rep(quantiles[k+1],length(grid2_2[j])),grid2_2[j]),1,min) - quantiles[k]), 1, max)
          }

          alpha<-object$estimate[(n.cons+1)] # association
          kappa<-object$estimate[(n.cons+2)] # association
          beta<-object$estimate[(n.cons+3):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        denumerator[i,j] <- S1- (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }


  condition <- surv2%*%solve(denumerator)

  if (cond_margin == 2){
    output <- list(grid1=grid1, grid2=grid2, condition = condition) # vector grid is greater than cond_time, smaller than u
  }
  if (cond_margin == 1){
    output <- list(grid1=grid2, grid2=grid1, condition = condition) # vector grid is greater than cond_time, smaller than u
  }
  return(output)


}



### joint survival, ic, sieve margins ###
surv2_sieve <- function(object, grid.length1, grid.length2, margin1, margin2, evalTimes1 = NULL, evalTimes2 = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  l <- object$l
  u <- object$u
  copula <- object$copula
  r <- object$r
  m <- object$m
  p <- dim(object$x1)[2]

  # create grids for two events; assume same grid for both events
  if (is.null(evalTimes1)) {grid1 = seq(l,u,grid.length1); grid2 = seq(l,u,grid.length2)}
  if (!is.null(evalTimes1)) {grid1 = evalTimes1; grid2 = evalTimes2}

  # generate bernstein polynomials based on grids
  BL=matrix(0,nrow = length(grid1),ncol = m+1)
  BR=matrix(0,nrow = length(grid2),ncol = m+1)

  for (i in 0:m) {
    BL[,(i+1)] = bern(i,m,l,u,grid1)
    BR[,(i+1)] = bern(i,m,l,u,grid2)
  }

  # joint survival probability
  surv2 <- matrix(NA,nrow=length(grid1),ncol=length(grid2))
  rownames(surv2) <- grid1
  colnames(surv2) <- grid2


  if (copula == "Clayton") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }

  if (copula == "Gumbel") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }

  if (copula == "Frank") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- S1*S2/(1-eta*(1-S1)(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        eta <- object$estimate[p+1+m+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- 1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){
        beta <- object$estimate[1:p]
        phi <- object$estimate[(p+1):(p+1+m)]
        alpha <- object$estimate[p+1+m+1]
        kappa <- object$estimate[p+1+m+1+1]
        ep<-cumsum(exp(phi))
        H1<-G(exp(margin1%*%beta)*(BL[i,]%*%ep), r)
        H2<-G(exp(margin2%*%beta)*(BR[j,]%*%ep), r)
        S1 = exp(-H1)
        S2 = exp(-H2)
        surv2[i,j] <- (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }

  output <- list(grid1=grid1, grid2=grid2, surv2 = t(surv2))
  return(output)

}


### joint survival, ic, parametric margins ###

surv2_ic <- function(object, grid.length1, grid.length2, margin1, margin2, evalTimes1 = NULL, evalTimes2 = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  l1 <- min(object$indata1[,"Left"], object$indata1[,"Right"])
  u1 <- max(object$indata1[,"Left"], object$indata1[,"Right"][is.finite(object$indata1[,"Right"])])
  l2 <-  min(object$indata2[,"Left"], object$indata2[,"Right"])
  u2 <- max(object$indata2[,"Left"], object$indata2[,"Right"][is.finite(object$indata2[,"Right"])])

  copula <- object$copula
  m.dist <- object$m.dist
  p <- dim(object$x1)[2]

  # create grids for two events; assume same grid for both events
  if (is.null(evalTimes1)) {grid1 = seq(l1,u1,grid.length1); grid2 = seq(l2,u2,grid.length2)}
  if (!is.null(evalTimes1)) {grid1 = evalTimes1; grid2 = evalTimes2}


  # joint survival probability
  surv2 <- matrix(NA,nrow=length(grid1),ncol=length(grid2))
  rownames(surv2) <- grid1
  colnames(surv2) <- grid2


  if (copula == "Clayton") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))

        }

        surv2[i,j] <- (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }



  if (copula == "Gumbel") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }



  if (copula == "Frank") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1*S2/(1-eta*(1-S1)(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))

        }

        surv2[i,j] <- 1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    alpha <- object$estimate[(length(object$estimate)-1)]
    kappa <- object$estimate[(length(object$estimate))]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        surv2[i,j] <- (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }

  output <- list(grid1=grid1, grid2=grid2, surv2 = t(surv2))
  return(output)

}



### joint survival, rc, parametric margins ###

surv2_rc <- function(object, grid.length1, grid.length2, margin1, margin2, evalTimes1 = NULL, evalTimes2 = NULL){ # length of margin1/margin2 is the same as p = # of covariates in object



  # data
  l1 <- min(object$indata1[,"obs_time"])
  u1 <- max(object$indata1[,"obs_time"][is.finite(object$indata1[,"obs_time"])])
  l2 <-  min(object$indata2[,"obs_time"])
  u2 <- max(object$indata2[,"obs_time"][is.finite(object$indata2[,"obs_time"])])

  copula <- object$copula
  m.dist <- object$m.dist
  n.cons <- object$n.cons # NULL if not piecewise
  quantiles <- object$quantiles # NULL if not piecewise
  p <- dim(object$x1)[2]

  # create grids for two events; assume same grid for both events
  if (is.null(evalTimes1)) {grid1 = seq(l1,u1,grid.length1); grid2 = seq(l2,u2,grid.length2)}
  if (!is.null(evalTimes1)) {grid1 = evalTimes1; grid2 = evalTimes2}

  # joint survival probability
  surv2 <- matrix(NA,nrow=length(grid1),ncol=length(grid2))
  rownames(surv2) <- grid1
  colnames(surv2) <- grid2


  if (copula == "Clayton") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }


        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }


        surv2[i,j] <- (S1^(-eta)+S2^(-eta)-1)^(-1/eta)
      }
    }
  }



  if (copula == "Gumbel") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- exp(-((-log(S1))^eta + (-log(S2))^eta)^(1/eta))

      }
    }
  }



  if (copula == "Frank") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- (1/eta) * log(1 + (exp(eta*S1)-1)*(exp(eta*S2)-1)/(exp(eta)-1)) # frank
      }
    }
  }

  if (copula == "AMH") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- S1*S2/(1-eta*(1-S1)*(1-S2)) # AMH
      }
    }
  }


  if (copula == "Joe") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    eta <- object$estimate[length(object$estimate)]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))

        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          eta<-object$estimate[n.cons+1] # association
          beta<-object$estimate[(n.cons+2):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- 1 - ((1-S1)^eta + (1-S2)^eta - ((1-S1)^eta)*((1-S2)^eta) )^(1/eta) # Joe
      }
    }
  }


  if (copula == "Copula2") {

    baseline <- object$estimate[1:2]
    beta <- object$estimate[3:(2+p)]
    alpha <- object$estimate[(length(object$estimate)-1)]
    kappa <- object$estimate[(length(object$estimate))]

    for (i in 1:nrow(surv2)){
      for (j in 1:ncol(surv2)){

        if (m.dist == "Weibull") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1<-exp(-(grid1[i]/lambda)^k*exp(margin1%*%beta))
          S2<-exp(-(grid2[j]/lambda)^k*exp(margin2%*%beta))
        }

        if (m.dist == "Loglogistic") {
          lambda <- baseline[1]
          k <- baseline[2]
          S1 <- (1+((grid1[i]/lambda)^k)*exp(margin1%*%beta))^(-1)
          S2 <- (1+((grid2[j]/lambda)^k)*exp(margin2%*%beta))^(-1)
        }

        if (m.dist == "Gompertz") {
          a <- baseline[1]
          b <- baseline[2]
          S1<-exp(b/a*(1-exp(a*grid1[i]))*exp(margin1%*%beta))
          S2<-exp(b/a*(1-exp(a*grid2[j]))*exp(margin2%*%beta))
        }

        if (m.dist == "Piecewise") {

          Lambda1<-0 # left eye
          Lambda2<-0 # right eye

          # different parameterization
          for (k in 1:n.cons)
          {
            Lambda1<-Lambda1+ object$estimate[k] * apply(cbind(rep(0,length(grid1[i])), apply(cbind(rep(quantiles[k+1],length(grid1[i])),grid1[i]),1,min) - quantiles[k]), 1, max)
            Lambda2<-Lambda2+ object$estimate[k] * apply(cbind(rep(0,length(grid2[j])), apply(cbind(rep(quantiles[k+1],length(grid2[j])),grid2[j]),1,min) - quantiles[k]), 1, max)
          }

          alpha<-object$estimate[(n.cons+1)] # association
          kappa<-object$estimate[(n.cons+2)] # association
          beta<-object$estimate[(n.cons+3):length(object$estimate)] # coefficients

          S1<-exp(-Lambda1*exp(margin1%*%beta))
          S2<-exp(-Lambda2*exp(margin2%*%beta))
        }

        surv2[i,j] <- (1+(((S1)^(-1/kappa)-1)^(1/alpha) + ((S2)^(-1/kappa)-1)^(1/alpha))^(alpha))^(-kappa)
      }
    }
  }

  output <- list(grid1=grid1, grid2=grid2, surv2 = t(surv2))
  return(output)

}













