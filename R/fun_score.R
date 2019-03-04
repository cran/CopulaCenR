#' @importFrom stats pchisq

############## calcualte score statistics and perform generalized score test ##############

### rc, parametric margins
rc_par_copula_score <- function(object, var_score){

  obj <- object
  copula <- obj$copula
  m.dist <- obj$m.dist
  indata1 <- obj$indata1 # raw dataset containing NULL and Alternative variables
  indata2 <- obj$indata2
  var_list <- obj$var_list

  # for piecewise
  quantiles <- obj$quantiles
  n.cons <- obj$n.cons

  x1.1 <- obj$x1
  x2.1 <- obj$x2

  tmp1 <- get_covariates_rc(indata1, var_score)
  x1 <- cbind(x1.1, tmp1$x)
  tmp2 <- get_covariates_rc(indata2, var_score)
  x2 <- cbind(x2.1, tmp2$x)
  x1 <- as.matrix(x1,dim(x1)[1])
  x2 <- as.matrix(x2,dim(x2)[1])


  p <- dim(x1)[2]
  p1 <- length(var_list)
  p2 <- dim(x1)[2] - p1



  if (m.dist != "Piecewise") {


    score = grad(rc_copula_log_lik, c(obj$estimates[1:(p1+2)],rep(0,p2),obj$estimates[(p1+2+1):length(obj$estimates)]),
                 x1=x1, x2=x2,indata1=indata1,indata2=indata2,
                 copula = copula, m.dist = m.dist)

    hes = hessian(rc_copula_log_lik, c(obj$estimates[1:(p1+2)],rep(0,p2),obj$estimates[(p1+2+1):length(obj$estimates)]),
                  x1=x1, x2=x2,indata1=indata1,indata2=indata2,
                  copula = copula, m.dist = m.dist)

    stat = t(score) %*% pseudoinverse(hes) %*% score
    pvalue = pchisq(stat,p2,lower.tail=F)
    output = c(stat, pvalue)
    names(output) = c("stat", "pvalue")

  }


  if (m.dist == "Piecewise") {

    score = grad(rc_copula_log_lik,c(obj$estimates,rep(0,p2)),
                 x1=x1, x2=x2, indata1=indata1,indata2=indata2,quantiles=quantiles,
                 copula = copula, m.dist = m.dist)

    hes = hessian(rc_copula_log_lik,c(obj$estimates,rep(0,p2)),
                  x1=x1, x2=x2, indata1=indata1,indata2=indata2,quantiles=quantiles,
                  copula = copula, m.dist = m.dist)

    stat = t(score) %*% pseudoinverse(hes) %*% score
    pvalue = pchisq(stat,p2,lower.tail=F)
    output = c(stat, pvalue)
    names(output) = c("stat", "pvalue")

  }

  return(output)
}



### ic, parametric margins
ic_par_copula_score <- function(object, var_score){



  copula <- object$copula
  m.dist <- object$m.dist
  indata1 <- object$indata1 # raw dataset containing NULL and Alternative variables
  indata2 <- object$indata2
  var_list <- object$var_list
  x1.1 <- object$x1
  x2.1 <- object$x2


  # split and calculate time ranges
  t1_left<-indata1[,"Left"]
  t1_right<-indata1[,"Right"]
  t2_left<-indata2[,"Left"]
  t2_right<-indata2[,"Right"]

  # categorical to factor
  tmp1 <- var_transform_score(indata1, var_score, x1.1)
  tmp2 <- var_transform_score(indata2, var_score, x2.1)
  x1 <- tmp1$x
  x2 <- tmp2$x

  # transform data frame into matrix for matrix operations
  x1 <- as.matrix(x1,dim(x1)[1])
  x2 <- as.matrix(x2,dim(x2)[1])

  # NEW design matrix = var_list (Null model) + length of tested variables
  p1 <- length(var_list)
  p2 <- dim(x1)[2] - p1
  p <- dim(x1)[2]


  if (copula != "Copula2") {
    score = grad(ic_copula_log_lik_param,x0=c(object$estimates[1:(p1+2)],rep(0,p2),object$estimates[p1+2+1]),
                 p=p, x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                 copula = copula, m.dist = m.dist)

    hes = hessian(ic_copula_log_lik_param,x0=c(object$estimates[1:(p1+2)],rep(0,p2),object$estimates[p1+2+1]),
                  p=p, x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                  copula = copula, m.dist = m.dist)
  }


  if (copula == "Copula2") {
    score = grad(ic_copula_log_lik_param,x0=c(object$estimates[1:(p1+2)],rep(0,p2),object$estimates[(p1+2+1):(p1+2+2)]),
                 p=p, x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                 copula = copula, m.dist = m.dist)

    hes = hessian(ic_copula_log_lik_param,x0=c(object$estimates[1:(p1+2)],rep(0,p2),object$estimates[(p1+2+1):(p1+2+2)]),
                  p=p, x1=x1, x2=x2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right, indata1=indata1, indata2=indata2,
                  copula = copula, m.dist = m.dist)
  }

  stat = t(score) %*% pseudoinverse(hes) %*% score
  pvalue = pchisq(stat,p2,lower.tail=F)
  output = c(stat, pvalue)
  names(output) = c("stat", "pvalue")


  return(output)
}




### ic, sieve margins
ic_sp_copula_score <- function(object, var_score){



  # data
  copula <- object$copula
  m <- object$m
  r <- object$r
  indata1 <- object$indata1 # raw dataset containing NULL and Alternative variables
  indata2 <- object$indata2
  var_list <- object$var_list
  x1.1 <- object$x1
  x2.1 <- object$x2

  # BP
  bl1 <- object$bl1 #left eye, left end
  br1 <- object$br1 #left eye, right end
  bl2 <- object$bl2 #right eye, left end
  br2 <- object$br2 #right eye, right end

  # split and calculate time ranges
  t1_left<-indata1[,"Left"]
  t1_right<-indata1[,"Right"]
  t2_left<-indata2[,"Left"]
  t2_right<-indata2[,"Right"]

  # categorical to factor
  tmp1 <- var_transform_score(indata1, var_score, x1.1)
  tmp2 <- var_transform_score(indata2, var_score, x2.1)
  x1 <- tmp1$x
  x2 <- tmp2$x


  # transform data frame into matrix for matrix operations
  x1 <- as.matrix(x1,dim(x1)[1])
  x2 <- as.matrix(x2,dim(x2)[1])

  # NEW design matrix = var_list (Null model) + colnames(x1_factor or x2_factor)
  p1 <- length(var_list)
  p2 <- dim(x1)[2] - p1
  p <- dim(x1)[2]


  if (copula != "Copula2") {
    score = grad(ic_copula_log_lik_sieve,x0=c(object$estimates[1:p1],rep(0,p2),object$estimates[(p1+1):(p1+m+1+1)]),
                 p=p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                 copula = copula)

    hes = hessian(ic_copula_log_lik_sieve,x0=c(object$estimates[1:p1],rep(0,p2),object$estimates[(p1+1):(p1+m+1+1)]),
                  p=p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                  copula = copula)
  }


  if (copula == "Copula2") {
    score = grad(ic_copula_log_lik_sieve,x0=c(object$estimates[1:p1],rep(0,p2),object$estimates[(p1+1):(p1+m+1+2)]),
                 p=p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                 copula = copula)

    hes = hessian(ic_copula_log_lik_sieve,x0=c(object$estimates[1:p1],rep(0,p2),object$estimates[(p1+1):(p1+m+1+2)]),
                  p=p, m=m, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, indata1=indata1, indata2=indata2, r=r,
                  copula = copula)
  }

  stat = t(score) %*% pseudoinverse(hes) %*% score
  pvalue = pchisq(stat,p2,lower.tail=F)
  output = c(stat, pvalue)
  names(output) = c("stat", "pvalue")

  return(output)
}
