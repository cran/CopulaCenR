
###### some data proprocessing functions before model fitting ######

### tranform categorical values into factors
var_transform <- function(indata, var_list){

  # scan class of each covariate
  list_all <- split(var_list,sapply(indata[1,var_list], function(x) paste(class(x), collapse=" ")))
  list_factor <- c(list_all$factor,list_all$character)
  list_num <- c(list_all$numeric,list_all$integer)

  # extract numeric and factor/character variables
  x_num <- data.frame(indata[,list_num])
  colnames(x_num) <- list_num
  x_factor <- data.frame(indata[,list_factor])
  colnames(x_factor) <- list_factor

  # transform categorical to factor
  if (ncol(x_factor)!=0) {
    for (i in 1:ncol(x_factor)){

      tmp <- unique(x_factor[,i])[order(unique(x_factor[,i]))]
      for(level in tmp[2:length(tmp)]){
        x_factor[paste(colnames(x_factor)[i], level, sep = "")] <- ifelse(x_factor[,i] == level, 1, 0)
      }
    }
    tmp <- colnames(x_factor)[-(1:length(list_factor))]
    x_factor <- data.frame(x_factor[,-(1:length(list_factor))])
    colnames(x_factor) <- tmp

  }

  # combine
  x <- cbind(x_num,x_factor)
  var_list <- colnames(x)

  return(list(x=x, var_list=var_list))
}


### tranform categorical values into factors for score test
var_transform_score <- function(indata, var_score, x0) {

  # scan class of each covariate in var_score
  list_all <- split(var_score,sapply(indata[1,var_score], function(x) paste(class(x), collapse=" ")))
  list_factor <- c(list_all$factor,list_all$character)
  list_num <- c(list_all$numeric,list_all$integer)

  # if var_score has numeric
  if (length(list_num)!=0) {
    x_num <- data.frame(indata[,list_num])
    colnames(x_num) <- list_num
    x0 <- cbind(x0, x_num) # generate new x
  }

  # if var_score has factor
  if (length(list_factor)!=0) { # generate dummy variables if x_factor is not empty
    x_factor <- data.frame(indata[,list_factor])
    colnames(x_factor) <- list_factor
    for (i in 1:ncol(x_factor)){
      for(level in unique(x_factor[,i])[2:length(unique(x_factor[,i]))]){ #use first level as reference
        x_factor[paste(colnames(x_factor)[i], level, sep = "")] <- ifelse(x_factor[,i] == level, 1, 0)
      }
    }

    tmp <- colnames(x_factor)[-(1:length(list_factor))]
    x_factor <- data.frame(x_factor[,-(1:length(list_factor))])
    colnames(x_factor) <- tmp
    x0 <- cbind(x0,x_factor) # generate new x
  }

  return(list(x=x0))


}



# ic, parametric
data_process <- function(data, var_list) {

  indata1 <- data[data[,"ind"]==1, ]
  indata2 <- data[data[,"ind"]==2, ]
  # indata1 <- subset(data, ind==1)
  # indata2 <- subset(data, ind==2)
  t1_left <- indata1[,"Left"]
  t1_right <- indata1[,"Right"]
  t2_left <- indata2[,"Left"]
  t2_right <- indata2[,"Right"]

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

  # for icenReg
  x <- data.frame(Left=c(t1_left,t2_left),Right=c(t1_right,t2_right),rbind(x1,x2))

  return(list(indata1=indata1, indata2=indata2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right,
              n=n, p=p, x1=x1, x2=x2, x=x, var_list = var_list))
}


# ic, sieve
data_process_sieve <- function(data, l, u, var_list, m) {

  # replace Inf by constant u
  data$Right[data$status==0] <- u

  # split and calculate time ranges
  indata1 <- data[data[,"ind"]==1, ]
  indata2 <- data[data[,"ind"]==2, ]
  # indata1 <- subset(data, ind==1)
  # indata2 <- subset(data, ind==2)

  t1_left<-indata1[,"Left"]
  t1_right<-indata1[,"Right"]
  t2_left<-indata2[,"Left"]
  t2_right<-indata2[,"Right"]
  t_left <- data[,"Left"]
  t_right <- data[,"Right"]

  dim_m<-dim(as.matrix(indata1[,var_list]))
  n <- dim_m[1]

  # x1, x2, new var_list
  tmp1 <- var_transform(indata1, var_list)
  tmp2 <- var_transform(indata2, var_list)

  x1 <- tmp1$x
  x2 <- tmp2$x
  var_list <- tmp1$var_list


  # matrix
  x1 <- as.matrix(x1,dim_m[1])
  x2 <- as.matrix(x2,dim_m[1])
  p <- dim(x1)[2]


  # BP
  bl<-matrix(0,nrow = 2*n,ncol = m+1)
  br<-matrix(0,nrow = 2*n,ncol = m+1)
  for (i in 0:m) {
    bl[,(i+1)] <- bern(i,m,l,u,t_left)
    br[,(i+1)] <- bern(i,m,l,u,t_right)
  }
  odd_index <- seq(1,2*n,by=2)
  even_index <- seq(2,2*n,by=2)
  bl1 <- bl[odd_index,] #left eye, left end
  br1 <- br[odd_index,] #left eye, right end
  bl2 <- bl[even_index,] #right eye, left end
  br2 <- br[even_index,] #right eye, right end


  return(list(indata1=indata1, indata2=indata2, t1_left=t1_left, t1_right=t1_right, t2_left=t2_left, t2_right=t2_right,
              n=n, p=p, x1=x1, x2=x2, var_list = var_list, bl1=bl1, br1=br1, bl2=bl2, br2=br2))

}


# rc
get_covariates_rc <- function(data, var_list){

  # scan class of each covariate
  list_all <- split(var_list,sapply(data[1,var_list], function(x) paste(class(x), collapse=" ")))
  list_factor <- c(list_all$factor,list_all$character)
  list_num <- c(list_all$numeric,list_all$integer)

  # For data: extract numeric and factor/character variables
  x1_num <- data.frame(data[,list_num])
  colnames(x1_num) <- list_num
  x1_factor <- data.frame(data[,list_factor])
  colnames(x1_factor) <- list_factor
  if (ncol(x1_factor)!=0) {# generate dummy variables if x1_factor is not empty
    for (i in 1:ncol(x1_factor)){
      # for(level in unique(x1_factor[,i])[2:length(unique(x1_factor[,i]))]){ #use first level as reference
      tmp <- unique(x1_factor[,i])[order(unique(x1_factor[,i]))]
      for(level in tmp[2:length(tmp)]){ #use first level as reference
        x1_factor[paste(colnames(x1_factor)[i], level, sep = "")] <- ifelse(x1_factor[,i] == level, 1, 0)
      }
    }
    tmp <- colnames(x1_factor)[-(1:length(list_factor))]
    x1_factor <- data.frame(x1_factor[,-(1:length(list_factor))])
    colnames(x1_factor) <- tmp
    # x1_factor <- x1_factor[,-(1:length(list_factor))]
  }
  x1 <- cbind(x1_num,x1_factor) # generate new x1
  var_list <- colnames(x1) # definite new column names (with dummy variable names)

  return(list(x=x1,var_list=var_list))

}


data_preprocess_rc <- function(data, var_list){

  indata1 <- data[data[,"ind"]==1, ]
  indata2 <- data[data[,"ind"]==2, ]
  # indata1 <- subset(data, ind==1)
  # indata2 <- subset(data, ind==2)

  dim_m <- dim(as.matrix(indata1[,var_list]))
  n <- dim_m[1]

  return(list(indata1=indata1, indata2=indata2, n = n))

}


