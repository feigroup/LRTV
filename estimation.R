
library(fda)


############################################################
# Functions 
############################################################
# beta functions (ROI-wise, polynomial): generate beta_j(Z_ij) for each i,j
beta = function(coef,Z){
if(dim(coef)[1]!=dim(Z)[2]) stop("Dimension error.")
p0 = dim(coef)[2]
d = dim(coef)[1]
beta.Z0 = matrix(NA,dim(Z)[1],dim(Z)[2])
for(i in 1:d){
  beta.Z0[,i] = cbind(1,poly(Z[,i],p0-1,raw = T))%*%coef[i,]
}
beta.Z0
}

# estimate bG.true using B-spline regression
bG.true = function(Z,p,N,coef){
  beta.Z = beta(coef,Z)
  ## Use B-spline regression to estimate beta_j for each j
  do.call(rbind, lapply(1:dim(Z)[2], function(j){as.numeric(NA.adjust(coef(lm(beta.Z[,j]~bsplineS(Z[,j],breaks = seq(min(Z),max(Z),length.out = N),norder = p,nderiv = 0)-1))))} ))
}

# soft-thresholding SVD
SVT = function(D,c){
  singular.decomp = svd(D)
  U = singular.decomp$u
  V = singular.decomp$v
  sing.d = singular.decomp$d
  D_svt = tcrossprod(tcrossprod(U,diag(ifelse(sing.d>=c,sing.d-c,0))),V)
  return(D_svt)
}

# proximal gradient estimation
proximal = function(Y,Z,X,p,N,lambdaL,alpha,thres=1e-4,smooth.penalty=FALSE,lambdaP=0,D=NULL,verbose=FALSE){
  if(length(unique(list(dim(Y),dim(Z),dim(X))))!=1) stop("Y, Z, X have different dimensions.")
  if(smooth.penalty==FALSE) lambdaP = 0
  n = dim(Y)[1]
  d = dim(Y)[2]
  # Number of basis functions
  r = N + p - 2
  # Create B spline basis
  bbase = create.bspline.basis(c(min(Z), max(Z)), nbasis = r, norder = p)
  # Iterative estimation
  epsilon = 100
  bG = matrix(0,d,r)
  i = 0
  while(epsilon >= thres){
    time1 = Sys.time()
    # 1st part of nabla.f
    nabla1.1 = 2*Reduce('+', lapply(1:d, function(j){
      E = matrix(0,d,d)
      E[j,j] = 1
      E%*%bG%*%Reduce('+', lapply(1:n, function(i){t(eval.basis(bbase, Z[i,j])*X[i,j])%*%(eval.basis(bbase, Z[i,j])*X[i,j])}))
    }))
    nabla1.2 = -2*Reduce('+', lapply(1:d, function(j){
      ej = rep(0,d)
      ej[j] = 1
      ej%*%Reduce('+', lapply(1:n, function(i){eval.basis(bbase, Z[i,j])*X[i,j]*Y[i,j]}))
    }))
    nabla1 = nabla1.1 + nabla1.2
    # 2nd part of nabla.f
    nabla2 = if (smooth.penalty) lambdaP*(bG%*%(t(D)+D)) else 0
    nabla = nabla1 + nabla2
    # Soft-thresholding SVT
    bGnew = SVT(bG-alpha*nabla,alpha*lambdaL)
    epsilon = sum((bGnew-bG)^2)/sum(bG^2)
    time2 = Sys.time()
    t_temp = difftime(time2, time1, unit="s")
    i = i + 1
    if(verbose) print(paste("Iteration ",i,", epsilon = ",round(epsilon,5),", time = ",round(as.numeric(t_temp),2)," s.",sep=""))
    bG = bGnew
  }
  return(bG)
}

# replace NA by 0
NA.adjust = function(x){
  x[is.na(x)] = 0
  x
}

# return Y_hat from bG_hat
Y_hat = function(bG,Z,p,N){
  n = dim(Z)[1]
  d = dim(Z)[2]
  Y_hat = matrix(NA,n,d)
  for(i in 1:n){
    for(j in 1:d){
      Y_hat[i,j] = (bG%*%t(bsplineS(Z[i,j],breaks = seq(min(Z),max(Z),length.out = N),norder = p,nderiv = 0)*X[i,j]))[j,]
    }
  }
  Y_hat
}

# Auto-save regardless of inexistence of objects: there can be inexistent objects.
autosave <- function(...,file){
  .Internal(save(unlist(sapply(...,function(x) if(exists(x)) x )),file=file,ascii=F,version=NULL,envir=parent.frame(),eval.promises=T))
}

# Penalty matrix D
# eval.BBT = function(Z,N,p,x){t(bsplineS(x,breaks = seq(min(Z),max(Z),length.out = N),norder = p,nderiv = 2))%*%bsplineS(x,breaks = seq(min(Z),max(Z),length.out = N),norder = p,nderiv = 2)}
# integrate.matrix <- function(Z,N,p) {
# sdelta = (max(Z)-min(Z))/1e4
# slist = seq(min(Z),max(Z),length.out = 1e4+1)
# D = matrix(0,N + p - 2,N + p - 2)
# for (i in 1:(length(slist)-1)) {
#   D = D + (eval.BBT(Z,N,p,slist[i])+eval.BBT(Z,N,p,slist[i+1]))/2*sdelta
# }
# D
# }






############################################################
# Simulations
############################################################
# Fix Y,Z,X,p,coef
n = 50
d = 100
set.seed(1)
# Generate X in matrix form
X = matrix(rnorm(n*d),n,d)
# Generate Z in matrix form
Z = matrix(rnorm(n*d),n,d)
# Assume beta_j is a cubic function with coef[j,] as the coefficient
coef = matrix(runif(4*d,-2,2),d,4)

set.seed(2)
err = matrix(rnorm(n*d),n,d)
Y = beta(coef,Z)*X + err
p = 4


# Vary: N,lambdaL,alpha
N.list = 3:10
lambdaL.list = exp(seq(log(1), log(1000), length.out = 20))
alpha.list = exp(seq(log(1e-5), log(1e-2), length.out = 10))

for (i in 1:length(N.list)) {
  N = N.list[i]
  r = N + p -2
  temp = array(0,dim=c(d,r,length(lambdaL.list),length(alpha.list)))
  t.N1 = Sys.time()
  for (j in 1:length(lambdaL.list)) {
    lambdaL = lambdaL.list[j]
    t.lambda1 = Sys.time()
    for (k in 1:length(alpha.list)) {
      alpha = alpha.list[k]
      t.alpha1 = Sys.time()
      bG = proximal(Y,Z,X,p=4,N=N,lambdaL=lambdaL,alpha=alpha,thres=1e-4,verbose = F)
      t.alpha2 = Sys.time()
      temp[,,j,k] = bG
      cat(paste("N = ",N,", lambdaL = ",round(lambdaL,2),", alpha = ",formatC(alpha, format = "e", digits = 2),", time = ",round(as.numeric(difftime(t.alpha2,t.alpha1, unit="s")),2)," s.\n",sep=""))
    }
    t.lambda2 = Sys.time()
    cat(paste("Sub-iteration for N = ",N,", lambdaL = ",round(lambdaL,2), ", time = ",round(as.numeric(difftime(t.lambda2,t.lambda1, unit="min")),2)," min.\n\n",sep=""))
  }
  t.N2 = Sys.time()
  assign(paste("bG_hat_array.",N,sep=""),temp)
  autosave(paste("bG_hat_array.",N.list,sep=""),file="bG_hat.RData")
  cat(paste("Whole iteration for N = ",N, ", time = ",round(as.numeric(difftime(t.N2, t.N1, unit="h")),2)," h.\n\n\n",sep=""))
}










############################################################
# Testing functions 
############################################################

# n = 50
# d = 100
# 
# set.seed(1)
# # Generate X in matrix form
# X = matrix(rnorm(n*d),n,d)
# # Generate Z in matrix form
# Z = matrix(rnorm(n*d),n,d)
# # Assume beta_j is a cubic function with coef[j,] as the coefficient
# coef = matrix(runif(4*d,-2,2),d,4)
# 
# set.seed(2)
# err = matrix(rnorm(n*d),n,d)
# Y = beta(coef,Z)*X + err



#########################################
# Estimate bG from low-rank model
# Given Y,Z,X,p,N,lambdaL,alpha
# Test RMSE(\wh Y)

# t1=Sys.time()
# bG = proximal(Y,Z,X,p=4,N=10,lambdaL=1,alpha=5e-3,thres=1e-4,verbose = T)
# t2=Sys.time()
# t2-t1
# 
# bG1 = proximal(Y,Z,X,p=4,N=10,lambdaL=10,alpha=0.001,thres=1e-4,verbose = T)
# bG2 = proximal(Y,Z,X,p=4,N=10,lambdaL=500,alpha=5e-3,thres=1e-4,verbose = T)
# 
# Y_hat = Y_hat(bG,Z,p,N)
# sqrt(sum((Y_hat1-Y)^2/d/n))
# sqrt(sum((Y)^2/d/n))



#########################################
# Estimate bG from B-spline regression
# Given Z,p,N,coef
# Test low-rank assumption

# bG0 = bG.true(Z,p,N,coef)
# rankMatrix(bG0)




















