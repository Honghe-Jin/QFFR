library(gss)
source("funciton.R")
source("estimate.R")
xdomain <- c(0,1)
nquad <- 20 #(T)
xquad <- gauss.quad(nquad,xdomain)
SNR = 5
n=100

nu <- 2 


K <- 50
phi <- function(t,k) {#for x
  if (k==1) val <- rep(1,length(t))
  else val <- sqrt(2)*cos((k-1)*pi*t)
  val
}
iiphi <- function(t,k) {#for z
  if (k==1) val <- t*t/2
  else val <- sqrt(2)*(1-cos((k-1)*pi*t))/(((k-1)*pi)^.5)
  val
}


beta0 <- function(t,r,s) {
  exp(-t)*(r*s^2+r^2*cos(s*pi))
}

true=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      true[i,j,k]=beta0(xquad$pt[i],xquad$pt[j],xquad$pt[k])
persp(xquad$pt,xquad$pt,true[10,,],ticktype="detailed",xlab="s",ylab="r",zlab="",theta=300,phi=30,main = "true")


zeta <- sapply(1:K,function(k)(-1)^(k+1)*k^(-nu/2))
xmat <- zmat <- mu.true <- ymat <- NULL
Zmat <- matrix(runif(K*n,-sqrt(3),sqrt(3)),n,K)
sd.mu <- rep(0,n)
set.seed(15)
for (i in 1:n) {
  val <- iival <- 0
  for (k in 1:K) {
    val <- val + zeta[k]*Zmat[i,k]*phi(xquad$pt,k)
    iival <- iival + zeta[k]*Zmat[i,k]*iiphi(xquad$pt,k)
  }
  xmat <- cbind(xmat,val)
  zmat <- cbind(zmat,iival)
  
  muwk=rep(NA,nquad)
  for (nq in 1:nquad){
    sum=0
    for (nqq in 1:nquad)
      for(nqq2 in 1:nquad) sum=sum+xquad$wt[nqq]*val[nqq]*xquad$wt[nqq2]*iival[nqq2]*beta0(xquad$pt[nq],xquad$pt[nqq],xquad$pt[nqq2])
      muwk[nq]=sum
  }
  sd.mu[i] <- sd(muwk)
  mu.true <- cbind(mu.true,muwk)
}
ymat <- mu.true + matrix(rnorm(nquad*n,sd=mean(sd.mu)/sqrt(SNR)),ncol=n)

fl.fit=flr2d(ymat,xmat,zmat,quad=xquad,id.basis = 1:10)

est=estimate.flr2d(fl.fit,xquad$pt,xquad$pt,xquad$pt,se=F)
betahat=est$fit