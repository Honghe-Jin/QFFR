library(gss)

partK <- function(x1,x2)
{# x1 and x2 vectors of same length
 # return value = k2(x1)*k2(x2)-k4(abs(x1 - x2)).
  val <- ((x1-0.5)^2-1/12)/2*((x2-0.5)^2-1/12)/2
  wk <- abs(x1-x2)
  val <- val-((wk-0.5)^4-(wk-0.5)^2/2+7/240)/24
}

perprod <- function(A,B){# return a matrix with columns (a1b1,a2b1,..,anb1,...,anbn)
	dA=dim(A)
	dB=dim(B)
	if(dA[1]!=dB[1]) error ("wrong dimensions")
	C=matrix(nrow=dA[1],ncol=dA[2]*dB[2])
	for (j in 1:dB[2])
		for (i in 1:dA[2])
			C[,i+(j-1)*(dA[2])]=A[,i]*B[,j]
	return(C)
}

flr2d <- function (yqmat,xqmat,zqmat,quad,method="v",alpha=1.4,id.basis=NULL,xdomain=c(0,1))
{## assume (r_1,...,r_T)=(t_1,...,t_T)=(s_1,...,s_T). yqmat,xqmat,zqmat\in R^(nobs*T)
  nobs <- dim(xqmat)[2]
  nmesh <- dim(xqmat)[1]
  if (is.null(id.basis)) id.basis <- 1:nobs
  nbasis <- length(id.basis)
  # need to center xqmat and yqmat first
  xmean <- apply(xqmat,1,mean)
	zmean <- apply(zqmat,1,mean)
  ymean <- apply(yqmat,1,mean)
  #xqmat <- xqmat - xmean
	#zqmat <- zqmat - zmean
  #yqmat <- yqmat - ymean
  # Kmat (nmesh*nmesh): Kmat[i,j] = K(t_i,t_j)
  # wmat (nmesh*nobs)*(nmesh*nobs): diagonal matrix with sqrt(a_1),...,sqrt(a_T)
  #    repeating nobs times on the diagonal.
	# Kxmat (nmesh*nobs): Kxmat[j,i]=int(K(t_j, .)*x_i(.))=(Kx_i)(t_j)
	# Kzmat (nmesh*nobs): Kzmat[j,i]=int(K(t_j, .)*z_i(.))=(Kz_i)(t_j)
  # Kymat (nmesh*nobs): Kymat[j,i]=int(K(t_j, .)*y_i(.))=(Ky_i)(t_j)
  qpt <- (quad$pt-xdomain[1])/(xdomain[2]-xdomain[1])
  wvec <- rep(sqrt(quad$wt),nobs)
  Kmat <- outer(qpt,qpt,partK)
  Kxmat <- Kmat%*%(quad$wt*xqmat)
  Kymat <- Kmat%*%(quad$wt*yqmat)
	Kzmat <- Kmat%*%(quad$wt*zqmat)
  # Before multiplying wvec,
	#####s
  # s matrix has row (i-1)*T+j as 1, intint(x_i(r)z_i(s)), intint(x_i(r)*(r-0.5)z_i(s)),
	#intint(x_i(r)*(s-0.5)z_i(s)),intint(x_i(r)*(r-0.5)z_i(s)(s-.5))
  #intint((t_j-0.5)x_i(r)z_i(s)), intint((t_j-0.5)x_i(r)*(r-0.5)z_i(s)),
	#intint((t_j-0.5)x_i(r)*(s-0.5)z_i(s)),intint((t_j-0.5)x_i(r)*(r-0.5)z_i(s)(s-.5))
  #####r
  # r matrix has row (i-1)*T+j as intint(Kx_1(r)*x_i(r)*z_i(s)), ..., intint(Kx_n(r)*x_i(r)*z_i(s));
	#   intint(Kx_1(r)*x_i(r)*z_i(s)(s-0.5)), ..., intint(Kx_n(r)*x_i(r)*z_i(s)(s-0.5));
	#   intint(Kx_1(r)*x_i(r)*z_i(s)(t_j-0.5)), ..., intint(Kx_n(r)*x_i(r)*z_i(s)(t_j-0.5));
	#   intint(Kx_1(r)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5)), ..., intint(Kx_n(r)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5));  ##(4n)

	#   intint(Kz_1(s)*x_i(r)*z_i(s)), ..., intint(Kz_n(s)*x_i(r)*z_i(s));
	#   intint(Kz_1(s)*x_i(r)*z_i(s)(r-0.5)), ..., intint(Kz_n(s)*x_i(r)*z_i(s)(r-0.5));
	#   intint(Kz_1(s)*x_i(r)*z_i(s)(t_j-0.5)), ..., intint(Kz_n(s)*x_i(r)*z_i(s)(t_j-0.5));
	#   intint(Kz_1(s)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5)), ..., intint(Kz_n(s)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5));  ##(4n)

	#   intint(Ky_1(t_j)*x_i(r)*z_i(s)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s));
	#   intint(Ky_1(t_j)*x_i(r)*z_i(s)(r-0.5)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s)(r-0.5));
	#   intint(Ky_1(t_j)*x_i(r)*z_i(s)(s-0.5)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s)(s-0.5));
	#   intint(Ky_1(t_j)*x_i(r)*z_i(s)(s-0.5)(r-0.5)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s)(s-0.5)(r-0.5));  ##(4n)

  #   intint(Kx_1(r)*x_i(r)*Kz_1(s)z_i(s)),...,intint(Kx_1(r)*x_i(r)*Kz_n(s)z_i(s)),...,intint(Kx_n(r)*x_i(r)*Kz_n(s)z_i(s)) #n^2
	#   (t_j-0.5)intint(Kx_1(r)*x_i(r)*Kz_1(s)z_i(s)),...,(t_j-0.5)intint(Kx_1(r)*x_i(r)*Kz_n(s)z_i(s)),...,(t_j-0.5)intint(Kx_n(r)*x_i(r)*Kz_n(s)z_i(s))
	#   intint(Ky_1(t_j)*x_i(r)*Kz_1(s)z_i(s)),...,intint(Ky_1(t_j)*x_i(r)*Kz_n(s)z_i(s)),...,intint(Ky_n(t_j)*x_i(r)*Kz_n(s)z_i(s))
	#   intint(Ky_1(t_j)*x_i(r)*Kz_1(s)z_i(s)(r-0.5)),...,intint(Ky_1(t_j)*x_i(r)*Kz_n(s)z_i(s)(r-0.5)),...,intint(Ky_n(t_j)*x_i(r)*Kz_n(s)z_i(s)(r-0.5))
	#   intint(Ky_1(t_j)*x_i(r)*Kx_1(s)z_i(s)),...,intint(Ky_1(t_j)*x_i(r)*Kx_n(s)z_i(s)),...,intint(Ky_n(t_j)*x_i(r)*Kx_n(s)z_i(s))
	#   intint(Ky_1(t_j)*x_i(r)*Kx_1(s)z_i(s)(s-0.5)),...,intint(Ky_1(t_j)*x_i(r)*Kx_n(s)z_i(s)(s-0.5)),...,intint(Ky_n(t_j)*x_i(r)*Kx_n(s)z_i(s)(s-0.5))

  #   intint(Ky_1(t_j)Kx_1(r)Kz_1(s)x_i(r)z_i(s)),...,intint(Ky_n(t_j)Kx_n(r)Kz_n(s)x_i(r)z_i(s)) #n^3

  #   dim=(nobs*nmesh)*(nobs^3+6*nobs^2+12*nobs)
  # q matrix is (nobs^3+6*nobs^2+12*nobs)*(nobs^3+6*nobs^2+12*nobs) block diagonal matrix
  #   diag(rxwk,rxwk,rywk,rywk,rxwk%x%rywk), %x% is Kronecker product.
  y <- wvec*c(yqmat)
	##############################################################################################################
	# swkx is nobs*2 with ith row int(x_i(r)), int(x_i(r)*(r-0.5))
	# swkz is nobs*2 with ith row int(z_i(s)), int(z_i(s)*(s-0.5))
  # swk is nobs*4 with ith row int(x_i(r)*z_i(s)), int(x_i(r)*(r-0.5)*z_i(s)),
	# int(x_i(r)*z_i(s)(s-0.5)),int(x_i(r)*(r-0.5)*z_i(s)(s-0.5))
	swkx <- cbind(apply(quad$wt*xqmat,2,sum),apply(quad$wt*(quad$pt-0.5)*xqmat,2,sum))
	swkz <- cbind(apply(quad$wt*zqmat,2,sum),apply(quad$wt*(quad$pt-0.5)*zqmat,2,sum))
  swk <- cbind(swkx[,1]*swkz[,1],swkx[,2]*swkz[,1],swkx[,1]*swkz[,2],swkx[,2]*swkz[,2])
	#swk=matrix(nrow=nobs,ncol4)
	#for (i in 1:nobs){
		#swk[i,1]=sum((quad$wt*xqmat[i,])%*%t(quad$wt*zqmat[i,]))
		#swk[i,2]=sum((quad$wt*(quad$pt-0.5)*xqmat[i,])%*%t(quad$wt*zqmat[i,]))
		#swk[i,3]=sum((quad$wt*xqmat[i,])%*%t(quad$wt*(quad$pt-0.5)*zqmat[i,]))
		#swk[i,4]=sum((quad$wt*(quad$pt-0.5)*xqmat[i,])%*%t(quad$wt*(quad$pt-0.5)*zqmat[i,]))
	#}

	##s <- cbind(swk%x%rep(1,nmesh),swk%x%(quad$pt-0.5)) This should be kept for any dimenional cases. Because this step is to do the tensor product for y
  s <- cbind(swk%x%rep(1,nmesh),swk%x%(quad$pt-0.5))


  # rxwk is nobs*nobs with i-th row int(x_1*(Kx_i)),...,int(x_n*(Kx_i))
	# rzwk is nobs*nobs with i-th row int(z_1*(Kz_i)),...,int(z_n*(Kz_i))

  rxwk <- t(Kxmat)%*%(quad$wt*xqmat)
	rzwk <- t(Kzmat)%*%(quad$wt*zqmat)

	#1:12n
  R1=cbind(swkz[,1]*rxwk[,id.basis],swkz[,2]*rxwk[,id.basis])
  R1=cbind(R1%x%rep(1,nmesh),R1%x%(quad$pt-0.5))

  R2=cbind(swkx[,1]*rzwk[,id.basis],swkx[,2]*rzwk[,id.basis])
  R2=cbind(R2%x%rep(1,nmesh),R2%x%(quad$pt-0.5))

  R3=perprod(swkx,swkz)
  R3=R3%x%Kymat[,id.basis]

  #12:12n+6n^2
  R4=perprod(rxwk[,id.basis],rzwk[,id.basis]) #yes
  R4=cbind(R4%x%rep(1,nmesh),R4%x%(quad$pt-0.5))

  R5=rzwk[,id.basis]%x%Kymat[,id.basis] #y first from 1 to n
  R5=cbind(R5*swkx[,1],R5*swkx[,2])

  R6=rxwk[,id.basis]%x%Kymat[,id.basis] #y first from 1 to n
  R6=cbind(R6*swkz[,1],R6*swkz[,2])

  #12:12n+6n^2
  R7=perprod(rxwk[,id.basis],rzwk[,id.basis])%x%Kymat[,id.basis] #y,x,z

	r=cbind(R1,R2,R3,R4,R5,R6,R7)

  s <- wvec*s
  r <- wvec*r
	##############################################################################################################
  q <- matrix(0,nbasis^3+6*nbasis^2+12*nbasis,nbasis^3+6*nbasis^2+12*nbasis)
  rywk <- t(Kymat)%*%(quad$wt*yqmat)
	#1:12n
  q[1:nbasis,1:nbasis] <- rxwk[id.basis,id.basis]
  q[(nbasis+1):(2*nbasis),(nbasis+1):(2*nbasis)] <- rxwk[id.basis,id.basis]/12
  q[(2*nbasis+1):(3*nbasis),(2*nbasis+1):(3*nbasis)] <- rxwk[id.basis,id.basis]/12
  q[(3*nbasis+1):(4*nbasis),(3*nbasis+1):(4*nbasis)] <- rxwk[id.basis,id.basis]/144
	q[(4*nbasis+1):(5*nbasis),(4*nbasis+1):(5*nbasis)] <- rzwk[id.basis,id.basis]
  q[(5*nbasis+1):(6*nbasis),(5*nbasis+1):(6*nbasis)] <- rzwk[id.basis,id.basis]/12
  q[(6*nbasis+1):(7*nbasis),(6*nbasis+1):(7*nbasis)] <- rzwk[id.basis,id.basis]/12
  q[(7*nbasis+1):(8*nbasis),(7*nbasis+1):(8*nbasis)] <-rzwk[id.basis,id.basis]/144
	q[(8*nbasis+1):(9*nbasis),(8*nbasis+1):(9*nbasis)] <- rywk[id.basis,id.basis]
	q[(9*nbasis+1):(10*nbasis),(9*nbasis+1):(10*nbasis)] <- rywk[id.basis,id.basis]/12
	q[(10*nbasis+1):(11*nbasis),(10*nbasis+1):(11*nbasis)] <- rywk[id.basis,id.basis]/12
	q[(11*nbasis+1):(12*nbasis),(11*nbasis+1):(12*nbasis)] <- rywk[id.basis,id.basis]/144
	#12n:12n+6n^2
  q[(12*nbasis+1):(12*nbasis+nbasis^2),(12*nbasis+1):(12*nbasis+nbasis^2)] <- rzwk[id.basis,id.basis]%x%rxwk[id.basis,id.basis]
	q[(12*nbasis+nbasis^2+1):(12*nbasis+2*nbasis^2),(12*nbasis+nbasis^2+1):(12*nbasis+2*nbasis^2)] <- rzwk[id.basis,id.basis]%x%rxwk[id.basis,id.basis]/12
	q[(12*nbasis+2*nbasis^2+1):(12*nbasis+3*nbasis^2),(12*nbasis+2*nbasis^2+1):(12*nbasis+3*nbasis^2)] <- rzwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
	q[(12*nbasis+3*nbasis^2+1):(12*nbasis+4*nbasis^2),(12*nbasis+3*nbasis^2+1):(12*nbasis+4*nbasis^2)] <- rzwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]/12
	q[(12*nbasis+4*nbasis^2+1):(12*nbasis+5*nbasis^2),(12*nbasis+4*nbasis^2+1):(12*nbasis+5*nbasis^2)] <- rxwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
	q[(12*nbasis+5*nbasis^2+1):(12*nbasis+6*nbasis^2),(12*nbasis+5*nbasis^2+1):(12*nbasis+6*nbasis^2)] <- rxwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]/12

	q[(12*nbasis+6*nbasis^2+1):(12*nbasis+6*nbasis^2+nbasis^3),(12*nbasis+6*nbasis^2+1):(12*nbasis+6*nbasis^2+nbasis^3)] <- rzwk[id.basis,id.basis]%x%rxwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
  #browser()
  wt <- NULL
  if (qr(s)$rank < dim(s)[2])
    stop("gss error in flr2d: fixed effects are linearly dependent")
  z <- mysspreg(s,r,q,y,wt,method,alpha,varht=1,random=NULL)
  # d = (d1, d2, d3, d4), c = (c1, ..., cn),
  # 10^(nlambda-theta)/nobs = lambda
  c(list(id.basis=id.basis,alpha=alpha,xdomain=xdomain,xmean=xmean,ymean=ymean,
    xqmat=xqmat,yqmat=yqmat,zqmat=zqmat,Kxmat=Kxmat,Kymat=Kymat,Kzmat=Kzmat,quad=quad),z)
}






#########################################################
regaux <- function(s,r,q,nlambda,fit)
{
  nnull <- dim(s)[2]
  nn <- nnull +  dim(q)[1]
  zzz <- eigen(q,symmetric=TRUE)
  rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
  val <- zzz$val[1:rkq]
  vec <- zzz$vec[,1:rkq,drop=FALSE]
  if (nnull) {
    wk1 <- qr(s)
    wk1 <- (qr.qty(wk1,r%*%vec))[-(1:nnull),]
  }
  else wk1 <- r%*%vec
  wk2 <- t(t(wk1)/sqrt(val))
  wk2 <- t(wk2)%*%wk2
  wk2 <- solve(wk2+diag(10^nlambda,dim(wk2)[1]),wk2)
  wk2 <- (wk2+t(wk2))/2
  wk2 <- t(wk2/sqrt(val))/sqrt(val)
  wk2 <- diag(1/val,dim(wk2)[1])-wk2
  z <- .Fortran("regaux",
                as.double(fit$chol), as.integer(nn),
                as.integer(fit$jpvt), as.integer(fit$rkv),
                drcr=as.double(t(cbind(s,r))%*%r%*%vec), as.integer(rkq),
                sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                PACKAGE="gss")[c("drcr","sms")]
  drcr <- matrix(z$drcr,nn,rkq)
  dr <- drcr[1:nnull,,drop=FALSE]
  sms <- 10^nlambda*matrix(z$sms,nnull,nnull)
  wk1 <- matrix(0,nnull+rkq,nnull+rkq)
  wk1[1:nnull,1:nnull] <- sms
  wk1[1:nnull,nnull+(1:rkq)] <- -t(t(dr)/val)
  wk1[nnull+(1:rkq),nnull+(1:rkq)] <- wk2
  z <- chol(wk1,pivot=TRUE)
  wk1 <- z
  rkw <- attr(z,"rank")
  while (wk1[rkw,rkw]<wk1[1,1]*sqrt(.Machine$double.eps)) rkw <- rkw-1
  wk1[row(wk1)>col(wk1)] <- 0
  if (rkw<nnull+rkq)
    wk1[(rkw+1):(nnull+rkq),(rkw+1):(nnull+rkq)] <- diag(0,nnull+rkq-rkw)
  hfac <- wk1
  hfac[,attr(z,"pivot")] <- wk1
  list(vec=vec,hfac=hfac)
}

mysspreg <- function(s,r,q,y,wt,method,alpha,varht,random)
{
  qr.trace <- FALSE
  if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
  alpha <- abs(alpha)
  ## get dimensions
  nobs <- nrow(r)
  nxi <- ncol(r)
  if (!is.null(s)) {
    if (is.vector(s)) nnull <- 1
    else nnull <- ncol(s)
  }
  else nnull <- 0
  if (!is.null(random)) nz <- ncol(as.matrix(random$z))
  else nz <- 0
  nxiz <- nxi + nz
  nn <- nxiz + nnull
  if (!is.null(wt)) {
    y <- wt*y
    s <- wt*s
    r <- wt*r
    if (!is.null(random)) random$z <- wt*random$z
  }
  ## cv function
  cv <- function(lambda) {
    if (is.null(random)) q.wk <- 10^(lambda+theta)*q
    else {
      q.wk <- matrix(0,nxiz,nxiz)
      q.wk[1:nxi,1:nxi] <- 10^(lambda[1]+theta)*q
      q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
        10^(2*ran.scal)*random$sigma$fun(lambda[-1],random$sigma$env)
    }
    if (qr.trace) {
      qq.wk <- chol(q.wk,pivot=TRUE)
      sr <- cbind(s,10^theta*r[,attr(qq.wk,"pivot")])
      sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
      sr <- qr(sr,tol=0)
      rss <- mean(qr.resid(sr,c(y,rep(0,nxiz)))[1:nobs]^2)
      trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
      if (method=="u") score <- rss + alpha*2*varht*trc
      if (method=="v") score <- rss/(1-alpha*trc)^2
      alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
      alpha.wk <- min(alpha.wk,3)
      if (alpha.wk>alpha) {
        if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
        if (method=="v") score <- rss/(1-alpha.wk*trc)^2
      }
      if (return.fit) {
        z <- .Fortran("reg",
                      as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxiz), as.double(y),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        z$score <- score
        assign("fit",z[c(1:5,7)],inherits=TRUE)
      }
    }
    else {
      z <- .Fortran("reg",
                    as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                    as.double(q.wk), as.integer(nxiz), as.double(y),
                    as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                    as.double(alpha), varht=as.double(varht),
                    score=double(1), dc=double(nn),
                    as.double(.Machine$double.eps),
                    chol=double(nn*nn), double(nn),
                    jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                    wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                    PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
      if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
      assign("fit",z[c(1:5,7)],inherits=TRUE)
      score <- z$score
      alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
      alpha.wk <- min(alpha.wk,3)
      if (alpha.wk>alpha) {
        if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
        if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
      }
    }
    score
  }
  cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
  ## initialization
  tmp <- sum(r^2)
  if (is.null(s)) theta <- 0
  else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
  log.la0 <- log10(tmp/sum(diag(q))) + theta
  if (!is.null(random)) {
    ran.scal <- theta - log10(sum(random$z^2)/nz/tmp*nxi) / 2
    r <- cbind(r,10^(ran.scal-theta)*random$z)
  }
  else ran.scal <- NULL
  ## lambda search
  return.fit <- FALSE
  fit <- NULL
  if (is.null(random)) la <- log.la0
  else la <- c(log.la0,random$init)
  if (length(la)-1) {
    counter <- 0
    ## scale and shift cv
    tmp <- abs(cv(la))
    cv.scale <- 1
    cv.shift <- 0
    if (tmp<1&tmp>10^(-4)) {
      cv.scale <- 10/tmp
      cv.shift <- 0
    }
    if (tmp<10^(-4)) {
      cv.scale <- 10^2
      cv.shift <- 10
    }
    repeat {
      zz <- nlm(cv.wk,la,stepmax=1,ndigit=7)
      if (zz$code<=3) break
      la <- zz$est
      counter <- counter + 1
      if (counter>=5) {
        warning("gss warning in ssanova: iteration for model selection fails to converge")
        break
      }
    }
  }
  else {
    repeat {
      mn <- la-1
      mx <- la+1
      if (mx>log.la0+6) break
      zz <- nlm0(cv,c(mn,mx))
      if (min(zz$est-mn,mx-zz$est)>=1e-3) break
      else la <- zz$est
    }
  }
  ## return
  return.fit <- TRUE
  jk1 <- cv(zz$est)
  if (is.null(random)) q.wk <- 10^theta*q
  else {
    q.wk <- matrix(0,nxiz,nxiz)
    q.wk[1:nxi,1:nxi] <- 10^theta*q
    q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
      10^(2*ran.scal-zz$est[1])*random$sigma$fun(zz$est[-1],random$sigma$env)
  }
  se.aux <- regaux(s,10^theta*r,q.wk,zz$est[1],fit)
  c <- fit$dc[nnull+(1:nxi)]
  if (nnull) d <- fit$dc[1:nnull]
  else d <- NULL
  if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
  else b <- NULL
  c(list(method=method,theta=theta,ran.scal=ran.scal,c=c,d=d,b=b,
         nlambda=zz$est[1],zeta=zz$est[-1]),fit[-3],list(se.aux=se.aux))
}

kerf<- function (yqmat,xqmat,zqmat,quad,method="v",alpha=1.4,id.basis=NULL,xdomain=c(0,1))
{## assume (r_1,...,r_T)=(t_1,...,t_T)=(s_1,...,s_T). yqmat,xqmat,zqmat\in R^(nobs*T)
  nobs <- dim(xqmat)[2]
  nmesh <- dim(xqmat)[1]
  if (is.null(id.basis)) id.basis <- 1:nobs
  nbasis <- length(id.basis)
  # need to center xqmat and yqmat first
  xmean <- apply(xqmat,1,mean)
  zmean <- apply(zqmat,1,mean)
  ymean <- apply(yqmat,1,mean)
  #xqmat <- xqmat - xmean
  #zqmat <- zqmat - zmean
  #yqmat <- yqmat - ymean
  # Kmat (nmesh*nmesh): Kmat[i,j] = K(t_i,t_j)
  # wmat (nmesh*nobs)*(nmesh*nobs): diagonal matrix with sqrt(a_1),...,sqrt(a_T)
  #    repeating nobs times on the diagonal.
  # Kxmat (nmesh*nobs): Kxmat[j,i]=int(K(t_j, .)*x_i(.))=(Kx_i)(t_j)
  # Kzmat (nmesh*nobs): Kzmat[j,i]=int(K(t_j, .)*z_i(.))=(Kz_i)(t_j)
  # Kymat (nmesh*nobs): Kymat[j,i]=int(K(t_j, .)*y_i(.))=(Ky_i)(t_j)
  qpt <- (quad$pt-xdomain[1])/(xdomain[2]-xdomain[1])
  wvec <- rep(sqrt(quad$wt),nobs)
  Kmat <- outer(qpt,qpt,partK)
  Kxmat <- Kmat%*%(quad$wt*xqmat)
  Kymat <- Kmat%*%(quad$wt*yqmat)
  Kzmat <- Kmat%*%(quad$wt*zqmat)
  # Before multiplying wvec,
  #####s
  # s matrix has row (i-1)*T+j as 1, intint(x_i(r)z_i(s)), intint(x_i(r)*(r-0.5)z_i(s)),
  #intint(x_i(r)*(s-0.5)z_i(s)),intint(x_i(r)*(r-0.5)z_i(s)(s-.5))
  #intint((t_j-0.5)x_i(r)z_i(s)), intint((t_j-0.5)x_i(r)*(r-0.5)z_i(s)),
  #intint((t_j-0.5)x_i(r)*(s-0.5)z_i(s)),intint((t_j-0.5)x_i(r)*(r-0.5)z_i(s)(s-.5))
  #####r
  # r matrix has row (i-1)*T+j as intint(Kx_1(r)*x_i(r)*z_i(s)), ..., intint(Kx_n(r)*x_i(r)*z_i(s));
  #   intint(Kx_1(r)*x_i(r)*z_i(s)(s-0.5)), ..., intint(Kx_n(r)*x_i(r)*z_i(s)(s-0.5));
  #   intint(Kx_1(r)*x_i(r)*z_i(s)(t_j-0.5)), ..., intint(Kx_n(r)*x_i(r)*z_i(s)(t_j-0.5));
  #   intint(Kx_1(r)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5)), ..., intint(Kx_n(r)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5));  ##(4n)
  
  #   intint(Kz_1(s)*x_i(r)*z_i(s)), ..., intint(Kz_n(s)*x_i(r)*z_i(s));
  #   intint(Kz_1(s)*x_i(r)*z_i(s)(r-0.5)), ..., intint(Kz_n(s)*x_i(r)*z_i(s)(r-0.5));
  #   intint(Kz_1(s)*x_i(r)*z_i(s)(t_j-0.5)), ..., intint(Kz_n(s)*x_i(r)*z_i(s)(t_j-0.5));
  #   intint(Kz_1(s)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5)), ..., intint(Kz_n(s)*x_i(r)*z_i(s)(s-0.5)(t_j-0.5));  ##(4n)
  
  #   intint(Ky_1(t_j)*x_i(r)*z_i(s)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s));
  #   intint(Ky_1(t_j)*x_i(r)*z_i(s)(r-0.5)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s)(r-0.5));
  #   intint(Ky_1(t_j)*x_i(r)*z_i(s)(s-0.5)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s)(s-0.5));
  #   intint(Ky_1(t_j)*x_i(r)*z_i(s)(s-0.5)(r-0.5)), ..., intint(Ky_n(t_j)*x_i(r)*z_i(s)(s-0.5)(r-0.5));  ##(4n)
  
  #   intint(Kx_1(r)*x_i(r)*Kz_1(s)z_i(s)),...,intint(Kx_1(r)*x_i(r)*Kz_n(s)z_i(s)),...,intint(Kx_n(r)*x_i(r)*Kz_n(s)z_i(s)) #n^2
  #   (t_j-0.5)intint(Kx_1(r)*x_i(r)*Kz_1(s)z_i(s)),...,(t_j-0.5)intint(Kx_1(r)*x_i(r)*Kz_n(s)z_i(s)),...,(t_j-0.5)intint(Kx_n(r)*x_i(r)*Kz_n(s)z_i(s))
  #   intint(Ky_1(t_j)*x_i(r)*Kz_1(s)z_i(s)),...,intint(Ky_1(t_j)*x_i(r)*Kz_n(s)z_i(s)),...,intint(Ky_n(t_j)*x_i(r)*Kz_n(s)z_i(s))
  #   intint(Ky_1(t_j)*x_i(r)*Kz_1(s)z_i(s)(r-0.5)),...,intint(Ky_1(t_j)*x_i(r)*Kz_n(s)z_i(s)(r-0.5)),...,intint(Ky_n(t_j)*x_i(r)*Kz_n(s)z_i(s)(r-0.5))
  #   intint(Ky_1(t_j)*x_i(r)*Kx_1(s)z_i(s)),...,intint(Ky_1(t_j)*x_i(r)*Kx_n(s)z_i(s)),...,intint(Ky_n(t_j)*x_i(r)*Kx_n(s)z_i(s))
  #   intint(Ky_1(t_j)*x_i(r)*Kx_1(s)z_i(s)(s-0.5)),...,intint(Ky_1(t_j)*x_i(r)*Kx_n(s)z_i(s)(s-0.5)),...,intint(Ky_n(t_j)*x_i(r)*Kx_n(s)z_i(s)(s-0.5))
  
  #   intint(Ky_1(t_j)Kx_1(r)Kz_1(s)x_i(r)z_i(s)),...,intint(Ky_n(t_j)Kx_n(r)Kz_n(s)x_i(r)z_i(s)) #n^3
  
  #   dim=(nobs*nmesh)*(nobs^3+6*nobs^2+12*nobs)
  # q matrix is (nobs^3+6*nobs^2+12*nobs)*(nobs^3+6*nobs^2+12*nobs) block diagonal matrix
  #   diag(rxwk,rxwk,rywk,rywk,rxwk%x%rywk), %x% is Kronecker product.
  y <- wvec*c(yqmat)
  ##############################################################################################################
  # swkx is nobs*2 with ith row int(x_i(r)), int(x_i(r)*(r-0.5))
  # swkz is nobs*2 with ith row int(z_i(s)), int(z_i(s)*(s-0.5))
  # swk is nobs*4 with ith row int(x_i(r)*z_i(s)), int(x_i(r)*(r-0.5)*z_i(s)),
  # int(x_i(r)*z_i(s)(s-0.5)),int(x_i(r)*(r-0.5)*z_i(s)(s-0.5))
  swkx <- cbind(apply(quad$wt*xqmat,2,sum),apply(quad$wt*(quad$pt-0.5)*xqmat,2,sum))
  swkz <- cbind(apply(quad$wt*zqmat,2,sum),apply(quad$wt*(quad$pt-0.5)*zqmat,2,sum))
  swk <- cbind(swkx[,1]*swkz[,1],swkx[,2]*swkz[,1],swkx[,1]*swkz[,2],swkx[,2]*swkz[,2])
  #swk=matrix(nrow=nobs,ncol4)
  #for (i in 1:nobs){
  #swk[i,1]=sum((quad$wt*xqmat[i,])%*%t(quad$wt*zqmat[i,]))
  #swk[i,2]=sum((quad$wt*(quad$pt-0.5)*xqmat[i,])%*%t(quad$wt*zqmat[i,]))
  #swk[i,3]=sum((quad$wt*xqmat[i,])%*%t(quad$wt*(quad$pt-0.5)*zqmat[i,]))
  #swk[i,4]=sum((quad$wt*(quad$pt-0.5)*xqmat[i,])%*%t(quad$wt*(quad$pt-0.5)*zqmat[i,]))
  #}
  
  ##s <- cbind(swk%x%rep(1,nmesh),swk%x%(quad$pt-0.5)) This should be kept for any dimenional cases. Because this step is to do the tensor product for y
  s <- cbind(swk%x%rep(1,nmesh),swk%x%(quad$pt-0.5))
  
  
  # rxwk is nobs*nobs with i-th row int(x_1*(Kx_i)),...,int(x_n*(Kx_i))
  # rzwk is nobs*nobs with i-th row int(z_1*(Kz_i)),...,int(z_n*(Kz_i))
  
  rxwk <- t(Kxmat)%*%(quad$wt*xqmat)
  rzwk <- t(Kzmat)%*%(quad$wt*zqmat)
  
  #1:12n
  R1=cbind(swkz[,1]*rxwk[,id.basis],swkz[,2]*rxwk[,id.basis])
  R1=cbind(R1%x%rep(1,nmesh),R1%x%(quad$pt-0.5))
  
  R2=cbind(swkx[,1]*rzwk[,id.basis],swkx[,2]*rzwk[,id.basis])
  R2=cbind(R2%x%rep(1,nmesh),R2%x%(quad$pt-0.5))
  
  R3=perprod(swkx,swkz)
  R3=R3%x%Kymat[,id.basis]
  
  #12:12n+6n^2
  R4=perprod(rxwk[,id.basis],rzwk[,id.basis]) #yes
  R4=cbind(R4%x%rep(1,nmesh),R4%x%(quad$pt-0.5))
  
  R5=rzwk[,id.basis]%x%Kymat[,id.basis] #y first from 1 to n
  R5=cbind(R5*swkx[,1],R5*swkx[,2])
  
  R6=rxwk[,id.basis]%x%Kymat[,id.basis] #y first from 1 to n
  R6=cbind(R6*swkz[,1],R6*swkz[,2])
  
  #12:12n+6n^2
  R7=perprod(rxwk[,id.basis],rzwk[,id.basis])%x%Kymat[,id.basis] #y,x,z
  
  r=cbind(R1,R2,R3,R4,R5,R6,R7)
  
  s <- wvec*s
  r <- wvec*r
  ##############################################################################################################
  q <- matrix(0,nbasis^3+6*nbasis^2+12*nbasis,nbasis^3+6*nbasis^2+12*nbasis)
  rywk <- t(Kymat)%*%(quad$wt*yqmat)
  #1:12n
  q[1:nbasis,1:nbasis] <- q[(nbasis+1):(2*nbasis),(nbasis+1):(2*nbasis)] <- rxwk[id.basis,id.basis]
  q[(2*nbasis+1):(3*nbasis),(2*nbasis+1):(3*nbasis)] <- q[(3*nbasis+1):(4*nbasis),(3*nbasis+1):(4*nbasis)] <- rxwk[id.basis,id.basis]
  q[(4*nbasis+1):(5*nbasis),(4*nbasis+1):(5*nbasis)] <- rzwk[id.basis,id.basis]
  q[(5*nbasis+1):(6*nbasis),(5*nbasis+1):(6*nbasis)] <- q[(6*nbasis+1):(7*nbasis),(6*nbasis+1):(7*nbasis)] <- rzwk[id.basis,id.basis]
  q[(7*nbasis+1):(8*nbasis),(7*nbasis+1):(8*nbasis)] <-rzwk[id.basis,id.basis]
  q[(8*nbasis+1):(9*nbasis),(8*nbasis+1):(9*nbasis)] <- rywk[id.basis,id.basis]
  q[(9*nbasis+1):(10*nbasis),(9*nbasis+1):(10*nbasis)] <- q[(10*nbasis+1):(11*nbasis),(10*nbasis+1):(11*nbasis)] <- rywk[id.basis,id.basis]
  q[(11*nbasis+1):(12*nbasis),(11*nbasis+1):(12*nbasis)] <- rywk[id.basis,id.basis]
  #12n:12n+6n^2
  q[(12*nbasis+1):(12*nbasis+nbasis^2),(12*nbasis+1):(12*nbasis+nbasis^2)] <- rzwk[id.basis,id.basis]%x%rxwk[id.basis,id.basis]
  q[(12*nbasis+nbasis^2+1):(12*nbasis+2*nbasis^2),(12*nbasis+nbasis^2+1):(12*nbasis+2*nbasis^2)] <- rzwk[id.basis,id.basis]%x%rxwk[id.basis,id.basis]
  q[(12*nbasis+2*nbasis^2+1):(12*nbasis+3*nbasis^2),(12*nbasis+2*nbasis^2+1):(12*nbasis+3*nbasis^2)] <- rzwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
  q[(12*nbasis+3*nbasis^2+1):(12*nbasis+4*nbasis^2),(12*nbasis+3*nbasis^2+1):(12*nbasis+4*nbasis^2)] <- rzwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
  q[(12*nbasis+4*nbasis^2+1):(12*nbasis+5*nbasis^2),(12*nbasis+4*nbasis^2+1):(12*nbasis+5*nbasis^2)] <- rxwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
  q[(12*nbasis+5*nbasis^2+1):(12*nbasis+6*nbasis^2),(12*nbasis+5*nbasis^2+1):(12*nbasis+6*nbasis^2)] <- rxwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
  
  q[(12*nbasis+6*nbasis^2+1):(12*nbasis+6*nbasis^2+nbasis^3),(12*nbasis+6*nbasis^2+1):(12*nbasis+6*nbasis^2+nbasis^3)] <- rzwk[id.basis,id.basis]%x%rxwk[id.basis,id.basis]%x%rywk[id.basis,id.basis]
  
  c(list(id.basis=id.basis,xdomain=xdomain,S=s,R=r,Q=q,y=y,
         xqmat=xqmat,yqmat=yqmat,Kxmat=Kxmat,Kymat=Kymat,quad=quad))
}
