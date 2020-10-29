estimate.flr2d <- function(object,tgrid,rgrid,sgrid,se.fit = FALSE)
{
  nnull <- length(object$d)
  nt <- length(tgrid)
  nr <- length(sgrid)
  ns <- length(sgrid)
  quad <- object$quad
  xdomain <- object$xdomain
  qpt <- (quad$pt-xdomain[1])/(xdomain[2]-xdomain[1])
  twk <- (tgrid-xdomain[1])/(xdomain[2]-xdomain[1])
  rwk <- (rgrid-xdomain[1])/(xdomain[2]-xdomain[1])
  swk <- (sgrid-xdomain[1])/(xdomain[2]-xdomain[1])
  Ktmat <- outer(twk,qpt,partK)
  Krmat <- outer(rwk,qpt,partK)
  Ksmat <- outer(swk,qpt,partK)
  Kymat <- Ktmat%*%(quad$wt*object$yqmat)
  Kxmat <- Krmat%*%(quad$wt*object$xqmat)
  Kzmat <- Ksmat%*%(quad$wt*object$zqmat)
  # use subset basis
  Kxmat <- Kxmat[,object$id.basis]
  Kzmat <- Kzmat[,object$id.basis]
  Kymat <- Kymat[,object$id.basis]
  q <- length(object$id.basis)
  #Sd
	btest=array(dim=c(nt,nr,ns))
	for(i in 1:nt)
    for(j in 1:nr)
      for(k in 1:ns)
		    btest[i,j,k] <-object$d[1]+object$d[2]*(rgrid[j]-0.5)+object$d[3]*(sgrid[k]-0.5)+object$d[4]*(sgrid[k]-0.5)*(rgrid[j]-0.5)+
                      object$d[5]*(tgrid[i]-0.5)+object$d[6]*(tgrid[i]-0.5)*(rgrid[j]-0.5)+object$d[7]*(tgrid[i]-0.5)*(sgrid[k]-0.5)+object$d[8]*(tgrid[i]-0.5)*(sgrid[k]-0.5)*(rgrid[j]-0.5)


  #Qc
  btestrk <- array(dim=c(nt,nr,ns))
  c=object$c*10^object$theta


  for(i in 1:nt)
    for(j in 1:nr)
      for(k in 1:ns){
        btestrk[i,j,k]=Kxmat[j,]%*%c[1:q]+Kxmat[j,]%*%c[(q+1):(2*q)]*(sgrid[k]-0.5)+(tgrid[i]-0.5)*Kxmat[j,]%*%c[(2*q+1):(3*q)]+
                      (tgrid[i]-0.5)*Kxmat[j,]%*%c[(3*q+1):(4*q)]*(sgrid[k]-0.5)+Kzmat[k,]%*%c[(4*q+1):(5*q)]+
                      Kzmat[k,]%*%c[(5*q+1):(6*q)]*(rgrid[j]-0.5)+(tgrid[i]-0.5)*Kzmat[k,]%*%c[(6*q+1):(7*q)]+
                      (tgrid[i]-0.5)*Kzmat[k,]%*%c[(7*q+1):(8*q)]*(rgrid[j]-0.5)+Kymat[i,]%*%c[(8*q+1):(9*q)]+
                      (rgrid[j]-0.5)*Kymat[i,]%*%c[(9*q+1):(10*q)]+(sgrid[k]-0.5)*Kymat[i,]%*%c[(10*q+1):(11*q)]+
                      (rgrid[j]-0.5)*(sgrid[k]-0.5)*Kymat[i,]%*%c[(11*q+1):(12*q)]+ ###first 12q
                      as.vector(outer(Kxmat[j,],Kzmat[k,]))%*%c[(12*q+1):(12*q+q^2)]+
                      (tgrid[i]-0.5)*as.vector(outer(Kxmat[j,],Kzmat[k,]))%*%c[(12*q+q^2+1):(12*q+2*q^2)]+
                      as.vector(outer(Kymat[i,],Kzmat[k,]))%*%c[(12*q+2*q^2+1):(12*q+3*q^2)]+
                      (rgrid[j]-0.5)*as.vector(outer(Kzmat[j,],Kymat[i,]))%*%c[(12*q+3*q^2+1):(12*q+4*q^2)]+
                      as.vector(outer(Kymat[i,],Kxmat[j,]))%*%c[(12*q+4*q^2+1):(12*q+5*q^2)]+
                      (sgrid[k]-0.5)*as.vector(outer(Kxmat[j,],Kymat[i,]))%*%c[(12*q+5*q^2+1):(12*q+6*q^2)]+ ###6n^2
                      outer(outer(Kymat[i,],Kxmat[j,]),Kzmat[k,])%*%c[(12*q+6*q^2+1):(12*q+6*q^2+q^3)]
      }


  if (se.fit) {
    b <- object$varht/10^object$nlambda
    pse <- NULL
  }


  btestsum <- btest+btestrk
  if (se.fit) list(fit = btest, se.fit = pse)
  else return(list(fit=btestsum,btest=btest,btestrk=btestrk))
}
