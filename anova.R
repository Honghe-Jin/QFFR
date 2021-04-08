k1 <- function(x) x-.5
k2 <- function(x) x^2/2-x/2+1/12
k4 <- function(x) ((k1(x))^4-(k1(x))^2/2+7/240)/24
K2 <- function(x1,x2) k2(x1)*k2(x2)-k4(abs(x1-x2))

K2_2 <- function(x1,x2) k2(x2)-(x1-x2+.5*sign(x2-x1))^2/2+1/24 #d^2K2/dx1^2

############################################################
##integral
quad=gauss.quad(20,c(0,1))
f=function(x,y,z) x*y*z
wt3=outer(outer(quad$wt,quad$wt),quad$wt)

f_value <- array(dim = c(20,20,20))
for(i in 1:20)
  for(j in 1:20)
    for(k in 1:20) f_value[i,j,k] <- f(quad$pt[i],quad$pt[j],quad$pt[k])
sum(f_value*wt3)

wt2=outer(quad$wt,quad$wt)
sum(f_value[1,,]*wt2)

###########################################################
beta0 <- function(t,r,s){exp(-t)*(s^2*r+r^2*cos(s*pi)) }
true=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      true[i,j,k]=beta0(xquad$pt[i],xquad$pt[j],xquad$pt[k])

#true intercept
sum(true*wt3)
est0$beta_anova$intercept
#true x01 
p_beta_x <- function(t,r,s) { exp(-t)*(s^2+2*r*cos(pi*s))}
p_beta_x_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_x_value[i,j,k]=p_beta_x(xquad$pt[i],xquad$pt[j],xquad$pt[k])

#sum(p_beta_x_value*wt3)*(quad$pt-.5)

true_x01 <- sum(p_beta_x_value*wt3)*(quad$pt-.5)
plot(quad$pt,true_x01)
plot(quad$pt,est0$beta_anova$x01[1,,1])

#true z01
p_beta_z <- function(t,r,s) {exp(-t)*(2*r*s-pi*r^2*sin(pi*s))}
p_beta_z_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_z_value[i,j,k]=p_beta_z(xquad$pt[i],xquad$pt[j],xquad$pt[k])
#sum(p_beta_z_value*wt3)*(quad$pt-.5)
true_z01 <- sum(p_beta_z_value*wt3)*(quad$pt-.5)
plot(quad$pt,true_z01)
plot(quad$pt,est0$beta_anova$z01[1,1,])

#true y01
p_beta_y <- function(t,r,s) { -beta0(t,r,s)}
p_beta_y_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_y_value[i,j,k]=p_beta_y(xquad$pt[i],xquad$pt[j],xquad$pt[k])
true_y01 <- sum(p_beta_y_value*wt3)*(quad$pt-.5)
plot(quad$pt,true_y01)
plot(quad$pt,est0$beta_anova$y01[,1,1])

#true x01z01
p_beta_xz <- function(t,r,s) { exp(-t)*(2*s-2*pi*r*sin(pi*s))}
p_beta_xz_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xz_value[i,j,k]=p_beta_xz(xquad$pt[i],xquad$pt[j],xquad$pt[k])
x01z01 <- sum(p_beta_xz_value*wt3)*outer(quad$pt-.5,quad$pt-.5)
persp(xquad$pt,xquad$pt,x01z01,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$x01z01[1,,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

#true x01y01
p_beta_xy <- function(t,r,s) { -p_beta_x(t,r,s)}
p_beta_xy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xy_value[i,j,k]=p_beta_xy(xquad$pt[i],xquad$pt[j],xquad$pt[k])
x01y01 <- sum(p_beta_xy_value*wt3)*outer(quad$pt-.5,quad$pt-.5)
persp(xquad$pt,xquad$pt,x01y01,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$x01y01[,,1],ticktype="detailed",xlab="r",ylab="s",zlab="beta",theta=300,phi=30)

#true z01y01
p_beta_zy <- function(t,r,s) { -p_beta_z(t,r,s)}
p_beta_zy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_zy_value[i,j,k]=p_beta_zy(xquad$pt[i],xquad$pt[j],xquad$pt[k])
z01y01 <- sum(p_beta_zy_value*wt3)*outer(quad$pt-.5,quad$pt-.5)
persp(xquad$pt,xquad$pt,z01y01 ,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$z01y01[,1,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true x01z01y01
p_beta_xzy <- function(t,r,s) { 0}
p_beta_xzy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xzy_value[i,j,k]=p_beta_xzy(xquad$pt[i],xquad$pt[j],xquad$pt[k])
#sum(p_beta_xyz_value*wt3)*outer(outer(quad$pt-.5,quad$pt-.5),quad$pt-.5)

# true x1
p_beta_xx <- function(t,r,s) { 2*exp(-t)*cos(pi*s) }
p_beta_xx_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xx_value[i,j,k]=p_beta_xx(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Kx2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(x1,x2)/dx1^2
int_x <- numeric(20) #intint p_beta_xx dzdy
for(i in 1:20) int_x[i] <- sum(p_beta_xx_value[,i,]*wt2)
true_x1 <- colSums(Kx2*int_x*quad$wt)
plot(xquad$pt,true_x1)
plot(xquad$pt,est0$beta_anova$x1[1,,1])

# true z1
p_beta_zz <- function(t,r,s) { 2*r*exp(-t)-pi^2*r^2*cos(pi*s)*exp(-t) }
p_beta_zz_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_zz_value[i,j,k]=p_beta_zz(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Kz2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(z1,z2)/dz1^2
int_z <- numeric(20) #intint p_beta_xx dxdy
for(i in 1:20) int_z[i] <- sum(p_beta_zz_value[,,i]*wt2)
true_z1 <- colSums(Kz2*int_z*quad$wt)
plot(xquad$pt,true_z1)
plot(xquad$pt,est0$beta_anova$z1[1,1,])

# true y1
p_beta_yy <- function(t,r,s) beta0(t,r,s)
p_beta_yy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_yy_value[i,j,k]=p_beta_yy(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Ky2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(y1,y2)/dy1^2
int_y <- numeric(20) #intint p_beta_xx dxdz
for(i in 1:20) int_y[i] <- sum(p_beta_yy_value[i,,]*wt2)
true_y1 <- colSums(Ky2*int_y*quad$wt)
plot(xquad$pt,true_y1)
plot(xquad$pt,est0$beta_anova$y1[,1,1])

# true x1z01
p_beta_xxz <- function(t,r,s) -2*pi*exp(-t)*sin(pi*s)
p_beta_xxz_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xxz_value[i,j,k]=p_beta_xxz(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Kx2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(y1,y2)/dy1^2
int_xxz <- numeric(20) #intint p_beta_xxz dzdy
for(i in 1:20) int_xxz[i] <- sum(p_beta_xxz_value[,i,]*wt2)

x1z01 <- outer(colSums(Kx2*int_xxz),quad$pt-.5)
persp(xquad$pt,xquad$pt,outer(colSums(Kx2*int_xxz*quad$wt),quad$pt-.5),ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$x1z01[1,,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true x1y01
p_beta_xxy <- function(t,r,s) -2*exp(-t)*cos(pi*s)
p_beta_xxy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xxy_value[i,j,k]=p_beta_xxy(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Kx2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(y1,y2)/dy1^2
int_xxy <- numeric(20) #intint p_beta_xxz dzdy
for(i in 1:20) int_xxy[i] <- sum(p_beta_xxy_value[,i,]*wt2)

x1y01 <- outer(quad$pt-.5,colSums(Kx2*int_xxy*quad$wt))
persp(xquad$pt,xquad$pt,x1y01,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$x1y01[,,1],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true z1y01
p_beta_zzy <- function(t,r,s) -2*r*exp(-t) + pi^2*r^2*cos(pi*s)*exp(-t)
p_beta_zzy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_zzy_value[i,j,k]=p_beta_zzy(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Kz2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(y1,y2)/dy1^2
int_zzy <- numeric(20) #intint p_beta_xxz dzdy
for(i in 1:20) int_zzy[i] <- sum(p_beta_zzy_value[,,i]*wt2)

z1y01 <- outer(quad$pt-.5,colSums(Kx2*int_zzy*quad$wt))
persp(xquad$pt,xquad$pt,z1y01,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$z1y01[,1,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true x01y1
p_beta_xyy <- function(t,r,s) p_beta_x(t,r,s)
p_beta_xyy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xyy_value[i,j,k]=p_beta_xyy(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Ky2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(y1,y2)/dy1^2
int_xyy <- numeric(20) #intint p_beta_xxz dzdy
for(i in 1:20) int_xyy[i] <- sum(p_beta_xyy_value[i,,]*wt2)

x01y1 <- outer(colSums(Ky2*int_xyy*quad$wt),quad$pt-.5)
persp(xquad$pt,xquad$pt,x01y1,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$x01y1[,,1],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true z01y1
p_beta_zyy <- function(t,r,s) p_beta_z(t,r,s)
p_beta_zyy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_zyy_value[i,j,k]=p_beta_zyy(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Ky2 <- outer(quad$pt,quad$pt,K2_2) #d^2K2(y1,y2)/dy1^2
int_zyy <- numeric(20) #intint p_beta_xxz dzdy
for(i in 1:20) int_zyy[i] <- sum(p_beta_zyy_value[i,,]*wt2)

z01y1 <- outer(colSums(Ky2*int_zyy*quad$wt),quad$pt-.5)
persp(xquad$pt,xquad$pt,z01y1,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$z01y1[,1,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true x1y1
p_beta_xxyy <- function(t,r,s) p_beta_xx(t,r,s)
p_beta_xxyy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_xxyy_value[i,j,k]=p_beta_xxyy(xquad$pt[i],xquad$pt[j],xquad$pt[k])

Kx2 <- outer(quad$pt,quad$pt,K2_2)
Ky2 <- outer(quad$pt,quad$pt,K2_2)
Kx2y2 <- outer(Kx2,Ky2) #(i,j,k,l)=Kx[i,j]*Ky[k,l]

int_xy <- matrix(nrow=20,ncol=20) #int p_beta_xxyy dz
for(i in 1:20)
  for(j in 1:20){
    int_xy[i,j] <- sum(quad$wt*p_beta_xxyy_value[i,j,])
  }

x1y1 <- matrix(nrow=20,ncol=20)
for(i in 1:20)
  for(j in 1:20){
    x1y1[i,j] <- sum(wt2*int_xy*Kx2y2[j,,i,])
  }
persp(xquad$pt,xquad$pt,x1y1,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(xquad$pt,xquad$pt,est0$beta_anova$x1y1[,,1],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

# true z1y1
p_beta_zzyy <- function(t,r,s) p_beta_zz(t,r,s)
p_beta_zzyy_value=array(dim=c(20,20,20))
for (i in 1:20)
  for(j in 1:20)
    for(k in 1:20)
      p_beta_zzyy_value[i,j,k]=p_beta_zzyy(quad$pt[i],quad$pt[j],quad$pt[k])

Kz2 <- outer(quad$pt,quad$pt,K2_2)
Ky2 <- outer(quad$pt,quad$pt,K2_2)
Kz2y2 <- outer(Kz2,Ky2) #(i,j,k,l)=Kz[i,j]*Ky[k,l]

int_zy <- matrix(nrow=20,ncol=20) #int p_beta_zzyy dx
for(i in 1:20)
  for(j in 1:20){
    int_zy[i,j] <- sum(quad$wt*p_beta_zzyy_value[i,,j])
  }

z1y1 <- matrix(nrow=20,ncol=20)
for(i in 1:20)
  for(j in 1:20){
    z1y1[i,j] <- sum(wt2*int_zy*Kz2y2[i,,j,])
  }
persp(quad$pt,quad$pt,z1y1,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(quad$pt,quad$pt,est0$beta_anova$z1y1[,1,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)









##############################################################################
#XY <- x01y01+x1y01+x1y1+x01y1
#Xmain <- x01+x1

XY <- x01y01+x1y01+x1y1+x01y1
persp(quad$pt,quad$pt,XY,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(quad$pt,quad$pt,est0$beta_anova$XY[,,1],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)


ZY <- z01y01+z01y1+z1y01+z1y1
persp(quad$pt,quad$pt,ZY,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(quad$pt,quad$pt,est0$beta_anova$ZY[,1,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)


Xmain <- true_x01+true_x1
Ymain <- true_y01+true_y1
Zmain <- true_z01+true_z1

plot(quad$pt,Xmain)
plot(quad$pt,est0$beta_anova$Xmain[1,,1])
plot(quad$pt,Zmain)
plot(quad$pt,est0$beta_anova$Zmain[1,1,])
plot(quad$pt,Ymain)
plot(quad$pt,est0$beta_anova$Ymain[,1,1])

persp(quad$pt,quad$pt,XY+Xmain,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(quad$pt,quad$pt,est0$beta_anova$XY[,,1]+est0$beta_anova$Xmain[1,,1],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)

persp(quad$pt,quad$pt,ZY+Zmain,ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)
persp(quad$pt,quad$pt,est0$beta_anova$ZY[,1,]+est0$beta_anova$Zmain[1,1,],ticktype="detailed",xlab="r",ylab="s",zlab="",theta=300,phi=30)


pdf("~/Documents/researches/function on function/anova/beta_xy_fit.pdf")
persp(quad$pt,quad$pt,est0$beta_anova$XY[,,1]+est0$beta_anova$Xmain[1,,1],ticktype="detailed",xlab="",ylab="",zlab="",theta=300,phi=30)
dev.off()

pdf("~/Documents/researches/function on function/anova/beta_xy_true.pdf")
persp(quad$pt,quad$pt,XY+Xmain,ticktype="detailed",xlab="",ylab="",zlab="",theta=300,phi=30)
dev.off()

pdf("~/Documents/researches/function on function/anova/beta_zy_fit.pdf")
persp(quad$pt,quad$pt,est0$beta_anova$ZY[,1,]+est0$beta_anova$Zmain[1,1,],ticktype="detailed",xlab="",ylab="",zlab="",theta=300,phi=30)
dev.off()

pdf("~/Documents/researches/function on function/anova/beta_zy_true.pdf")
persp(quad$pt,quad$pt,ZY+Zmain,ticktype="detailed",xlab="",ylab="",zlab="",theta=300,phi=30)
dev.off()

