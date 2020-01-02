## -----------------------------------------------------------------------------
sseR=function(n){
f=function(x){x^n*exp(-x)*(x>0)*(x<1)}
s=numeric(4)
for(i in 1:4){
  s[i]=mean(f(runif(4,(i-1)*0.25,i*0.25)))
}
mean(s)
}
sseR(0);sseR(1);sseR(1.5)

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('double sseC(double n){
  NumericVector s(4);
  for(int i=0;i<4;i++){
  NumericVector x=runif(4,i*0.25,(i+1)*0.25);
  NumericVector t(4);
  for(int j=0;j<4;j++){
  t[j]=pow(x[j],n)*exp(-x[j]);
  }
  s[i]=mean(t);}
return mean(s);
  
}')
sseC(0);sseC(1);sseC(1.5)

## -----------------------------------------------------------------------------
library(microbenchmark)
tm2 <- microbenchmark(
  vR = sseR(5),
  vC = sseC(5)
)
knitr::kable(summary(tm2)[,c(1,3,5,6)])

## -----------------------------------------------------------------------------
xtR=function(N,lambda){
x=rexp(N,lambda)
Tn=vector(length=N)
Tn[1]=x[1]
for(i in 2:N){
  Tn[i]=Tn[i-1]+x[i]
} #the total time
t=10
i=1
for(i in 1:N){
  if(Tn[i]<=t&&Tn[i+1]>t){
    xt=i
  }
}
xt
}
t=numeric(1e3)
for(i in 1:1e3){
  t[i]=xtR(100,1)
}
plot(t,ylab = "xt")
round(mean(t))

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction(' NumericVector xtC(int N,double lambda){
NumericVector x=rexp(N,lambda);
NumericVector Tn(N);
NumericVector xt(1);
Tn[0]=x[0];
for(int i=1;i<N;i++){
  Tn[i]=Tn[i-1]+x[i];
} 
for(int i=0;i<N;i++){
 int t=10;
  if(Tn[i]<=t&&Tn[i+1]>t){
   xt=i;}}
return xt;
}')
t=numeric(1e3)
for(i in 1:1e3){
  t[i]=xtC(100,1)
}
plot(t,ylab = "xt")
round(mean(t))

## -----------------------------------------------------------------------------
tm2 <- microbenchmark(
  vR = xtR(100, 1),
  vC = xtC(100, 1)
)
knitr::kable(summary(tm2)[,c(1,3,5,6)])

## -----------------------------------------------------------------------------
library(e1071) 
Brown=function(t,n){
N=1000 
m=n-1
time=seq(0,t,length=N) 
path=rwiener(end=1,frequency=N) 
plot(time,path,xlab="time",type="l",main=" Standard BM")
for(i in 1:m){ 
  path=rwiener(end=1,frequency=N) 
  lines(time,path)
}
}
Brown(1,6)

## -----------------------------------------------------------------------------
# 3.4
n=1000
u=runif(n)

Sigma1=.05
x=(-2*(log(1-u))*Sigma1^2)^(1/2)
hist(x,prob=TRUE,main = "sigma=0.05")
y=seq(0,1,.01)
lines(y,(y/Sigma1^2)*exp(-y^2/(2*(Sigma1^2))))

Sigma2=.1
x=(-2*(log(1-u))*Sigma2^2)^(1/2)
hist(x,prob=TRUE,main = "sigma=0.1")
y=seq(0,1,.01)
lines(y,(y/Sigma2^2)*exp(-y^2/(2*(Sigma2^2))))

Sigma3=.2
x=(-2*(log(1-u))*Sigma3^2)^(1/2)
hist(x,prob=TRUE,main = "sigma=0.2")
y=seq(0,1,.01)
lines(y,(y/Sigma3^2)*exp(-y^2/(2*(Sigma3^2))))

Sigma4=.5
x=(-2*(log(1-u))*Sigma4^2)^(1/2)
hist(x,prob=TRUE,main = "sigma=0.5")
y=seq(0,1,.01)
lines(y,(y/Sigma4^2)*exp(-y^2/(2*(Sigma4^2))))
# 3.11
n=1000
x1=rnorm(n,0,1)
x2=rnorm(n,3,1)
u=runif(n)
p=runif(1)
k=as.integer(u < p)
y=k * x1 + (1-k) * x2 
hist(y,freq = FALSE,main = "y")
lines(density(y))

u1=runif(n)
k1=as.integer(u1<.75)
y1=k1* x1 + (1-k1) * x2 
hist(y1,freq = FALSE,main = "p=0.75")
lines(density(y1),lwd=2)
u2=runif(n)
k2=as.integer(u1<.6)
y2=k2* x1 + (1-k2) * x2 
hist(y2,freq = FALSE,main = "p=0.6")
lines(density(y2),lwd=2)
# 3.18
X=function(n,d,Sigma){
  mu=rep(0,d)
  ev=eigen(Sigma,symmetric = TRUE)
  lambda = ev$values
  V = ev$vectors
  C = V %*% diag(sqrt(lambda)) %*% t(V)
  Z = matrix(rnorm(n*d), nrow = n, ncol = d)
  X = Z %*% C + matrix(mu, n, d, byrow = TRUE)
  X
}
t=X(3,3,matrix(c(.1,.2,.3,.4,.5,.6,.7,.8,.9),nrow = 3,ncol = 3))
print(t)
cov(t)

## -----------------------------------------------------------------------------
# 5.1
n=1e3
m=1e5
t1=runif(n,min = 0,max = pi*(1/3))
t2=runif(m,min = 0,max = pi*(1/3))
theta.hat1=mean(sin(t1))*pi*(1/3)
theta.hat2=mean(sin(t2))*pi*(1/3)
print(c(theta.hat1,theta.hat2,cos(0)-cos(pi*(1/3))))
# 5.10
#Monte Carlo estimate
n=1e4
x=runif(n,min = 0,max = 1)
theta.hat=mean(exp(-x)/(1+x^2))
t=integrate(function(x)exp(-x)/(1+x^2),0,1)
theta.hat;t
#antithetic variables
MC.Phi=function(x,n=1e4,antithetic=TRUE){
  u=runif(n/2)
  if(!antithetic)v=runif(n/2)else
    v=1-u
  u=c(u,v)
  cdf=numeric(length(x))
  for(i in 1:length(x)){
    g=x[i]*exp(-u*x[i])/(1+(u*x[i])^2)
    cdf[i]=mean(g)
  }
  cdf
}
x=1
theta.hat2=integrate(function(y)exp(-y)/(1+y^2),0,x)
Phi=theta.hat2
MC1=MC2=numeric(n)
for(i in 1:n){
MC1[i]=MC.Phi(x,antithetic = FALSE)
MC2[i]=MC.Phi(x)}
print(c(mean(MC1),mean(MC2)))
c(sd(MC1),sd(MC2),(sd(MC1)-sd(MC2))/sd(MC1))
# 5.15
m =10000
theta.hat =se =numeric(5)
g =function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
x=runif(m)
fg=g(x)
theta.hat[1]= mean(fg)
se[1]=sd(fg)
x=rexp(m, 1)
fg=g(x) / exp(-x)
theta.hat[2]=mean(fg)
se[2]=sd(fg)
x=rcauchy(m)
i=c(which(x > 1), which(x < 0))
x[i]=2 
fg=g(x) / dcauchy(x)
theta.hat[3]=mean(fg)
se[3]=sd(fg)
u =runif(m) 
x=- log(1 - u * (1 - exp(-1)))
fg=g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[4]=mean(fg)
se[4]=sd(fg)
u=runif(m) 
x=tan(pi * u / 4)
fg=g(x) / (4 / ((1 + x^2) * pi))
theta.hat[5]=mean(fg)
se[5]=sd(fg)
rbind(theta.hat,se)
#stratified importance sampling
m=1e4
u=runif(m)
k=5 
r=m/k
n=50
t=numeric(k)
est=matrix(0, n, 2)
x=-log(1-u*(1-exp(-1)))
g=function(x){exp(-x)/(1+x^2)*(x>0)*(x<1)}
for (i in 1:n) {
est[i, 1]=mean(g(runif(m)))
for(j in 1:k)t[j]=mean(g(runif(m/k,(j-1)/k,j/k)))
est[i, 2]=mean(t)
}
apply(est,2,mean)
apply(est,2,var)


## -----------------------------------------------------------------------------
#6.5
n=20
alpha=.05
theta.hat1=replicate(10000,expr ={
  x=rchisq(n,2)
  y=sum(x)
  sigma=var(x)^(1/2)
  theta.hat=mean(x)-(y^(1/2)*sigma*qt(alpha,df=40))/(800^(1/2))
})
print(c(mean(theta.hat1>=2)),4)
theta.hat2=replicate(10000,expr ={
  x=rchisq(n,2)
  s=var(x)^(1/2)
  xbat=mean(x)
  xbat+abs((s*qt(alpha,df=n-1))/n^(1/2))
})
print(c(mean(theta.hat2>=2)),4)
ucl=replicate(10000,expr ={
  x=rchisq(n,2)
  (n-1)*var(x)/qchisq(alpha,df=n-1)
})
mean(ucl>4)

# 6.6
m=1e4
n=1e3
p=c(0.025,0.05,0.95,0.975)
#skewness under normality
ske=function(x){ 
  x=rnorm(n)
  x.bat=mean(x)
  m1=mean((x-x.bat)^3)
  m2=mean((x-x.bat)^2)
  return (m1/m2^1.5)
}
k=replicate(m,ske(x))
c1=quantile(k,p)
c1
z=(p*(1-p))/(n*dnorm(c1,0,sqrt(6/n))^2)
sqrt(z)
ske.hat=rnorm(m,0,sqrt(6/n))
c2=quantile(ske.hat,p)
c2
c3=qnorm(p,0,sqrt(6/n))
c=data.frame(c1,c2,c3,z)
c

## -----------------------------------------------------------------------------
# 6.7
alpha=0.05
n=20
m=1e4
d=c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N=length(d)
pwr=pwr1=numeric(N)
cv=qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sk=function(x){ 
  x.bat=mean(x)
  m1=mean((x-x.bat)^3)
  m2=mean((x-x.bat)^2)
  return (m1/m2^1.5)
}
for(j in 1:N){
  e=d[j]
  sktests=numeric(m)
  for(i in 1:m){
  sigma=sample(c(1,100),replace = TRUE,n,c(1-e,e))
  x=rbeta(n,sigma,sigma)
  sktests[i]=as.integer(abs(sk(x))>=cv)
  }
 pwr[j]=mean(sktests)
}
mean(pwr)
for(j in 1:N){
  e=d[j]
  sktests=numeric(m)
  for(i in 1:m){
  sigma=sample(c(1,10),replace = TRUE,n,c(1-e,e))
  x=rt(n,df=sigma)
  sktests[i]=as.integer(abs(sk(x))>=cv)
  }
  pwr1[j]=mean(sktests)
}
mean(pwr1)
par(mfrow=c(1,2))
plot(d,pwr,type = "b",xlab = bquote(d),ylim = c(0,1),pch="￥",col="red3")
se=sqrt(pwr*(1-pwr)/m)
lines(d,pwr+se,col="green",lwd=1)
lines(d,pwr-se,col="blue",lwd=1)
plot(d,pwr1,type = "b",xlab = bquote(d),ylim = c(0,1),pch="￥",col="red3")
se1=sqrt(pwr1*(1-pwr1)/m)
lines(d,pwr1+se1,col="green",lwd=1)
lines(d,pwr1-se1,col="blue",lwd=1)
# 6.A
alpha=0.05
mu0=1
n=20
m=1e4
p1=p2=p3=numeric(m)
for(j in 1:m){
  x=rchisq(n,1)
  y=runif(n,0,2)
  z=rexp(n,1)
  t1=t.test(x,alternative = "greater",mu=mu0)
  t2=t.test(y,alternative = "greater",mu=mu0)
  t3=t.test(z,alternative = "greater",mu=mu0)
  p1[j]=t1$p.value
  p2[j]=t2$p.value
  p3[j]=t3$p.value
}
p.hat=c(mean(p1<alpha),mean(p2<alpha),mean(p3<alpha))
se.hat=sqrt(p.hat*(1-p.hat)/m)
data.frame(p.hat,se.hat,row.names = c("(i)χ2(1)","(ii)Uniform(0,2)","(iii)Exponential(1)"))

## ----warning=FALSE------------------------------------------------------------
# 7.6
library(bootstrap)
pairs(scor[1:5],pch=20,col="blue")
COR=cor(scor)
COR
N=1e4
n=88
r1=r2=r3=r4=numeric(N)
for(j in 1:N){
  i=sample(1:n,n,replace = TRUE)
  x1=scor$mec[i]
  x2=scor$vec[i]
  x3=scor$alg[i]
  x4=scor$ana[i]
  x5=scor$sta[i]
  r1[j]=cor(x1,x2)
  r2[j]=cor(x3,x4)
  r3[j]=cor(x3,x5)
  r4[j]=cor(x4,x5)
}
s.e=c(se.r1=sd(r1),se.r2=sd(r2),se.r3=sd(r3),se.r4=sd(r4))
print(s.e,4)
# 7.B
m=1000
n=100
mu=0
t=numeric(n)
set.seed(1234)
library(boot)
ci.norm=ci.basic=ci.perc=matrix(NA,m,2)
sk=function(x){ 
  x.bat=mean(x)
  m1=mean((x-x.bat)^3)
  m2=mean((x-x.bat)^2)
  return (m1/m2^1.5)
}
boot.median=function(x,j) mean(sk(x[j]))
# normal populations
for(i in 1:m){
for(j in 1:n){
  x=rnorm(n,0,1)
  t[j]=sk(x)
}  
  de=boot(data=t,statistic=boot.median,R=100)
  ci=boot.ci(de,type = c("norm","basic","perc","bca"))
  ci.norm[i,]=ci$norm[2:3]
  ci.basic[i,]=ci$basic[4:5]
  ci.perc[i,]=ci$percent[4:5]
}
p1.norm=c(mean(ci.norm[,1]<=mu&ci.norm[,2]>=mu),mean(ci.norm[,1]>mu),mean(ci.norm[,2]<mu))
p1.basic=c(mean(ci.basic[,1]<=mu&ci.basic[,2]>=mu),mean(ci.basic[,1]>mu),mean(ci.basic[,2]<mu))
p1.perc=c(mean(ci.perc[,1]<=mu&ci.perc[,2]>=mu),mean(ci.perc[,1]>mu),mean(ci.perc[,2]<mu))
data.frame(p1.norm,p1.basic,p1.perc,row.names = c("coverage","on the left","on the right"))
# χ2(5) distributions
library(moments)
mu1=(8/5)^(1/2)
t=numeric(n)
library(boot)
set.seed(12345)
for(i in 1:m){
for(j in 1:n){
  x=rchisq(n,5)
  t[j]=sk(x)
}  
  de=boot(data=t,statistic=boot.median,R=1000)
  ci=boot.ci(de,type = c("norm","basic","perc","bca"))
  ci.norm[i,]=ci$norm[2:3]
  ci.basic[i,]=ci$basic[4:5]
  ci.perc[i,]=ci$percent[4:5]
}
p2.norm=c(mean(ci.norm[,1]<=mu1&ci.norm[,2]>=mu1),mean(ci.norm[,1]>mu1),mean(ci.norm[,2]<mu1))
p2.basic=c(mean(ci.basic[,1]<=mu1&ci.basic[,2]>=mu1),mean(ci.basic[,1]>mu1),mean(ci.basic[,2]<mu1))
p2.perc=c(mean(ci.perc[,1]<=mu1&ci.perc[,2]>=mu1),mean(ci.perc[,1]>mu1),mean(ci.perc[,2]<mu1))
data.frame(p2.norm,p2.basic,p2.perc,row.names = c("coverage","on the left","on the right"))

## ----warning=FALSE------------------------------------------------------------
# 7.8
n=5;m=88
library(bootstrap)
COV.0=cov(scor)
COV=(87/88)*COV.0
eigen=c(eigen(COV,only.values=TRUE))
e=c(eigen$values[1:1],eigen$values[2:2],eigen$values[3:3],eigen$values[4:4],eigen$values[5:5])
theta.hat=max(e)/sum(e)
theta.jack=numeric(m)
e.hat=numeric(n)
for(i in 1:m){
 COV.hat0=cov(scor[-i,])
 COV.hat=(87/88)*COV.hat0
 eigen.hat=eigen(COV.hat)
 for (j in 1:n) {
   e.hat[j]=eigen.hat$values[j:j]
 }
 theta.jack[i]=max(e.hat)/sum(e.hat)
}
bias.jack=(m-1)*(mean(theta.jack)-theta.hat)
se.jack=sqrt((m-1)*mean((theta.jack-theta.hat)^2))
bias.jack;se.jack
# 7.10
library(DAAG); attach(ironslag)
# estimate
a=seq(10, 40, .1)
## linear
J10=lm(magnetic ~ chemical)
yhat10=J10$coef[1] + J10$coef[2] *a
R10=summary(J10)
## quadratic
J20=lm(magnetic~ chemical + I(chemical^2))
yhat20=J20$coef[1] + J20$coef[2] * a+
J20$coef[3] * a^2
R20=summary(J20)
## exponential
J30=lm(log(magnetic) ~ chemical)
logyhat30=J30$coef[1] + J30$coef[2] * a
R30=summary(J30)
## log-log
J40=lm(log(magnetic) ~ log(chemical))
logyhat40=J40$coef[1] + J40$coef[2] * log(a)
R40=summary(J40)
## cubic
J50=lm(magnetic~chemical+I(chemical^2)+I(chemical^3))
yhat5=J50$coef[1]+J50$coef[2] * a+J50$coef[3] * a^2+J50$coef[4] * a^3
R50=summary(J50)
R1=summary(J10);R2=summary(J20);R3=summary(J30);R4=summary(J40);R5=summary(J50)
R.squard=c(R1$adj.r.squared,R2$adj.r.squared,R3$adj.r.squared,R4$adj.r.squared,R5$adj.r.squared)
# compute the error
n =length(magnetic)
e1=e2=e3=e4=e5=numeric(n)
for (k in 1:n) {
y=magnetic[-k]
x=chemical[-k]
## linear
J1=lm(y ~ x)
yhat1=J1$coef[1] + J1$coef[2] * chemical[k]
e1[k]=magnetic[k] - yhat1

## quadratic
J2=lm(y ~ x + I(x^2))
yhat2=J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k]=magnetic[k] - yhat2

## exponential
J3=lm(log(y) ~ x)
logyhat3=J3$coef[1] + J3$coef[2] * chemical[k]
yhat3=exp(logyhat3)
e3[k]=magnetic[k] - yhat3

## log-log
J4=lm(log(y) ~ log(x))
logyhat4=J4$coef[1] + J4$coef[2] * log(chemical[k])
yhat4=exp(logyhat4)
e4[k]=magnetic[k] - yhat4

## cubic
J5=lm(y~x+I(x^2)+I(x^3))
yhat5=J5$coef[1]+J5$coef[2] * chemical[k]+J5$coef[3] * chemical[k]^2+J5$coef[4] * chemical[k]^3
e5[k]=magnetic[k]-yhat5

}
error=c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2),mean(e5^2))
R.squard=c(R1$adj.r.squared,R2$adj.r.squared,R3$adj.r.squared,R4$adj.r.squared,R5$adj.r.squared)
data.frame(error,R.squard,row.names = c("linear","quadratic","exponential","log-log","cubic"))
J20

## -----------------------------------------------------------------------------
# 8.3
n=1e3
K=1:50
m=20
power=x=y=z=numeric(n)
set.seed(12345)
# bulid a function
count5test=function(x, y) {
X=x-mean(x)
Y=x-mean(y)
outx=sum(X > max(Y)) + sum(X < min(Y))
outy=sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(outx, outy)) > 5))
}
# permutation test method
p=replicate(n,expr={
  x=rnorm(20,0,1)
  y=rnorm(30,0,1)
  x=x-mean(x)
  y=y-mean(y)
  z=c(x,y)
for(i in 1:n){
  k=sample(K,m,replace = FALSE)
  x1=z[k]
  y1=z[-k]
  power[i]=count5test(x1,y1)
}
  mean(power)
  })
p1=mean(p)
p1
# page 31
dCov=function(x, y) {
x=as.matrix(x); y=as.matrix(y)
n=nrow(x); m=nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl=function(x) {
d=as.matrix(dist(x))
m=rowMeans(d); M=mean(d)
a=sweep(d, 1, m); b=sweep(a, 2, m)
b + M
}
A=Akl(x); B=Akl(y)
sqrt(mean(A * B))
}
ndCov2=function(z, ix, dims) {
#dims contains dimensions of x and y
p=dims[1]
q=dims[2]
d=p + q
x=z[ , 1:p] #leave x as is
y=z[ix, -(1:p)] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}
set.seed(123)
library(boot)
library(MASS)
library(Ball)
m=10
sigma=matrix(c(1,0,0,1),2,2)
p.cor=p.ball=numeric(m)
# model 1
p=function(n){
for(i in 1:m){
x=mvrnorm(n,rep(0,2),sigma)
e=mvrnorm(n,rep(0,2),sigma)
y=(1/4)*x+e
z=cbind(x,y)
boot.obj= boot(data = z, statistic = ndCov2, R = 99,sim = "permutation", dims = c(2, 2))
tb = c(boot.obj$t0, boot.obj$t)
p.cor[i] = mean(tb>=tb[1])
p.ball[i] = bcov.test(z[,1:2],z[,3:4],R=99)$p.value
}
c(1-mean(p.cor),1-mean(p.ball))}
t=rbind(p(5),p(10),p(20),p(30),p(40),p(50))
plot(c(5,seq(10,50,10)),t[,1],'l',ylim = c(0,1),xlab='n',main = 'model 1',ylab = 'power',col="blue")
lines(c(5,seq(10,50,10)),t[,2],col='red')
legend(35,0.4,c('p.cor','p.ball'),lty=c(1,1),
       col=c('blue','red'),cex = 0.7)
# model 2
f=function(n){
for(i in 1:m){
x=mvrnorm(n,rep(0,2),sigma)
e=mvrnorm(n,rep(0,2),sigma)
y=(1/4)*x*e
z=cbind(x,y)
boot.obj= boot(data = z, statistic = ndCov2, R = 99,sim = "permutation", dims = c(2, 2))
tb = c(boot.obj$t0, boot.obj$t)
p.cor[i] = mean(tb>=tb[1])
p.ball[i] = bcov.test(z[,1:2],z[,3:4],R=99)$p.value
}
c(1-mean(p.cor),1-mean(p.ball))}
t=rbind(f(5),f(10),f(20),f(30),f(40),f(50))
plot(c(5,seq(10,50,10)),t[,1],'l',ylim = c(0,1),xlab='n',main = 'model 2',ylab = 'power',col="blue")
lines(c(5,seq(10,50,10)),t[,2],col='red')
legend(35,0.4,c('p.cor','p.ball'),lty=c(1,1),
       col=c('blue','red'),cex = 0.7)

## -----------------------------------------------------------------------------
N=1e3
df=function(x)(1/2)*exp(-abs(x))
rw.Metropolis = function(sigma, x0, N) {
x=numeric(N)
x[1] = x0
u=runif(N)
k = 0
for (i in 2:N) {
y = rnorm(1, x[i-1], sigma)
if (u[i] <= (df(y) / df(x[i-1])))
x[i] = y else {
x[i] = x[i-1]
k = k + 1
}
}
return(list(x=x, k=k))
}
# generated for different variances σ2 of the proposal distribution
sigma=c(0.05,0.5,1,10)
m=4;p=numeric(m)
x0=10
rw1 = rw.Metropolis(sigma[1], x0, N)
rw2 = rw.Metropolis(sigma[2], x0, N)
rw3 = rw.Metropolis(sigma[3], x0, N)
rw4 = rw.Metropolis(sigma[4], x0, N)
rw=c(rw1$k, rw2$k, rw3$k, rw4$k)
# rates
for(i in 1:m){
        p[i]=rw[i]/N
}
q=1-p
 rw
 p
 q
 # plots
 index=1:1000
 y1=rw1$x[index]
 plot(index, y1, type="l", main="sigma=0.05", ylab="x",col="red")
 y2=rw2$x[index]
 plot(index, y2, type="l", main="sigma=0.5", ylab="x",col="red")
 y3=rw3$x[index]
 plot(index, y3, type="l", main="sigma=1", ylab="x",col="red")
 y4=rw4$x[index]
 plot(index, y4, type="l", main="sigma=10", ylab="x",col="red")

## -----------------------------------------------------------------------------
# 11.1
x=10
f=function(x)log(exp(x))
g=function(x)exp(log(x))
print(f(x)-g(x))
# compare f and g
f(x)-g(x)==0
# nearly equal
isTRUE(all.equal(f(x),g(x)))
# 11.4
c.k=function(k,a){return(sqrt((k*(a^2))/(k+1-a^2)))}
f1=function(u){(1+(u^2)/(k-1))^(-k/2)}
f2=function(u){(1+(u^2)/k)^(-(k+1)/2)}
f=function(a){
(2*gamma(k/2))/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(f1,0,c.k(k-1,a))$value-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(f2,0,c.k(k,a))$value
}  
# compute roots
library(rootSolve)
t=c(4:25,100)
n=length(t)
root=root2=numeric(n)
for (i in 1:n) {
  k=t[i]
  root[i]=uniroot(f,c(0.05,sqrt(k)/2+1))$root
}
f2=function(a){
  pt(sqrt((k-1)*a^2/(k-a^2)),k-1)-pt(sqrt((k*a^2)/(k+1-a^2)),k)
}
f.4=function(a){
  pt(sqrt((k-1)*a^2/(k-a^2)),k-1)-pt(sqrt((k*a^2)/(k+1-a^2)),k)
}
K1=c(4:25,100,500,1000)
n=length(K1)
root.4=numeric(n)
for (i in 1:n) {
  k=K1[i]
  root.4[i]=uniroot(f.4,c(0.5,sqrt(k)/2+1))$root
}
#the roots of 11.5
root
#the roots of 11.4
root.4
# abo
N=1e3
# max. number of the iteration
n1=28;n2=24;n3=41;n4=70
L=c(.5,.4)
# initial estimates 
tol=.Machine$double.eps^0.5
L.old=L+1
E=numeric(N)
for(j in 1:N){
  E[j]=2*L[1]*n1*log(L[1])/(2-L[1]-2*L[2])+2*L[2]*n2*log(L[2])/(2-L[2]-2*L[1])+2*n3*log(1-L[1]-L[2])+n1*(2-2*L[1]-2*L[2])*log(2*L[1]*(1-L[1]-L[2]))/(2-L[1]-2*L[2])+n2*(2-2*L[1]-2*L[2])*log(2*L[2]*(1-L[1]-L[2]))/(2-L[2]-2*L[1])+n4*log(2*L[1]*L[2])
  model=function(x){
    F1=2*L[1]*n1/((2-L[1]-2*L[2])*x[1])-2*n3/(1-x[1]-x[2])+n1*(2-2*L[1]-2*L[2])*(1-2*x[1]-x[2])/((2-L[1]-2*L[2])*x[1]*(1-x[1]-x[2]))-n2*(2-2*L[1]-2*L[2])/((2-L[2]-2*L[1])*(1-x[1]-x[2]))+n4/x[1]
    F2=2*L[2]*n2/((2-L[2]-2*L[1])*x[2])-2*n3/(1-x[1]-x[2])-n1*(2-2*L[1]-2*L[2])/((2-L[1]-2*L[2])*(1-x[1]-x[2]))+n2*(2-2*L[1]-2*L[2])*(1-2*x[2]-x[1])/((2-L[2]-2*L[1])*x[2]*(1-x[1]-x[2]))+n4/x[2]
    c(F1=F1,F2=F2)
  }
  ss=multiroot(f=model,star=c(.1,.1))
  L=ss$root
  # update p and q
  if (sum(abs(L-L.old)/L.old)<tol) break
  L.old=L
}
L.old
plot(E,type = "l")


## ----warning=FALSE------------------------------------------------------------
# 3
mpg=mtcars$mpg
 disp=mtcars$disp
 wt=mtcars$wt
 formulas=list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt)
 n=4;mod=numeric(n)
for(i in 1:n){mod[i]=lm(formulas[[i]])}
mod
lapply(1:4,function(i)lm(formulas[[i]]))
# 4
bootstraps=lapply(1:10, function(i) {
  rows=sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
lapply(bootstraps,lm,formula=mpg~disp)
n=10;model=numeric(n)
for(i in 1:n){
  data=data.frame(bootstraps[[i]])
  model[i]=lm(data$mpg~data$disp)
}
model
# 5
mpg=mtcars$mpg
 disp=mtcars$disp
 wt=mtcars$wt
 formulas=list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt)
mod=lapply(1:4,function(i)lm(formulas[[i]]))
lapply(1:4,function(i)summary(mod[[i]])$r.squared)
bootstraps=lapply(1:10, function(i) {
  rows=sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
Mod=lapply(bootstraps,lm,formula=mpg~disp)
lapply(1:10,function(i)summary(Mod[[i]])$r.squared)
# 3
trials=replicate(100,t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
sapply(trials,"[[",3)
p=function(n){
  p.value=numeric(n)
  for(i in 1:n){
    result=trials[[i]]
    p.value[i]=result$p.value
  }
  p.value
}
p(100)
# 7
library(parallelsugar)
mpg=mtcars$mpg
 disp=mtcars$disp
 wt=mtcars$wt
 formulas=list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt)
mod=lapply(1:4,function(i)lm(formulas[[i]]))
mclapply(1:4,function(i)summary(mod[[i]])$r.squared,mc.cores = 4)

## ----warning=FALSE------------------------------------------------------------
N=1e3
df=function(x)(1/2)*exp(-abs(x))
rw.Metropolis = function(sigma, x0, N) {
x=numeric(N)
x[1] = x0
u=runif(N)
k = 0
for (i in 2:N) {
y = rnorm(1, x[i-1], sigma)
if (u[i] <= (df(y) / df(x[i-1])))
x[i] = y else {
x[i] = x[i-1]
k = k + 1
}
}
return(list(x=x, k=k))
}
# rcpp
library(Rcpp)
cppFunction('List RwMetropolis(double sigma, double x0, int N) {
NumericVector x(N);
x[0] = x0;
NumericVector u=runif(N);
int k = 0;
for (int i=1;i<N;i++) {
double y = rnorm(1, x[i-1], sigma)[0];
if (u[i] <= exp2(abs(x[i-1])-abs(y))){
x[i] = y; }
else {
x[i] = x[i-1];
k += 1;
}
}
return List::create(
    _["x"] = x,
    _["k"] = k
  );
}')
# qqplot
v1=rw.Metropolis(10,10,1e3)$x
v2=RwMetropolis(10,10,1e3)$x
library(qqplotr)
qqplot(v1,v2)
abline(0,1,col='blue',lwd=2)
# compare time
library(microbenchmark)
ts=microbenchmark(rw.Metropolis(10,10,1e3),RwMetropolis(10,10,1e3))
summary(ts)

