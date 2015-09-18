source("c:/o/ucam/measles/oup/data/ewpu4464.q")
#source("c:/o/ucam/measles/oup/data/ewmu4464.q")
#source("c:/o/ucam/measles/oup/data/ewbu4464.q")

#ewMu2w4464<-matrix(NA, nrow=546, ncol=954)
#table(floor(as.numeric(dimnames(ewMu4464)[[1]])-1/52+.01))
# table(floor(as.numeric(dimnames(ewMu4464)[[1]])[-c(209, 470, 783,1096)]-1/52+.01))

#for(i in 1:954){
#  tmp<-split(ewMu4464[-c(209, 470, 783,1096),i], sort(rep(1:546, 2)))
#  ewMu2w4464[,i]<-sapply(tmp, sum)
#  cat(i, "\n")
#}
#dimnames(ewMu2w4464)<-list(44+(0:545)/26, dimnames(ewMu4464)[[2]])
#dump("ewMu2w4464", "c:/o/ucam/measles/oup/data/ewmu2w4464.q")


#ewBu2w4464<-matrix(NA, nrow=546, ncol=954)
#for(i in 1:954){
#  ewBu2w4464[,i]<-(ewBu4464[,i]/26)[sort(rep(1:21,26))]
#  cat(i, "\n")
#}
#dump("ewBu2w4464", "c:/o/ucam/measles/oup/data/ewbu2w4464.q")

source("c:/o/ucam/measles/oup/data/ewbu2w4464.q")
source("c:/o/ucam/measles/oup/data/ewmu2w4464.q")

#urMu2w4464.df2.5<-matrix(NA, nrow=546, ncol=954)
#ewZu2w4464.df2.5<-matrix(NA, nrow=546, ncol=954)

#for(i in (1:954)[-c(911, 912)]){
#  tmp<-smooth.spline(cumsum(ewBu2w4464[,i]), cumsum(ewMu2w4464[,i]), df=2.5)
#  urMu2w4464.df2.5[,i]<-predict(tmp, deriv=1)$y
#  ewZu2w4464.df2.5[,i]<-tmp$y-tmp$yin
#  cat(i, "\n")
#}

#urMu2w4464.df1.1<-matrix(NA, nrow=546, ncol=954)
#ewZu2w4464.df1.1<-matrix(NA, nrow=546, ncol=954)

#for(i in (1:954)[-c(911,912)]){
#  tmp<-smooth.spline(cumsum(ewBu2w4464[,i]), cumsum(ewMu2w4464[,i]), df=1.1)
#  urMu2w4464.df1.1[,i]<-predict(tmp, deriv=1)$y
#  ewZu2w4464.df1.1[,i]<-tmp$y-tmp$yin
#  cat(i, "\n")
#}

#urMu2w4464.df4<-matrix(NA, nrow=546, ncol=954)
#ewZu2w4464.df4<-matrix(NA, nrow=546, ncol=954)

#for(i in (1:954)[-c(911,912)]){
#  tmp<-smooth.spline(cumsum(ewBu2w4464[,i]), cumsum(ewMu2w4464[,i]), df=4)
#  urMu2w4464.df4[,i]<-predict(tmp, deriv=1)$y
#  ewZu2w4464.df4[,i]<-tmp$y-tmp$yin
#  cat(i, "\n")
#}


#dump(c("ewZu2w4464.df1.1", "ewZu2w4464.df2.5","ewZu2w4464.df4"), "c:/o/ucam/measles/oup/data/ewZu2w4464.q")
#dump(c("urMu2w4464.df1.1", "urMu2w4464.df2.5","urMu2w4464.df4"), "c:/o/ucam/measles/oup/data/urMu2w4464.q")

source("c:/o/ucam/measles/oup/data/ewZu2w4464.q")
source("c:/o/ucam/measles/oup/data/urMu2w4464.q")

#TSIR!

Ic.df2.5 <- ewMu2w4464/urMu2w4464.df2.5

seas<-rep(1:26, 21)

measi<-data.frame(It = log(Ic.df2.5[2:546, i]), Ilag = log(Ic.df2.5[1:545, i]), Zlag = ewZu2w4464.df2.5[1:545, i], seas = seas[1:545])

measi[measi==-Inf]<-NA

measi<-na.omit(measi)

tdat2<-split(measi, measi$seas)

N<-dim(measi)[[1]]

tsirll<-function(par, data, n = N, Ps){
  ss<-0
  for(i in 1:26){
  y1bar <- par[i]+log(data[[i]]$Zlag+par[27]*Ps)+par[28]*data[[i]]$Ilag+par[29]/data[[i]]$Ilag
  ss<-ss+sum((data[[i]]$It-y1bar)^2)
 }

 return(-(1/n)*log(sum(ss)))
}

tsirll(par, data=tdat2, n=N, Ps= median(ewPu4464[,i]))


par<-c(rep(-9,26), 0.03, 1, 0)

tmp<-optim(par, tsirll, data=tdat2, n = N, Ps= median(ewPu4464[,i]), method="L-BFGS-B", lower=c(rep(-Inf, 26), 0.02, -Inf, -Inf), upper=c(rep(Inf, 26), 1, Inf, Inf), control=list(trace=2), hessian=T)




###




Ic.df2.5 <- ewMu2w4464/urMu2w4464.df2.5
seas<-rep(1:26, 21)

Ps954<-apply(ewPu4464, 2, median)

Shat.df2.5<-list(sgrid=matrix(NA, ncol=100, nrow=954), dev = matrix(NA, ncol=100, nrow=954))

for(i in c(1:909,913:954)){

  measi<-data.frame(It = Ic.df2.5[2:546, i], Ilag = Ic.df2.5[1:545, i], Zlag = ewZu2w4464.df2.5[1:545, i], seas = seas[1:545])

  tgrid<-seq(Ps954[i]*0.001-min(measi$Zlag), Ps954[i]*0.25, length=100)

  Shat.df2.5$sgrid[i,]<-tgrid

  for(j in 1:100){
    measi$slag<- measi$Zlag + tgrid[j]

    tmpres<-glm(log(It)~-1+as.factor(seas)+log(Ilag)+I(1/Ilag)+offset(log(slag)), data=measi, subset=Ilag>0&It>0)

    Shat.df2.5$dev[i,j]<-tmpres$deviance
  }

cat(i, "\n")

}


Sbar.df2.5<-list(estimate=rep(NA, 954), upper=rep(NA, 954), lower=rep(NA, 954))

for(i in c(1:909,913:954)){

  tmp<-smooth.spline(Shat.df2.5$sgrid[i,],Shat.df2.5$dev[i,])
  finex<-seq(min(Shat.df2.5$sgrid[i,]), max(Shat.df2.5$sgrid[i,]), length=1000)
  tmp2<-predict(tmp, x=finex)$y

  tmp3<-range(finex[tmp2-min(tmp2)<qchisq(0.95,1)/2])

  Sbar.df2.5$upper[i]<-tmp3[2]
  Sbar.df2.5$lower[i]<-tmp3[1]
  Sbar.df2.5$estimate[i]<-finex[order(tmp2)[1]]

  cat(i, "\n")
}


dump(c("Shat.df2.5", "Sbar.df2.5"), "c:/o/ucam/measles/oup/data/ewSudf2.5.q")


summary(lm(Sbar.df2.5$estimate/Ps954~1))
Call: lm(formula = Sbar.df2.5$estimate/Ps954 ~ 1)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.036423 -0.014410 -0.007286  0.002482  0.209166 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 0.040834   0.001076   37.95   <2e-16 ***
---
Signif. codes:  0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1 
Residual standard error: 0.03318 on 950 degrees of freedom


require(MASS)

summary(rlm(Sbar.df2.5$estimate/Ps954~1))

Call: rlm.formula(formula = Sbar.df2.5$estimate/Ps954 ~ 1)
Residuals:
      Min        1Q    Median        3Q       Max 
-0.030459 -0.008446 -0.001322  0.008446  0.215130 

Coefficients:
            Value   Std. Error t value
(Intercept)  0.0349  0.0004    81.1299
Residual standard error: 0.01252 on 950 degrees of freedom


 
alphabeta.Sbar.df2.5<-list(alpha=matrix(NA, ncol=1, nrow=954), beta = matrix(NA, ncol=26, nrow=954), 
   alpha.se=matrix(NA, ncol=1, nrow=954), beta.se = matrix(NA, ncol=26, nrow=954),
   beta.bar=matrix(NA, ncol=1, nrow=954), beta.sd=matrix(NA, ncol=1, nrow=954))

for(i in c(1:909,913:954)){

  measi<-data.frame(It = Ic.df2.5[2:546, i], Ilag = Ic.df2.5[1:545, i], Zlag = ewZu2w4464.df2.5[1:545, i], seas = seas[1:545])
  measi$slag<- measi$Zlag + Ps954[i]*0.040834

  tmpres<-glm(log(It)~-1+as.factor(seas)+log(Ilag)+I(1/Ilag)+offset(log(slag)), data=measi, subset=Ilag>0&It>0)
	
  regmat<-summary(tmpres)$coef
  nregmat<-dimnames(regmat)[[1]]

  alphabeta.Sbar.df2.5$beta[i,1]<-ifelse(length(regmat[nregmat=="as.factor(seas)1",1])==1, regmat[nregmat=="as.factor(seas)1",1], NA)
  alphabeta.Sbar.df2.5$beta[i,2]<-ifelse(length(regmat[nregmat=="as.factor(seas)2",1])==1, regmat[nregmat=="as.factor(seas)2",1], NA)
  alphabeta.Sbar.df2.5$beta[i,3]<-ifelse(length(regmat[nregmat=="as.factor(seas)3",1])==1, regmat[nregmat=="as.factor(seas)3",1], NA)
  alphabeta.Sbar.df2.5$beta[i,4]<-ifelse(length(regmat[nregmat=="as.factor(seas)4",1])==1, regmat[nregmat=="as.factor(seas)4",1], NA)
  alphabeta.Sbar.df2.5$beta[i,5]<-ifelse(length(regmat[nregmat=="as.factor(seas)5",1])==1, regmat[nregmat=="as.factor(seas)5",1], NA)
  alphabeta.Sbar.df2.5$beta[i,6]<-ifelse(length(regmat[nregmat=="as.factor(seas)6",1])==1, regmat[nregmat=="as.factor(seas)6",1], NA)
  alphabeta.Sbar.df2.5$beta[i,7]<-ifelse(length(regmat[nregmat=="as.factor(seas)7",1])==1, regmat[nregmat=="as.factor(seas)7",1], NA)
  alphabeta.Sbar.df2.5$beta[i,8]<-ifelse(length(regmat[nregmat=="as.factor(seas)8",1])==1, regmat[nregmat=="as.factor(seas)8",1], NA)
  alphabeta.Sbar.df2.5$beta[i,9]<-ifelse(length(regmat[nregmat=="as.factor(seas)9",1])==1, regmat[nregmat=="as.factor(seas)9",1], NA)
  alphabeta.Sbar.df2.5$beta[i,10]<-ifelse(length(regmat[nregmat=="as.factor(seas)10",1])==1, regmat[nregmat=="as.factor(seas)10",1], NA)
  alphabeta.Sbar.df2.5$beta[i,11]<-ifelse(length(regmat[nregmat=="as.factor(seas)11",1])==1, regmat[nregmat=="as.factor(seas)11",1], NA)
  alphabeta.Sbar.df2.5$beta[i,12]<-ifelse(length(regmat[nregmat=="as.factor(seas)12",1])==1, regmat[nregmat=="as.factor(seas)12",1], NA)
  alphabeta.Sbar.df2.5$beta[i,13]<-ifelse(length(regmat[nregmat=="as.factor(seas)13",1])==1, regmat[nregmat=="as.factor(seas)13",1], NA)
  alphabeta.Sbar.df2.5$beta[i,14]<-ifelse(length(regmat[nregmat=="as.factor(seas)14",1])==1, regmat[nregmat=="as.factor(seas)14",1], NA)
  alphabeta.Sbar.df2.5$beta[i,15]<-ifelse(length(regmat[nregmat=="as.factor(seas)15",1])==1, regmat[nregmat=="as.factor(seas)15",1], NA)
  alphabeta.Sbar.df2.5$beta[i,16]<-ifelse(length(regmat[nregmat=="as.factor(seas)16",1])==1, regmat[nregmat=="as.factor(seas)16",1], NA)
  alphabeta.Sbar.df2.5$beta[i,17]<-ifelse(length(regmat[nregmat=="as.factor(seas)17",1])==1, regmat[nregmat=="as.factor(seas)17",1], NA)
  alphabeta.Sbar.df2.5$beta[i,18]<-ifelse(length(regmat[nregmat=="as.factor(seas)18",1])==1, regmat[nregmat=="as.factor(seas)18",1], NA)
  alphabeta.Sbar.df2.5$beta[i,19]<-ifelse(length(regmat[nregmat=="as.factor(seas)19",1])==1, regmat[nregmat=="as.factor(seas)19",1], NA)
  alphabeta.Sbar.df2.5$beta[i,20]<-ifelse(length(regmat[nregmat=="as.factor(seas)20",1])==1, regmat[nregmat=="as.factor(seas)20",1], NA)
  alphabeta.Sbar.df2.5$beta[i,21]<-ifelse(length(regmat[nregmat=="as.factor(seas)21",1])==1, regmat[nregmat=="as.factor(seas)21",1], NA)
  alphabeta.Sbar.df2.5$beta[i,22]<-ifelse(length(regmat[nregmat=="as.factor(seas)22",1])==1, regmat[nregmat=="as.factor(seas)22",1], NA)
  alphabeta.Sbar.df2.5$beta[i,23]<-ifelse(length(regmat[nregmat=="as.factor(seas)23",1])==1, regmat[nregmat=="as.factor(seas)23",1], NA)
  alphabeta.Sbar.df2.5$beta[i,24]<-ifelse(length(regmat[nregmat=="as.factor(seas)24",1])==1, regmat[nregmat=="as.factor(seas)24",1], NA)
  alphabeta.Sbar.df2.5$beta[i,25]<-ifelse(length(regmat[nregmat=="as.factor(seas)25",1])==1, regmat[nregmat=="as.factor(seas)25",1], NA)
  alphabeta.Sbar.df2.5$beta[i,26]<-ifelse(length(regmat[nregmat=="as.factor(seas)26",1])==1, regmat[nregmat=="as.factor(seas)26",1], NA)

  alphabeta.Sbar.df2.5$alpha[i,]<-ifelse(length(regmat[nregmat=="log(Ilag)",1])==1, regmat[nregmat=="log(Ilag)",1], NA)

  alphabeta.Sbar.df2.5$beta.se[i,1]<-ifelse(length(regmat[nregmat=="as.factor(seas)1",2])==1, regmat[nregmat=="as.factor(seas)1",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,2]<-ifelse(length(regmat[nregmat=="as.factor(seas)2",2])==1, regmat[nregmat=="as.factor(seas)2",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,3]<-ifelse(length(regmat[nregmat=="as.factor(seas)3",2])==1, regmat[nregmat=="as.factor(seas)3",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,4]<-ifelse(length(regmat[nregmat=="as.factor(seas)4",2])==1, regmat[nregmat=="as.factor(seas)4",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,5]<-ifelse(length(regmat[nregmat=="as.factor(seas)5",2])==1, regmat[nregmat=="as.factor(seas)5",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,6]<-ifelse(length(regmat[nregmat=="as.factor(seas)6",2])==1, regmat[nregmat=="as.factor(seas)6",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,7]<-ifelse(length(regmat[nregmat=="as.factor(seas)7",2])==1, regmat[nregmat=="as.factor(seas)7",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,8]<-ifelse(length(regmat[nregmat=="as.factor(seas)8",2])==1, regmat[nregmat=="as.factor(seas)8",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,9]<-ifelse(length(regmat[nregmat=="as.factor(seas)9",2])==1, regmat[nregmat=="as.factor(seas)9",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,10]<-ifelse(length(regmat[nregmat=="as.factor(seas)10",2])==1, regmat[nregmat=="as.factor(seas)10",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,11]<-ifelse(length(regmat[nregmat=="as.factor(seas)11",2])==1, regmat[nregmat=="as.factor(seas)11",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,12]<-ifelse(length(regmat[nregmat=="as.factor(seas)12",2])==1, regmat[nregmat=="as.factor(seas)12",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,13]<-ifelse(length(regmat[nregmat=="as.factor(seas)13",2])==1, regmat[nregmat=="as.factor(seas)13",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,14]<-ifelse(length(regmat[nregmat=="as.factor(seas)14",2])==1, regmat[nregmat=="as.factor(seas)14",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,15]<-ifelse(length(regmat[nregmat=="as.factor(seas)15",2])==1, regmat[nregmat=="as.factor(seas)15",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,16]<-ifelse(length(regmat[nregmat=="as.factor(seas)16",2])==1, regmat[nregmat=="as.factor(seas)16",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,17]<-ifelse(length(regmat[nregmat=="as.factor(seas)17",2])==1, regmat[nregmat=="as.factor(seas)17",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,18]<-ifelse(length(regmat[nregmat=="as.factor(seas)18",2])==1, regmat[nregmat=="as.factor(seas)18",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,19]<-ifelse(length(regmat[nregmat=="as.factor(seas)19",2])==1, regmat[nregmat=="as.factor(seas)19",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,20]<-ifelse(length(regmat[nregmat=="as.factor(seas)20",2])==1, regmat[nregmat=="as.factor(seas)20",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,21]<-ifelse(length(regmat[nregmat=="as.factor(seas)21",2])==1, regmat[nregmat=="as.factor(seas)21",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,22]<-ifelse(length(regmat[nregmat=="as.factor(seas)22",2])==1, regmat[nregmat=="as.factor(seas)22",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,23]<-ifelse(length(regmat[nregmat=="as.factor(seas)23",2])==1, regmat[nregmat=="as.factor(seas)23",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,24]<-ifelse(length(regmat[nregmat=="as.factor(seas)24",2])==1, regmat[nregmat=="as.factor(seas)24",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,25]<-ifelse(length(regmat[nregmat=="as.factor(seas)25",2])==1, regmat[nregmat=="as.factor(seas)25",2], NA)
  alphabeta.Sbar.df2.5$beta.se[i,26]<-ifelse(length(regmat[nregmat=="as.factor(seas)26",2])==1, regmat[nregmat=="as.factor(seas)26",2], NA)

  alphabeta.Sbar.df2.5$alpha.se[i,]<-ifelse(length(regmat[nregmat=="log(Ilag)",2])==1, regmat[nregmat=="log(Ilag)",2], NA)

  alphabeta.Sbar.df2.5$beta.bar[i,]<-mean(alphabeta.Sbar.df2.5$beta[i,], na.rm=TRUE)
  alphabeta.Sbar.df2.5$beta.sd[i,]<-sd(alphabeta.Sbar.df2.5$beta[i,], na.rm=TRUE)

  cat(i, "\n")

}


tmp<-rlm(alphabeta.Sbar.df2.5$beta.bar~log(Ps954), weights= 1/(alphabeta.Sbar.df2.5$beta.sd^2)[,1])
summary(tmp)
Call: rlm.formula(formula = alphabeta.Sbar.df2.5$beta.bar ~ log(Ps954), weights = 1/(alphabeta.Sbar.df2.5$beta.sd^2)[, 1])
Residuals:
       Min         1Q     Median         3Q        Max 
-291.50424   -0.22893    0.04603    0.33694  245.03587 

Coefficients:
            Value    Std. Error t value 
(Intercept)   3.5628   0.1264    28.1798
log(Ps954)   -1.0422   0.0130   -80.3147

Residual standard error: 0.3931 on 949 degrees of freedom

######

alphabeta.Shat.df2.5<-list(alpha=matrix(NA, ncol=1, nrow=954), beta = matrix(NA, ncol=26, nrow=954), 
   alpha.se=matrix(NA, ncol=1, nrow=954), beta.se = matrix(NA, ncol=26, nrow=954),
   beta.bar=matrix(NA, ncol=1, nrow=954), beta.sd=matrix(NA, ncol=1, nrow=954))

for(i in c(1:909,913:954)){

  measi<-data.frame(It = Ic.df2.5[2:546, i], Ilag = Ic.df2.5[1:545, i], Zlag = ewZu2w4464.df2.5[1:545, i], seas = seas[1:545])
  measi$slag<- measi$Zlag + Sbar.df2.5$estimate[i]

  tmpres<-glm(log(It)~-1+as.factor(seas)+log(Ilag)+I(1/Ilag)+offset(log(slag)), data=measi, subset=Ilag>0&It>0)
	
  regmat<-summary(tmpres)$coef
  nregmat<-dimnames(regmat)[[1]]

  alphabeta.Shat.df2.5$beta[i,1]<-ifelse(length(regmat[nregmat=="as.factor(seas)1",1])==1, regmat[nregmat=="as.factor(seas)1",1], NA)
  alphabeta.Shat.df2.5$beta[i,2]<-ifelse(length(regmat[nregmat=="as.factor(seas)2",1])==1, regmat[nregmat=="as.factor(seas)2",1], NA)
  alphabeta.Shat.df2.5$beta[i,3]<-ifelse(length(regmat[nregmat=="as.factor(seas)3",1])==1, regmat[nregmat=="as.factor(seas)3",1], NA)
  alphabeta.Shat.df2.5$beta[i,4]<-ifelse(length(regmat[nregmat=="as.factor(seas)4",1])==1, regmat[nregmat=="as.factor(seas)4",1], NA)
  alphabeta.Shat.df2.5$beta[i,5]<-ifelse(length(regmat[nregmat=="as.factor(seas)5",1])==1, regmat[nregmat=="as.factor(seas)5",1], NA)
  alphabeta.Shat.df2.5$beta[i,6]<-ifelse(length(regmat[nregmat=="as.factor(seas)6",1])==1, regmat[nregmat=="as.factor(seas)6",1], NA)
  alphabeta.Shat.df2.5$beta[i,7]<-ifelse(length(regmat[nregmat=="as.factor(seas)7",1])==1, regmat[nregmat=="as.factor(seas)7",1], NA)
  alphabeta.Shat.df2.5$beta[i,8]<-ifelse(length(regmat[nregmat=="as.factor(seas)8",1])==1, regmat[nregmat=="as.factor(seas)8",1], NA)
  alphabeta.Shat.df2.5$beta[i,9]<-ifelse(length(regmat[nregmat=="as.factor(seas)9",1])==1, regmat[nregmat=="as.factor(seas)9",1], NA)
  alphabeta.Shat.df2.5$beta[i,10]<-ifelse(length(regmat[nregmat=="as.factor(seas)10",1])==1, regmat[nregmat=="as.factor(seas)10",1], NA)
  alphabeta.Shat.df2.5$beta[i,11]<-ifelse(length(regmat[nregmat=="as.factor(seas)11",1])==1, regmat[nregmat=="as.factor(seas)11",1], NA)
  alphabeta.Shat.df2.5$beta[i,12]<-ifelse(length(regmat[nregmat=="as.factor(seas)12",1])==1, regmat[nregmat=="as.factor(seas)12",1], NA)
  alphabeta.Shat.df2.5$beta[i,13]<-ifelse(length(regmat[nregmat=="as.factor(seas)13",1])==1, regmat[nregmat=="as.factor(seas)13",1], NA)
  alphabeta.Shat.df2.5$beta[i,14]<-ifelse(length(regmat[nregmat=="as.factor(seas)14",1])==1, regmat[nregmat=="as.factor(seas)14",1], NA)
  alphabeta.Shat.df2.5$beta[i,15]<-ifelse(length(regmat[nregmat=="as.factor(seas)15",1])==1, regmat[nregmat=="as.factor(seas)15",1], NA)
  alphabeta.Shat.df2.5$beta[i,16]<-ifelse(length(regmat[nregmat=="as.factor(seas)16",1])==1, regmat[nregmat=="as.factor(seas)16",1], NA)
  alphabeta.Shat.df2.5$beta[i,17]<-ifelse(length(regmat[nregmat=="as.factor(seas)17",1])==1, regmat[nregmat=="as.factor(seas)17",1], NA)
  alphabeta.Shat.df2.5$beta[i,18]<-ifelse(length(regmat[nregmat=="as.factor(seas)18",1])==1, regmat[nregmat=="as.factor(seas)18",1], NA)
  alphabeta.Shat.df2.5$beta[i,19]<-ifelse(length(regmat[nregmat=="as.factor(seas)19",1])==1, regmat[nregmat=="as.factor(seas)19",1], NA)
  alphabeta.Shat.df2.5$beta[i,20]<-ifelse(length(regmat[nregmat=="as.factor(seas)20",1])==1, regmat[nregmat=="as.factor(seas)20",1], NA)
  alphabeta.Shat.df2.5$beta[i,21]<-ifelse(length(regmat[nregmat=="as.factor(seas)21",1])==1, regmat[nregmat=="as.factor(seas)21",1], NA)
  alphabeta.Shat.df2.5$beta[i,22]<-ifelse(length(regmat[nregmat=="as.factor(seas)22",1])==1, regmat[nregmat=="as.factor(seas)22",1], NA)
  alphabeta.Shat.df2.5$beta[i,23]<-ifelse(length(regmat[nregmat=="as.factor(seas)23",1])==1, regmat[nregmat=="as.factor(seas)23",1], NA)
  alphabeta.Shat.df2.5$beta[i,24]<-ifelse(length(regmat[nregmat=="as.factor(seas)24",1])==1, regmat[nregmat=="as.factor(seas)24",1], NA)
  alphabeta.Shat.df2.5$beta[i,25]<-ifelse(length(regmat[nregmat=="as.factor(seas)25",1])==1, regmat[nregmat=="as.factor(seas)25",1], NA)
  alphabeta.Shat.df2.5$beta[i,26]<-ifelse(length(regmat[nregmat=="as.factor(seas)26",1])==1, regmat[nregmat=="as.factor(seas)26",1], NA)

  alphabeta.Shat.df2.5$alpha[i,]<-ifelse(length(regmat[nregmat=="log(Ilag)",1])==1, regmat[nregmat=="log(Ilag)",1], NA)

  alphabeta.Shat.df2.5$beta.se[i,1]<-ifelse(length(regmat[nregmat=="as.factor(seas)1",2])==1, regmat[nregmat=="as.factor(seas)1",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,2]<-ifelse(length(regmat[nregmat=="as.factor(seas)2",2])==1, regmat[nregmat=="as.factor(seas)2",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,3]<-ifelse(length(regmat[nregmat=="as.factor(seas)3",2])==1, regmat[nregmat=="as.factor(seas)3",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,4]<-ifelse(length(regmat[nregmat=="as.factor(seas)4",2])==1, regmat[nregmat=="as.factor(seas)4",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,5]<-ifelse(length(regmat[nregmat=="as.factor(seas)5",2])==1, regmat[nregmat=="as.factor(seas)5",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,6]<-ifelse(length(regmat[nregmat=="as.factor(seas)6",2])==1, regmat[nregmat=="as.factor(seas)6",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,7]<-ifelse(length(regmat[nregmat=="as.factor(seas)7",2])==1, regmat[nregmat=="as.factor(seas)7",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,8]<-ifelse(length(regmat[nregmat=="as.factor(seas)8",2])==1, regmat[nregmat=="as.factor(seas)8",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,9]<-ifelse(length(regmat[nregmat=="as.factor(seas)9",2])==1, regmat[nregmat=="as.factor(seas)9",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,10]<-ifelse(length(regmat[nregmat=="as.factor(seas)10",2])==1, regmat[nregmat=="as.factor(seas)10",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,11]<-ifelse(length(regmat[nregmat=="as.factor(seas)11",2])==1, regmat[nregmat=="as.factor(seas)11",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,12]<-ifelse(length(regmat[nregmat=="as.factor(seas)12",2])==1, regmat[nregmat=="as.factor(seas)12",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,13]<-ifelse(length(regmat[nregmat=="as.factor(seas)13",2])==1, regmat[nregmat=="as.factor(seas)13",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,14]<-ifelse(length(regmat[nregmat=="as.factor(seas)14",2])==1, regmat[nregmat=="as.factor(seas)14",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,15]<-ifelse(length(regmat[nregmat=="as.factor(seas)15",2])==1, regmat[nregmat=="as.factor(seas)15",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,16]<-ifelse(length(regmat[nregmat=="as.factor(seas)16",2])==1, regmat[nregmat=="as.factor(seas)16",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,17]<-ifelse(length(regmat[nregmat=="as.factor(seas)17",2])==1, regmat[nregmat=="as.factor(seas)17",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,18]<-ifelse(length(regmat[nregmat=="as.factor(seas)18",2])==1, regmat[nregmat=="as.factor(seas)18",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,19]<-ifelse(length(regmat[nregmat=="as.factor(seas)19",2])==1, regmat[nregmat=="as.factor(seas)19",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,20]<-ifelse(length(regmat[nregmat=="as.factor(seas)20",2])==1, regmat[nregmat=="as.factor(seas)20",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,21]<-ifelse(length(regmat[nregmat=="as.factor(seas)21",2])==1, regmat[nregmat=="as.factor(seas)21",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,22]<-ifelse(length(regmat[nregmat=="as.factor(seas)22",2])==1, regmat[nregmat=="as.factor(seas)22",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,23]<-ifelse(length(regmat[nregmat=="as.factor(seas)23",2])==1, regmat[nregmat=="as.factor(seas)23",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,24]<-ifelse(length(regmat[nregmat=="as.factor(seas)24",2])==1, regmat[nregmat=="as.factor(seas)24",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,25]<-ifelse(length(regmat[nregmat=="as.factor(seas)25",2])==1, regmat[nregmat=="as.factor(seas)25",2], NA)
  alphabeta.Shat.df2.5$beta.se[i,26]<-ifelse(length(regmat[nregmat=="as.factor(seas)26",2])==1, regmat[nregmat=="as.factor(seas)26",2], NA)

  alphabeta.Shat.df2.5$alpha.se[i,]<-ifelse(length(regmat[nregmat=="log(Ilag)",2])==1, regmat[nregmat=="log(Ilag)",2], NA)

  alphabeta.Shat.df2.5$beta.bar[i,]<-mean(alphabeta.Shat.df2.5$beta[i,], na.rm=TRUE)
  alphabeta.Shat.df2.5$beta.sd[i,]<-sd(alphabeta.Shat.df2.5$beta[i,], na.rm=TRUE)

  cat(i, "\n")

}


tmp<-rlm(alphabeta.Shat.df2.5$beta.bar~log(Ps954), weights= 1/(alphabeta.Shat.df2.5$beta.sd^2)[,1])
summary(tmp)

Call: rlm.formula(formula = alphabeta.Shat.df2.5$beta.bar ~ log(Ps954), 
    weights = 1/(alphabeta.Shat.df2.5$beta.sd^2)[, 1])
Residuals:
       Min         1Q     Median         3Q        Max 
 -7.443288  -0.353264  -0.000848   0.441144 211.625818 

Coefficients:
            Value    Std. Error t value 
(Intercept)   3.5976   0.1820    19.7636
log(Ps954)   -1.0265   0.0187   -54.9391

plot(alphabeta.Shat.df2.5$beta.bar~log(Ps954), ylim=c(-15,10))
abline(tmp)




plot(apply(alphabeta.Shat.df2.5$beta-alphabeta.Shat.df2.5$beta.bar[,1],2,median,na.rm=T))
points(apply(alphabeta.Shat.df2.5$beta-alphabeta.Shat.df2.5$beta.bar[,1],2,mean,na.rm=T), col='red')

dump(c("alphabeta.Shat.df2.5", "alphabeta.Sbar.df2.5"), "c:/o/ucam/measles/oup/data/alphabetadf2.5.q")

