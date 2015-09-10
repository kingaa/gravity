
###MLE in 1.9

Pu<-apply(ewPu4464, 2, median)

#library(ncf)

#distance from London
xdist <- rep(0, 954)
for(i in 1:954) {
      xdist[i] <- gcdist(ewxyu4464[442, 1],ewxyu4464[442, 2], ewxyu4464[i, 1], ewxyu4464[i, 2])
      }

sel<-xdist<25
sel[911:912]<-FALSE
MuSel<-round((ewMu2w4464/urMu2w4464.df2.5)[,sel])
ZuSel<-ewZu2w4464.df2.5[,sel]
xySel<-ewxyu4464[sel,]
PuSel<-Pu[sel]

distSel<-matrix(NA, ncol=length(PuSel), nrow=length(PuSel))
for(i in 1:length(PuSel)) {
for(j in i:length(PuSel)) {
      distSel[i,j] <- gcdist(xySel[i, 1],xySel[i, 2], xySel[j, 1], xySel[j, 2])
      distSel[j,i]<-distSel[i,j]
      }
      }

SuSel<-ZuSel
for(i in 1:length(PuSel)){
	SuSel[,i]<-ZuSel[,i]+0.04*PuSel[i]
}
SuSel[SuSel<0]<-1

betaseas<-apply(alphabeta.Shat.df2.5$beta-alphabeta.Shat.df2.5$beta.bar[,1],2,mean,na.rm=T)
seas<-rep(1:26, 21)

betaSel<-matrix(NA, ncol=dim(MuSel)[2], nrow=dim(MuSel)[1])
for(i in 1:length(PuSel)){
	betaSel[,i]<-(3.5976-1.0265*log(PuSel[i]) + betaseas)[seas]
}

alpha<-0.97

rho<- 0
tau1<-1
tau2<-1
theta<- -16

ItSel<-MuSel[-1,]
