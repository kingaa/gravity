
###########################################
load("/Users/onb1/Documents/UCam/Measles/oup/R/ch8_gravity/ch8.RData")

require(stats4)

par<-c(theta=-11.57, rho=0.045, tau1=0, tau2=0)

ySel=MuSel
for(i in 1:952) ySel[,i]<-MuSel[,i]/PuSel[i]

likfn<-function(theta=-11.57, rho=0.045, tau1=0, tau2=0){
iotaSel<-matrix(NA, ncol=dim(MuSel)[2], nrow=dim(MuSel)[1])
for(i in 1:length(PuSel)){
    tmp<-t(t(ySel*PuSel^exp(tau2))/(distSel[,i]^exp(rho)))
    iotaSel[,i]<-(PuSel[i]^exp(tau1))*apply(tmp[,-i],1,sum)
}

lambdaSel<-exp(betaSel[-546,])*SuSel[-546,]*(exp(theta)*iotaSel[-546,])^alpha

mu<-lambdaSel[MuSel[-546,]==0]
obs<-ItSel[MuSel[-546,]==0]


#-sum(dnbinom(x=obs, size=exp(par[1])*iotaSel[-546,], mu=mu, log=T))
-sum(dpois(obs, mu, log=T))
}

date();likfn(par);date()
test<-mle(likfn, start=list(theta=-11.57, rho=0.045, tau1=0, tau2=0))


par2<-c(theta=-11.57, rho=0.045, tau1=0, tau2=0, del=0)

likfncd<-function(par2){
iotaSel<-matrix(NA, ncol=dim(MuSel)[2], nrow=dim(MuSel)[1])
for(i in 1:length(PuSel)){
	for(j in 1:length(PuSel)){
	PuSel[k]^exp(par[4])/(distSel[,k])	
	}	
    tmp<-t(t(ySel*PuSel^exp(par[4]))/(distSel[,i]^exp(par[2])))
    iotaSel[,i]<-(PuSel[i]^exp(par[3]))*apply(tmp[,-i],1,sum)
}

lambdaSel<-exp(betaSel[-546,])*SuSel[-546,]*(exp(par[1])*iotaSel[-546,])^alpha

mu<-lambdaSel[MuSel[-546,]==0]
obs<-ItSel[MuSel[-546,]==0]


#-sum(dnbinom(x=obs, size=exp(par[1])*iotaSel[-546,], mu=mu, log=T))
-sum(dpois(obs, mu, log=T))
}


