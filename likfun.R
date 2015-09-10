##' comparison of gravity-like models on E&W measles data
##' using fade-out survival analysis

load("ch8.RData")

require(stats4)

par <- c(theta=-11.57, rho=0.045, tau1=0, tau2=0)

ySel <- t(MuSel)/PuSel
fadeout <- ySel==0
S <- t(SuSel)
beta <- t(betaSel)

relevant <- which(fadeout[,-546])
obs <- fadeout[,-1][relevant]

require(doMC)
require(foreach)
registerDoMC(4)

likfn <- function (theta, rho, tau1, tau2) {
  theta <- exp(theta)
  rho <- exp(rho)
  tau1 <- exp(tau1)
  tau2 <- exp(tau2)
  pd <- PuSel^tau2*distSel^(-rho)
  diag(pd) <- 0
  iotaSel <- theta*(PuSel^tau1)*(pd%*%ySel)
  lambda <- beta*S*iotaSel^alpha
  prob <- exp(-lambda[relevant])
  print(sum(!is.finite(prob)))
  -sum(log(ifelse(obs,prob,1-prob)))
}

with(as.list(par),likfn(theta,rho,tau1,tau2))

  lambdaSel<-exp(betaSel[-546,])*SuSel[-546,]*(exp(theta)*iotaSel[-546,])^alpha

  mu<-lambdaSel[MuSel[-546,]==0]
  obs<-ItSel[MuSel[-546,]==0]


  #-sum(dnbinom(x=obs, size=exp(par[1])*iotaSel[-546,], mu=mu, log=T))
  -sum(dpois(obs, mu, log=T))
}

date();likfn(par);date()
test<-mle(likfn, start=list(theta=-11.57, rho=0.045, tau1=0, tau2=0))

