##' comparison of gravity-like models on E&W measles data
##' using fade-out survival analysis

load("ch8.RData")

require(stats4)

par <- list(theta=-15, rho=0, tau1=0, tau2=0)

ySel <- t(MuSel)/PuSel
fadeout <- ySel==0
S <- t(SuSel)
beta <- t(exp(betaSel))

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
  #print(sum(!is.finite(prob)))
  -sum(log(ifelse(obs,prob,1-prob)))
}

system.time(do.call(likfn,par))

test <- mle(likfn,start=par,control=list(trace=4))

exp(test@coef)

summary(test)
vcov(test)

ctest <- confint(test)
ptest <- profile(test)
