options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8")

set.seed(594709947L)

library(plyr)
library(reshape2)
library(magrittr)
library(foreach)
library(iterators)
library(bbmle)
library(nloptr)
library(pomp)

## ----gravity-pre---------------------------------------------------------
readRDS("tsir-fits.rds") -> dat
dat %>% acast(town~year+biweek,value.var="ylag") -> ylag
dat %>% acast(town~year+biweek,value.var="Slag") -> slag
dat %>% acast(town~year+biweek,value.var="log.beta") %>% exp() -> beta
dat %>% acast(town~year+biweek,value.var="cases") -> obs
dat %>% daply(~town,function(x)unique(x$Ps)) -> N

readRDS("distances.rds") -> distances
dd <- 1/distances
diag(dd) <- 0

stopifnot(identical(rownames(ylag),rownames(slag)))
stopifnot(identical(rownames(ylag),rownames(beta)))
stopifnot(identical(rownames(ylag),rownames(obs)))
stopifnot(identical(rownames(ylag),names(N)))
stopifnot(identical(rownames(ylag),rownames(dd)))
stopifnot(identical(rownames(ylag),colnames(dd)))

ddscal <- exp(mean(log(dd[lower.tri(dd)])))
dd <- dd/ddscal

alpha <- mean(dat$alpha)

relevant <- which(ylag==0&slag>=0)
obs <- obs[relevant]
betaS <- beta[relevant]*slag[relevant]

## ----gravity-like--------------------------------------------------------
likfnGravity <- function (theta, phi, rho, psi, tau1, tau2) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  Q <- (N^tau2)*(dd^rho)
  iota <- (N^tau1)*(theta*crossprod(Q,ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}

## ----gravity-profile-2D,cache=FALSE--------------------------------------
bake(file="gravity.rds",{

  grd <- profileDesign(
    tau1=seq(0.1, 1.4, length=25),
    tau2=seq(-1, 1, length=25),
    lower=c(psi=log(5),theta=log(1e-8),phi=log(1e-6),rho=log(0.9)),
    upper=c(psi=log(7),theta=log(1),phi=log(1e-2),rho=log(1.7)),
    nprof=10
  )

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau1","tau2")
      formals(likfnGravity) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
             function(x){
               do.call(likfnGravity,
                       c(start[fixed],setNames(as.list(x),est)))
             },
             opts=list(
               algorithm="NLOPT_LN_SBPLX",
               ftol_abs=1000,
               maxeval=10000)
             ) -> fit
      fit$solution %>% setNames(est) %>%
        c(unlist(start[fixed]),loglik=-fit$objective) %>%
        as.list() %>% as.data.frame() %>%
        cbind(conv=fit$status) -> res
     }
  ) -> results

  results %>%
    extract(sapply(results,inherits,"try-error")) %>%
    sapply(as.character) %>%
    unique() %>%
    print()

  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    ddply(~tau1+tau2,subset,loglik==max(loglik,na.rm=TRUE)) -> grd

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau1","tau2")
      formals(likfnGravity) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
             function(x){
               do.call(likfnGravity,
                       c(start[fixed],setNames(as.list(x),est)))
             },
             opts=list(
               algorithm="NLOPT_LN_SBPLX",
               ftol_abs=1e-5,
               xtol_rel=1e-7,
               maxeval=10000)
      ) -> fit
      fit$solution %>% setNames(est) %>%
        c(unlist(start[fixed]),loglik=-fit$objective) %>%
        as.list() %>% as.data.frame() %>%
        cbind(conv=fit$status) -> res
    }
  ) -> results


  attr(results,"nproc") <- getDoParWorkers()
  results
}) -> results

## ----xia-like------------------------------------------------------------
likfnXia <- function (theta, phi, rho, psi, tau1, tau2) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  Q <- (N^tau2)*(dd^rho)
  iota <- (N^tau1)*(theta*crossprod(Q,ylag^tau2)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}

## ----xia-profile-2D,cache=FALSE------------------------------------------
bake(file="xia.rds",{

  grd <- expand.grid(
    tau1=seq(0.1, 1.5, length=25),
    tau2=seq(0.1, 1.0, length=25),
    psi=log(5),
    rho=log(1.0),
    theta=log(1e-6),
    phi=log(1e-6)
  )

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages="bbmle"
  ) %dopar% try(
    {
      mle2(likfnXia,
           method="Nelder-Mead",
           start=as.list(start[c("rho","theta","psi","phi")]),
           fixed=as.list(start[c("tau1","tau2")]),
           control=list(trace=0,maxit=10000),
           skip.hessian=TRUE) -> fit
      c(coef(fit),loglik=as.numeric(logLik(fit)),
        conv=fit@details$convergence)
    }
  ) -> results

  attr(results,"nproc") <- getDoParWorkers()

  results
}) -> results

## ----compdest-like-------------------------------------------------------
iii <- 1-diag(length(N))

likfnCompDest <- function (theta, phi, rho, psi, tau1, tau2, delta) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  Q <- (N^tau2)*(dd^rho)
  R <- crossprod(Q,iii)^delta
  iota <- (N^tau1)*(theta*crossprod(Q*R,ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}

## ----compdest-profile-2D,cache=FALSE-------------------------------------
bake(file="compdest.rds",{

  grd <- profileDesign(
    tau1=seq(0.1, 1.4, length=25),
    tau2=seq(-1, 1, length=25),
    lower=c(psi=log(5),theta=log(0.1),phi=log(1e-6),rho=log(1.5),delta=-3),
    upper=c(psi=log(7),theta=log(45),phi=log(1e-2),rho=log(30),delta=0),
    nprof=20
  )

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau1","tau2")
      formals(likfnCompDest) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
             function(x){
               do.call(likfnCompDest,
                       c(start[fixed],setNames(as.list(x),est)))
             },
             opts=list(
               algorithm="NLOPT_LN_SBPLX",
               ftol_abs=1000,
               maxeval=10000)
      ) -> fit
      fit$solution %>% setNames(est) %>%
        c(unlist(start[fixed]),loglik=-fit$objective) %>%
        as.list() %>% as.data.frame() %>%
        cbind(conv=fit$status) -> res
    }
  ) -> results

  results %>%
    extract(sapply(results,inherits,"try-error")) %>%
    sapply(as.character) %>%
    unique() %>%
    print()

  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    ddply(~tau1+tau2,subset,loglik==max(loglik,na.rm=TRUE)) -> grd

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau1","tau2")
      formals(likfnCompDest) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
             function(x){
               do.call(likfnCompDest,
                       c(start[fixed],setNames(as.list(x),est)))
             },
             opts=list(
               algorithm="NLOPT_LN_SBPLX",
               ftol_abs=1e-5,
               xtol_rel=1e-7,
               maxeval=10000)
      ) -> fit
      fit$solution %>% setNames(est) %>%
        c(unlist(start[fixed]),loglik=-fit$objective) %>%
        as.list() %>% as.data.frame() %>%
        cbind(conv=fit$status) -> res
    }
  ) -> results

  attr(results,"nproc") <- getDoParWorkers()

  results
}) -> results

## ----rankmat,cache=FALSE-------------------------------------------------
bake(file="rankmat.rds",{
  rr <- array(dim=dim(distances),dimnames=dimnames(distances))
  for (i in seq_along(N)) {
    for (j in seq_along(N)) {
      rr[i,j] <- N[j]/(sum(N[distances[i,]<=distances[i,j]])-N[i])
    }
  }
  diag(rr) <- 0
  rr
}) -> rr

## ----stouffer-like-------------------------------------------------------
likfnStouffer <- function (theta, phi, tau1, tau2, psi) {
  theta <- exp(theta)
  phi <- exp(phi)
  iota <- (N^tau1)*(theta*((rr^tau2)%*%ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}

## ----stouffer-profile-2D,cache=FALSE-------------------------------------
bake(file="stouffer.rds",{

  grd <- expand.grid(
    tau1=seq(0.5, 1.2, length=25),
    tau2=seq(0.5, 2.0, length=25),
    theta=log(0.2),
    phi=log(0.0001),
    psi=log(5)
  )

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages="bbmle"
  ) %dopar% try(
    {
      mle2(likfnStouffer,
           method="Nelder-Mead",
           start=as.list(start[c("theta","phi","psi")]),
           fixed=as.list(start[c("tau1","tau2")]),
           control=list(trace=0,maxit=10000),
           skip.hessian=TRUE) -> fit
      c(coef(fit),loglik=as.numeric(logLik(fit)),conv=fit@details$convergence)
    }
  ) -> results

  attr(results,"nproc") <- getDoParWorkers()

  results
}) -> results

## ----rankmat1,cache=FALSE------------------------------------------------
bake(file="rankmat1.rds",{
  rr <- array(dim=dim(distances),dimnames=dimnames(distances))
  for (i in seq_along(N)) {
    for (j in seq_along(N)) {
      rr[i,j] <- N[j]/sum(N[distances[i,]<=distances[i,j]])
    }
  }
  diag(rr) <- 0
  rr
}) -> rr

## ----stouffer1-like------------------------------------------------------
likfnStouffer1 <- function (theta, phi, tau1, tau2, psi) {
  theta <- exp(theta)
  phi <- exp(phi)
  iota <- (N^tau1)*(theta*(rr^tau2)%*%ylag+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}

## ----stouffer1-profile-2D,cache=FALSE------------------------------------
bake(file="stouffer1.rds",{

  grd <- expand.grid(
    tau1=seq(0.5, 1.2, length=25),
    tau2=seq(0.5, 2.0, length=25),
    theta=log(0.2),
    phi=log(0.0001),
    psi=log(5)
  )

  foreach(start=iter(grd,"row"),
          .errorhandling="pass",
          .inorder=FALSE,
          .packages="bbmle"
  ) %dopar% try(
    {
      mle2(likfnStouffer1,
           method="Nelder-Mead",
           start=as.list(start[c("theta","phi","psi")]),
           fixed=as.list(start[c("tau1","tau2")]),
           control=list(trace=0,maxit=10000),
           skip.hessian=TRUE) -> fit
      c(coef(fit),loglik=as.numeric(logLik(fit)),conv=fit@details$convergence)
    }
  ) -> results

  attr(results,"nproc") <- getDoParWorkers()

  results
}) -> results


## ----radmat,cache=FALSE--------------------------------------------------
bake(file="radmat.rds",{
  rr <- array(dim=dim(distances),dimnames=dimnames(distances))
  for (i in seq_along(N)) {
    for (j in seq_along(N)) {
      s <- sum(N[distances[i,]<=distances[i,j]])
      rr[i,j] <- N[i]*N[j]*N[j]/s/(s-N[i])
    }
  }
  diag(rr) <- 0
  rr
}) -> rr

## ----radmat1,cache=FALSE-------------------------------------------------
bake(file="radmat1.rds",{
  rr <- array(dim=dim(distances),dimnames=dimnames(distances))
  for (i in seq_along(N)) {
    for (j in seq_along(N)) {
      s <- sum(N[distances[i,]<=distances[i,j]])
      rr[i,j] <- N[i]*N[j]*N[j]/(s+N[i])/(s)
    }
  }
  diag(rr) <- 0
  rr
}) -> rr

## ----session-info-------------------------------------------------------
sessionInfo()
