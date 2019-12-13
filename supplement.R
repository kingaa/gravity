params <-
list(prefix = "coupling")

## ----prelims,cache=FALSE,purl=TRUE---------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  knitr.kable.NA="",
  encoding="UTF-8",
  aakmisc.dbname="ewmeasles",
  aakmisc.remotehost="kinglab.eeb.lsa.umich.edu",
  aakmisc.user="gravity")

set.seed(594709947L)

library(grid)
library(plyr)
library(reshape2)
library(magrittr)
library(foreach)
library(iterators)
library(bbmle)
library(nloptr)
library(pomp)
library(ggplot2)
library(scales)
library(sp)
library(ncf)
theme_set(theme_bw())


## ----clean-data,purl=TRUE------------------------------------------------
stew(file="clean_data.rda",{
  library(aakmisc)
  
  startTunnel()
  
  getQuery("select town,year,births,pop from demog where year>=1944 order by town,year") -> demog
  
  getQuery("select town,date,cases from measles where year>=1944 order by town,date") %>%
    mutate(year=as.integer(format(date+3,"%Y"))) %>%
    ddply(~town+year,mutate,week=seq_along(year),biweek=(week+1)%/%2) %>%
    subset(week<=52,select=-c(date,week)) %>%
    acast(town~year~biweek,value.var="cases",fun.aggregate=sum) %>%
    melt(varnames=c("town","year","biweek"),value.name="cases") %>%
    mutate(town=as.character(town)) %>%
    arrange(town,year,biweek) %>%
    join(demog,by=c("town","year")) %>%
    mutate(births=births/26) -> dat
  
  stopTunnel()
})




## ----coords,purl=TRUE----------------------------------------------------
bake(file="coords.rds",{
  library(aakmisc)
  startTunnel()
  getQuery("select * from coords") %>% arrange(town) -> coords
  stopTunnel()
  coords
}) -> coords


## ----distances,purl=TRUE-------------------------------------------------
bake(file="distances.rds",{
  library(aakmisc)
  library(geosphere)
  
  distm(
    coords %>% subset(select=c(long,lat)),
    coords %>% subset(select=c(long,lat))) %>%
    matrix(nrow=nrow(coords),ncol=nrow(coords),
           dimnames=list(coords$town,coords$town))
}) -> distances


## ----startDoMC,cache=FALSE,include=FALSE,purl=TRUE-----------------------
  library(doMC)
  registerDoMC()


## ----under-reporting,purl=TRUE-------------------------------------------
bake(file="under-reporting.rds",{
  foreach (d=dlply(dat,~town),
           .combine=rbind,.inorder=FALSE,
           .packages=c("plyr")
           ) %dopar% {
    cumbirths <- cumsum(d$births)
    cumcases <- cumsum(d$cases)
    fit <- smooth.spline(cumbirths,cumcases,df=2.5)
    mutate(d,
           ur=predict(fit,x=cumbirths,deriv=1)$y,
           I=cases/ur)
  } -> dat
}) -> dat


## ----susceptible-reconstruction,purl=TRUE--------------------------------
bake(file="susc-reconst.rds",{
  foreach (d=dlply(dat,~town),
           .combine=rbind,.inorder=FALSE,
           .packages=c("plyr")) %dopar% {
             cumbirths <- cumsum(d$births)
             cuminc <- cumsum(d$I)
             fit <- smooth.spline(cumbirths,cuminc,df=2.5)
             mutate(d,z=-residuals(fit))
           } -> dat
}) -> dat


## ----make-data-matrix,purl=TRUE------------------------------------------
bake(file="data-matrix.rds",{
  dat %>%
    ddply(~town,mutate,
          Ilag=c(NA,head(I,-1)),
          zlag=c(NA,head(z,-1)),
          Ps=median(pop),
          seas=factor(biweek,levels=1:26)) %>%
    ddply(~town,tail,-1) -> dat
}) -> dat


## ----exclude-weird-towns,purl=TRUE---------------------------------------
bake(file="exclusions.rds",{
  dat %>%
    ddply(~town,summarize,
          maxfrac=max(-z/Ps)) %>%
    subset(maxfrac>0.03) %>%
    extract2("town") -> excludes
}) -> excludes

print(length(excludes))


## ----sigma-profile,purl=TRUE---------------------------------------------
sigma1dev <- function (dat, sigma) {
  slag <- with(dat,sigma*Ps+zlag)
  fit <- glm(log(I)~-1+seas+log(Ilag)+I(1/Ilag)+offset(log(slag)),
             data=dat,subset=Ilag>0&I>0)
  fit$deviance
}

bake(file="sigma-profile.rds",{
  dat %>% subset(!(town%in%excludes)) -> dat1
  foreach (
    sigma=seq(0.03,0.25,length=100),
    .combine=rbind,.inorder=FALSE,
    .packages=c("plyr","magrittr")) %dopar% {
      dat1 %>%
        daply(~town,sigma1dev,sigma=sigma,.parallel=TRUE) %>%
        sum() -> dev
      data.frame(sigma=sigma,dev=dev)
    } -> sigmaProf
  
}) -> sigmaProf


## ----sigma-bar,purl=TRUE-------------------------------------------------
bake(file="sigma-bar.rds",{
  dat %>% subset(!(town%in%excludes)) -> dat1
  fit <- optim(par=0.037,lower=0.03,upper=0.10,
               method="Brent",hessian=TRUE,
               fn=function (sigma) {
                 dat1 %>%
                   daply(~town,sigma1dev,sigma=sigma,.parallel=TRUE) %>%
                   sum()
               })
  fit$par -> sigmaBar
}) -> sigmaBar

dat %>% mutate(Slag=sigmaBar*Ps+zlag) -> dat




## ----fit-tsir,purl=TRUE--------------------------------------------------
bake(file="tsir-fits.rds",{
  coefnames <- c(sprintf("seas%d",1:26),"log(Ilag)","I(1/Ilag)")
  newcoefnames <- c(sprintf("log.beta%02d",1:26),"alpha","m.alpha")
  
  tsirfit <- function (dat) {
    glm(log(I)~-1+seas+log(Ilag)+I(1/Ilag)+offset(log(Slag)),
      data=dat,subset=Ilag>0&I>0) %>% summary() %>%
      extract2("coefficients") -> fit
    fit[,"Estimate"] %>% extract(coefnames) %>% as.list() %>% as.data.frame() %>%
      set_names(newcoefnames) -> coefs
    fit[,"Std. Error"] %>% extract(coefnames) %>% as.list() %>% as.data.frame() %>%
      set_names(paste(newcoefnames,"se",sep=".")) -> se
    cbind(coefs,se,sigma=sigmaBar,
      town=unique(dat$town),Ps=unique(dat$Ps))
  }
  
  dat %>% ddply(~town,tsirfit,.parallel=TRUE) %>%
    melt(id=c("town","Ps")) %>%
    mutate(se=ifelse(grepl("\\.se$",variable),"se","est"),
      variable=sub("\\.se$","",variable)) %>%
    dcast(town+Ps+variable~se) -> tsirs
  
  dat %>% ddply(~town,summarize,Ps=unique(Ps)) -> tsircoef
  
  tsirs %>%
    na.omit() %>%
    ddply(~variable, function (d) {
      fit <- lm(est~log(Ps),data=d,weights=1/se^2)
      data.frame(town=tsircoef$town,value=predict(fit,newdata=tsircoef))
    },.parallel=TRUE) %>%
    dcast(town~variable) %>%
    melt(id=c("town","alpha","m.alpha"),value.name="log.beta") %>%
    mutate(biweek=as.integer(sub("log.beta","",as.character(variable)))) %>%
    arrange(town,biweek) %>%
    subset(select=-variable) -> tsircoef
  
  dat %>%
    join(tsircoef,by=c("town","biweek")) %>%
    mutate(ylag=Ilag/Ps) -> dat
}) -> dat




## ----gravity-pre,purl=TRUE-----------------------------------------------
readRDS("tsir-fits.rds") -> dat
dat %>% acast(town~year+biweek,value.var="ylag") -> ylag
dat %>% acast(town~year+biweek,value.var="Slag") -> slag
dat %>% acast(town~year+biweek,value.var="log.beta") %>% exp() -> beta
dat %>% acast(town~year+biweek,value.var="cases") -> obs
dat %>% daply(~town,function(x)unique(x$Ps)) -> N

readRDS("distances.rds") -> distances
dd <- 1/distances
diag(dd) <- 0
dd <- dd[rownames(ylag),rownames(ylag)]

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


## ----startMPI,cache=FALSE,include=FALSE,purl=TRUE------------------------
if (file.exists("CLUSTER")) {
  scan("CLUSTER",what=integer(0)) -> ncpu
  library(doMPI)
  cl <- startMPIcluster(ncpu,verbose=TRUE,logdir="/tmp")
  registerDoMPI(cl)
} else {
  library(doMC)
  registerDoMC()
}


## ----gravity-like,purl=TRUE----------------------------------------------
likfnGravity <- function (theta, phi, rho, psi, tau1, tau2) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  Q <- (N^tau2)*(dd^rho)
  iota <- (N^tau1)*(theta*crossprod(Q,ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----gravity-profile-2D,cache=FALSE,purl=TRUE----------------------------
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


## ----gravity-mle,cache=FALSE---------------------------------------------
bake("gravity-mle.rds",{
  
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnGravity,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.grav






## ----xia-like,purl=TRUE--------------------------------------------------
likfnXia <- function (theta, phi, rho, psi, tau1, tau2) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  Q <- (N^tau2)*(dd^rho)
  iota <- (N^tau1)*(theta*crossprod(Q,ylag^tau2)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----xia-profile-2D,cache=FALSE,purl=TRUE--------------------------------
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


## ----xia-mle,cache=FALSE-------------------------------------------------
bake("xia-mle.rds",{
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnXia,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.xia






## ----meanfield-like,purl=TRUE--------------------------------------------
likfnMeanField <- function (theta, phi, psi) {
  theta <- exp(theta)
  phi <- exp(phi)
  Q <- N*dd
  iota <- N*(theta*crossprod(Q,ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----meanfield-optim,cache=FALSE,purl=TRUE-------------------------------
bake(file="meanfield.rds",{

  grd <- sobolDesign(
    lower=c(psi=log(5),theta=log(1e-8),phi=log(1e-6)),
    upper=c(psi=log(7),theta=log(1),phi=log(1e-2)),
    nseq=250
  )
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c()
      formals(likfnMeanField) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
        function(x){
          do.call(likfnMeanField,
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
    subset(loglik==max(loglik,na.rm=TRUE)) -> grd
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c()
      formals(likfnMeanField) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
        function(x){
          do.call(likfnMeanField,
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




## ----meanfield-mle,cache=FALSE-------------------------------------------
bake("meanfield-mle.rds",{
  
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnMeanField,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.meanfield




## ----diffusion-like,purl=TRUE--------------------------------------------
likfnDiffusion <- function (theta, phi, rho, psi) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  Q <- dd^rho
  iota <- (theta*crossprod(Q,ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----diffusion-profile-2D,cache=FALSE,purl=TRUE--------------------------
bake(file="diffusion.rds",{

  grd <- profileDesign(
    rho = seq(log(1),log(3),length=50),
    lower=c(psi=log(5),theta=log(1e-8),phi=log(1e-6)),
    upper=c(psi=log(7),theta=log(1),phi=log(1e-2)),
    nprof=10
  )
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("rho")
      formals(likfnDiffusion) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
        function(x){
          do.call(likfnDiffusion,
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
    ddply(~rho,subset,loglik==max(loglik,na.rm=TRUE)) -> grd
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("rho")
      formals(likfnDiffusion) %>% names() %>% setdiff(fixed) -> est
      nloptr(unlist(start[est]),
        function(x){
          do.call(likfnDiffusion,
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


## ----diffusion-mle,cache=FALSE-------------------------------------------
bake("diffusion-mle.rds",{
  
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnDiffusion,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.diffusion






## ----compdest-like,purl=TRUE---------------------------------------------
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


## ----compdest-profile-2D,cache=FALSE,purl=TRUE---------------------------
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


## ----compdest-mle,cache=FALSE--------------------------------------------
bake("compdest-mle.rds",{
  
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnCompDest,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.compdest


## ----compdest-results,purl=TRUE------------------------------------------
results %>%
  extract(sapply(results,inherits,"try-error")) %>%
  sapply(as.character) %>%
  unique()

results %>%
  extract(!sapply(results,inherits,"try-error")) %>%
  ldply() %>%
  mutate(rho=exp(rho),theta=exp(theta),psi=exp(psi),phi=exp(phi)) -> results

results %>% count(~conv)






## ----compdest-profile-delta,cache=FALSE,purl=TRUE------------------------
bake(file="compdest_delta.rds",{
  
  grd <- profileDesign(
    delta=seq(-1.6,0,length=250),
    lower=c(tau1=0.1,tau2=-1,psi=log(5),theta=log(0.1),phi=log(1e-6),rho=log(1.5)),
    upper=c(tau1=1.4,tau2=1,psi=log(7),theta=log(45),phi=log(1e-2),rho=log(30)),
    nprof=50
  )
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("delta")
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
  ) -> res
  
  res %>%
    extract(sapply(res,inherits,"try-error")) %>%
    sapply(as.character) %>%
    unique() %>%
    print()
  
  res %>%
    extract(!sapply(res,inherits,"try-error")) %>%
    ldply() %>%
    ddply(~delta,subset,loglik==max(loglik,na.rm=TRUE)) -> grd
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("delta")
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
  ) -> res1
  
  attr(res1,"nproc") <- getDoParWorkers()
  
  res1
}) -> res1




## ----compdest-profile-tau1,cache=FALSE,purl=TRUE-------------------------
bake(file="compdest_tau1.rds",{
  
  grd <- profileDesign(
    tau1=seq(0.1,1.4,length=250),
    lower=c(tau2=-1,psi=log(5),theta=log(0.1),phi=log(1e-6),rho=log(1.5),delta=-3),
    upper=c(tau2=1,psi=log(7),theta=log(45),phi=log(1e-2),rho=log(30),delta=0),
    nprof=50
  )
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau1")
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
  ) -> res
  
  res %>%
    extract(sapply(res,inherits,"try-error")) %>%
    sapply(as.character) %>%
    unique() %>%
    print()
  
  res %>%
    extract(!sapply(res,inherits,"try-error")) %>%
    ldply() %>%
    ddply(~tau1,subset,loglik==max(loglik,na.rm=TRUE)) -> grd
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau1")
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
  ) -> res2
  
  attr(res2,"nproc") <- getDoParWorkers()
  
  res2
}) -> res2




## ----compdest-profile-tau2,cache=FALSE,purl=TRUE-------------------------
bake(file="compdest_tau2.rds",{
  
  grd <- profileDesign(
    tau2=seq(-1,1,length=250),
    lower=c(tau1=0.1,psi=log(5),theta=log(0.1),phi=log(1e-6),rho=log(1.5),delta=-3),
    upper=c(tau1=1.4,psi=log(7),theta=log(45),phi=log(1e-2),rho=log(30),delta=0),
    nprof=50
  )
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau2")
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
  ) -> res
  
  res %>%
    extract(sapply(res,inherits,"try-error")) %>%
    sapply(as.character) %>%
    unique() %>%
    print()
  
  res %>%
    extract(!sapply(res,inherits,"try-error")) %>%
    ldply() %>%
    ddply(~tau2,subset,loglik==max(loglik,na.rm=TRUE)) -> grd
  
  foreach(start=iter(grd,"row"),
    .errorhandling="pass",
    .inorder=FALSE,
    .packages=c("nloptr","magrittr")
  ) %dopar% try(
    {
      fixed <- c("tau2")
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
  ) -> res3
  
  attr(res3,"nproc") <- getDoParWorkers()
  
  res3
}) -> res3




## ----rankmat,cache=FALSE,purl=TRUE---------------------------------------
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


## ----stouffer-like,purl=TRUE---------------------------------------------
likfnStouffer <- function (theta, phi, tau1, tau2, psi) {
  theta <- exp(theta)
  phi <- exp(phi)
  iota <- (N^tau1)*(theta*((rr^tau2)%*%ylag)+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----stouffer-profile-2D,cache=FALSE,purl=TRUE---------------------------
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


## ----stouffer-mle,cache=FALSE--------------------------------------------
bake("stouffer-mle.rds",{
  
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnStouffer,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.stouffer






## ----rankmat1,cache=FALSE,purl=TRUE--------------------------------------
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


## ----stouffer1-like,purl=TRUE--------------------------------------------
likfnStouffer1 <- function (theta, phi, tau1, tau2, psi) {
  theta <- exp(theta)
  phi <- exp(phi)
  iota <- (N^tau1)*(theta*(rr^tau2)%*%ylag+phi)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----stouffer1-profile-2D,cache=FALSE,purl=TRUE--------------------------
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


## ----stouffer1-mle,cache=FALSE-------------------------------------------
bake("stouffer1-mle.rds",{
  
  results %>%
    extract(!sapply(results,inherits,"try-error")) %>%
    ldply() %>%
    subset(loglik==max(loglik),select=-c(loglik,conv)) -> start
  
  mle2(likfnStouffer1,
    method="Nelder-Mead",
    start=as.list(start),
    control=list(trace=0,maxit=10000),
    skip.hessian=TRUE) -> fit
  
  c(coef(fit),loglik=as.numeric(logLik(fit)),
    conv=fit@details$convergence) %>%
    as.list() %>%
    as.data.frame() -> mle
}) -> mle.stouffer1






## ----radmat,cache=FALSE,purl=TRUE----------------------------------------
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


## ----radiation-like------------------------------------------------------
likfnRadiation <- function (theta, psi) {
  theta <- exp(theta)
  iota <- theta*(rr%*%ylag)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----radiation-mle,cache=FALSE-------------------------------------------
bake(file="radiation-mle.rds",{
  mle2(likfnRadiation,
    method="Nelder-Mead",
    start=list(theta=log(1),psi=log(5)),
    control=list(trace=0,maxit=10000),
    skip.hessian=FALSE)
}) -> fit




## ----radmat1,cache=FALSE,purl=TRUE---------------------------------------
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


## ----radiation1-like-----------------------------------------------------
likfnRadiation1 <- function (theta, psi) {
  theta <- exp(theta)
  iota <- theta*(rr%*%ylag)
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----radiation1-mle,cache=FALSE------------------------------------------
bake(file="radiation1-mle.rds",{
  mle2(likfnRadiation1,
    method="Nelder-Mead",
    start=list(theta=log(1),psi=log(5)),
    control=list(trace=0,maxit=10000),
    skip.hessian=FALSE)
}) -> fit




## ----xradmat1,cache=FALSE,purl=TRUE--------------------------------------
readRDS("radmat1.rds") -> rr


## ----xrad1-like----------------------------------------------------------
likfnXRad1 <- function (theta, psi, phi) {
  theta <- exp(theta)
  phi <- exp(phi)
  iota <- theta*(rr%*%ylag)+phi*N
  lambda <- betaS*iota[relevant]^alpha
  -sum(dnbinom(x=obs,mu=lambda,size=exp(-psi),log=TRUE))
}


## ----xrad1-mle,cache=FALSE-----------------------------------------------
bake(file="xrad1-mle.rds",{
  mle2(likfnXRad1,
    method="Nelder-Mead",
    start=list(theta=log(1),psi=log(5),phi=log(0.00001)),
    control=list(trace=0,maxit=10000),
    skip.hessian=FALSE)
}) -> fit




























## ----session-info,cache=FALSE--------------------------------------------
sessionInfo()

