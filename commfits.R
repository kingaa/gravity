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
library(nloptr)
library(pomp)

## ----get-data------------------------------------------------------------
readRDS("us_counties.rds") -> counties
readRDS("us_commuter_flows.rds") -> flows

flows %>%
  arrange(src_fips,dest_fips) %>%
  acast(src_fips~dest_fips,value.var="distance") -> distances

rel <- (distances > 0) & (distances < 2e5) ## restrict to distances under 200km

dd <- 1/distances
diag(dd) <- 0
ddscal <- exp(mean(log(dd[lower.tri(dd)])))
dd <- dd/ddscal

flows %>%
  arrange(src_fips,dest_fips) %>%
  acast(src_fips~dest_fips,value.var="count") %>%
  extract(rel) -> obs

counties %>%
  arrange(fips) %>%
  extract2("pop") %>%
  set_names(counties$fips) -> N

## ----gravity-like----------------------------------------------------
matGravity <- function (theta, phi, rho, tau1, tau2, ...) {
  theta <- exp(theta)
  phi <- exp(phi)
  rho <- exp(rho)
  J <- tcrossprod(N^tau1,N^tau2)
  phi*N^tau1 + theta*J*dd^rho
}

likfnGravity <- function (theta, phi, rho, psi, tau1, tau2) {
  flow <- matGravity(theta,phi,rho,tau1,tau2)
  -sum(dnbinom(x=obs,mu=flow[rel],size=exp(-psi),log=TRUE))
}

## ----gravity-profile-2D,cache=FALSE--------------------------------------
bake(file="commgrav.rds",{

  grd <- profileDesign(
    tau1=seq(0.1, 1.0, length=25),
    tau2=seq(0.3, 1.2, length=25),
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

## ----rankmat,cache=FALSE-------------------------------------------------
bake(file="commrankmat.rds",{
  rr <- array(dim=dim(distances),dimnames=dimnames(distances))
  for (i in seq_along(N)) {
    for (j in seq_along(N)) {
      rr[i,j] <- N[j]/(sum(N[distances[i,]<=distances[i,j]])-N[i])
    }
  }
  diag(rr) <- 0
  rr
}) -> rr

## ----rankmat1,cache=FALSE------------------------------------------------
bake(file="commrankmat1.rds",{
  rr <- array(dim=dim(distances),dimnames=dimnames(distances))
  for (i in seq_along(N)) {
    for (j in seq_along(N)) {
      rr[i,j] <- N[j]/sum(N[distances[i,]<=distances[i,j]])
    }
  }
  diag(rr) <- 0
  rr
}) -> rr

## ----session-info-------------------------------------------------------
sessionInfo()
