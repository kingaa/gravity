options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  knitr.kable.NA="",
  encoding="UTF-8",
  aakmisc.dbname="ewmeasles",
  aakmisc.remotehost="kinglab.eeb.lsa.umich.edu",
  aakmisc.user="gravity")

set.seed(594709947L)

library(plyr)
library(reshape2)
library(magrittr)
library(foreach)
library(iterators)
library(pomp)
library(doMC)
registerDoMC()

## ----clean-data,cache=FALSE----------------------------------------------
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

## ----under-reporting-----------------------------------------------------
foreach (d=dlply(dat,~town),.combine=rbind,.inorder=FALSE) %dopar% {
  cumbirths <- cumsum(d$births)
  cumcases <- cumsum(d$cases)
  fit <- smooth.spline(cumbirths,cumcases,df=2.5)
  mutate(d,
         ur=predict(fit,x=cumbirths,deriv=1)$y,
         I=cases/ur)
} -> dat

## ----susceptible-reconstruction------------------------------------------
foreach (d=dlply(dat,~town),.combine=rbind,.inorder=FALSE) %dopar% {
  cumbirths <- cumsum(d$births)
  cuminc <- cumsum(d$I)
  fit <- smooth.spline(cumbirths,cuminc,df=2.5)
  mutate(d,z=-residuals(fit))
} -> dat

## ----make-data-matrix----------------------------------------------------
dat %>%
  ddply(~town,mutate,
        Ilag=c(NA,head(I,-1)),
        zlag=c(NA,head(z,-1)),
        Ps=median(pop),
        seas=factor(biweek,levels=1:26)) %>%
  ddply(~town,tail,-1) -> dat

## ----exclude-weird-towns-------------------------------------------------
dat %>%
  ddply(~town,summarize,
        maxfrac=max(-z/Ps)) %>%
  subset(maxfrac>0.03) %>%
  extract2("town") -> excludes

print(length(excludes))

## ----sigma-profile,cache=FALSE-------------------------------------------
dat %>% subset(!(town%in%excludes)) -> dat1

sigma1dev <- function (dat, sigma) {
  slag <- with(dat,sigma*Ps+zlag)
  fit <- glm(log(I)~-1+seas+log(Ilag)+I(1/Ilag)+offset(log(slag)),
             data=dat,subset=Ilag>0&I>0)
  fit$deviance
}

foreach (sigma=seq(0.03,0.25,length=100),
         .combine=rbind,.inorder=FALSE) %dopar%
         {
           dat1 %>%
             daply(~town,sigma1dev,sigma=sigma,.parallel=TRUE) %>%
             sum() -> dev
           data.frame(sigma=sigma,dev=dev)
         } -> sigmaProf

sigmaProf %>% saveRDS("sigma-profile.rds")

fit <- optim(par=0.037,lower=0.03,upper=0.10,
             method="Brent",hessian=TRUE,
             fn=function (sigma) {
               dat1 %>%
                 daply(~town,sigma1dev,sigma=sigma,.parallel=TRUE) %>%
                 sum()
             })
fit$par -> sigmaBar

rm(dat1)

sigmaBar %>% saveRDS("sigma-bar.rds")

dat %>% mutate(Slag=sigmaBar*Ps+zlag) -> dat

## ----fit-tsir,cache=FALSE------------------------------------------------
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

## ----distances,cache=FALSE-----------------------------------------------
bake(file="distances.rds",{
  library(aakmisc)

  startTunnel()
  getQuery("select * from coords") %>% arrange(town) -> coords
  stopTunnel()
  rownames(coords) <- coords$town

  library(geosphere)
  distm(coords %>% subset(select=c(long,lat)),
        coords %>% subset(select=c(long,lat))) %>%
    matrix(nrow=nrow(coords),ncol=nrow(coords),
           dimnames=list(coords$town,coords$town))
}) -> distances

## ----session-info-------------------------------------------------------
sessionInfo()
