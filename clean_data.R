## Data preprocessing

library(plyr)
library(reshape2)
library(magrittr)

options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  dplyr.summarise.inform=FALSE,
  encoding="UTF-8",
  aakmisc.dbname="ewmeasles",
  aakmisc.remotehost="kinglab.eeb.lsa.umich.edu",
  aakmisc.user="kingaa"
)

library(aakmisc)

startTunnel()

getQuery("select town,year,births,pop from demog where year>=1944 order by town,year") -> demog

getQuery("select * from coords") %>% arrange(town) -> coords

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

save(demog,coords,dat,file="clean_data.rda")
