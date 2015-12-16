##' Clean and process measles and demography data.

require(plyr)
require(reshape2)
require(magrittr)
options(stringsAsFactors=FALSE)

# ##' Source the original datafiles

# source("data/ewbu4464.q")
# source("data/ewpu4464.q")
# source("data/ewmu4464.q")

# ##' Process the measles data
# names(dimnames(ewMu4464)) <- c("time","town")
# ewMu4464 %>%
#   melt(value.name="cases") %>%
#   mutate(town=as.character(town),time=1900+time,
#          cases=as.integer(cases)) -> measles
#
# ##' Process the births data
# names(dimnames(ewBu4464)) <- c("year","town")
# ewBu4464 %>%
#   melt(value.name="births") %>%
#   mutate(town=as.character(town),year=1900+year,
#          births=as.integer(births)) -> births
#
# ##' Process the population-size data
# names(dimnames(ewPu4464)) <- c("year","town")
# ewPu4464 %>%
#   melt(value.name="pop") %>%
#   mutate(town=as.character(town),year=1900+year,
#          pop=as.integer(pop)) -> pop
#
# ##' Join annual data into a single data frame.
# join(births,pop,by=c("year","town")) -> demog

##' Retrieve data from database
require(aakmisc)
options(aakmisc.dbname="ewmeasles",aakmisc.remotehost="kinglab.eeb.lsa.umich.edu")

startTunnel()

getQuery("select town,year,births,pop from demog where year>=1944 order by town,year") -> demog

##' Measles data is weekly.
##' Associate each week's data with the Weds. of that week.
##' Discard any 53rd weeks (1947,1952,1958,1964)
##' Aggregate by biweek (26 biweeks/yr)
##' Join with demographic data
##' Scale birth rates to births/biweek
getQuery("select town,date,cases from measles where year>=1944 order by town,date") %>%
  mutate(year=as.integer(format(date+3,"%Y"))) %>%
  ddply(~town+year,mutate,week=seq_along(year),biweek=(week+1)%/%2) %>%
  subset(week<=52,select=-c(date,week)) %>%
  acast(town~year~biweek,value.var="cases",fun.aggregate=sum) %>%
  melt(varnames=c("town","year","biweek"),value.name="cases") %>%
  mutate(town=as.character(town)) %>%
  arrange(town,year,biweek) %>% 
  join(demog,by=c("town","year")) %>%
  mutate(births=births/26) -> measles

stopTunnel()

##' Write results to file.
saveRDS(measles,file="ew_measles_data.rds",compress='xz')
