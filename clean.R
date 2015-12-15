##' Clean and process measles and demography data.

require(plyr)
require(reshape2)
require(magrittr)
options(stringsAsFactors=FALSE)

##' Source the original datafiles

source("data/ewbu4464.q")
source("data/ewpu4464.q")
source("data/ewmu4464.q")

##' Process the measles data
names(dimnames(ewMu4464)) <- c("time","town")
ewMu4464 %>%
  melt(value.name="cases") %>%
  mutate(town=as.character(town),time=1900+time,
         cases=as.integer(cases)) -> measles

##' Process the births data
names(dimnames(ewBu4464)) <- c("year","town")
ewBu4464 %>%
  melt(value.name="births") %>%
  mutate(town=as.character(town),year=1900+year,
         births=as.integer(births)) -> births

##' Process the population-size data
names(dimnames(ewPu4464)) <- c("year","town")
ewPu4464 %>%
  melt(value.name="pop") %>%
  mutate(town=as.character(town),year=1900+year,
         pop=as.integer(pop)) -> pop

##' Join annual data into a single data frame.
join(births,pop,by=c("year","town")) -> demog

##' Write results to file.
save(measles,demog,file="ew_measles_data.rda",compress='xz')
