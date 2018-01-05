library(aakmisc)
library(readr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringi)

options(aakmisc.dbname="us_commute",
	aakmisc.remotehost="kinglab.eeb.lsa.umich.edu",
	aakmisc.user="kingaa",
  stringsAsFactors=FALSE)

## county data
read_csv("data/counties_vivek.csv",col_types="iiddic") %>%
  select(fips=FIPS.id,lat,lon=long,pop=pop_wkr,name=long.name) %>%
  mutate(fips=sprintf("%05d",fips),
    st=stri_replace_first_regex(name,"(.+)\\s([A-Z]{2})","$2"),
    co=stri_replace_first_regex(name,"(.+)\\s([A-Z]{2})","$1")) %>%
  select(fips,name,co,st,lat,lon,pop) %>%
  arrange(fips) -> y

read_csv("data/workmatrix.csv") %>%
  gather(dest,count,-Counties) %>%
  select(src=Counties,dest,count) %>%
  left_join(y %>% select(src_fips=fips,src=name),by="src") %>%
  left_join(y %>% select(dest_fips=fips,dest=name),by="dest") %>%
  select(src_fips,dest_fips,count) %>%
  arrange(src_fips,dest_fips) -> x

## distances between county centroids
library(geosphere)
distm(y %>% select(lon,lat), y %>% select(lon,lat)) %>%
  matrix(nrow=nrow(y),ncol=nrow(y),
    dimnames=list(src_fips=y$fips,dest_fips=y$fips)) %>%
  melt(value.name="distance",as.is=TRUE) %>%
  right_join(x,by=c("src_fips","dest_fips")) -> x

startTunnel()
writeDBTable("counties",y)
writeDBTable("flows",x)
stopTunnel()

x %>% saveRDS("us_commuter_flows.rds")
y %>% saveRDS("us_counties.rds")
