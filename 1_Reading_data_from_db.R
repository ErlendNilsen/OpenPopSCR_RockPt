
library(RPostgreSQL)
library(tidyverse)
library(lubridate)
source("0_credentials.R")

## 1: Connecting to the database on nina gisdata-server
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname=paste(dbname1), host=paste(host1), port=paste(port1), user=paste(user1), password=paste(password1))


## 2: Which tables are in the db schema
dbGetQuery(con,
           "SELECT table_name FROM information_schema.tables
                   WHERE table_schema='fjellrype'")


## 3: fetching the main data for this purpose; 
## The DNA sample data
d <- as_tibble(fetch(dbSendQuery(con, "SELECT * from fjellrype.biological_sample ;"), -1))
d1c <- d %>% distinct(transect_id)


## The sampling grid
d2 <- as_tibble(fetch(dbSendQuery(con, "SELECT * from fjellrype.box_area ;"), -1))

## The track log data - compiled by grid cell nr
d3 <- as_tibble(fetch(dbSendQuery(con, "SELECT * from fjellrype.track_log ;"), -1)) 
d3b <- d3 %>%  distinct(filename)

## The track log data - compiled by grid cell nr. 
## This file is LARGE - so might need to loop or something. 
d4 <- as_tibble(fetch(dbSendQuery(con, "SELECT * from fjellrype.track_log_points ;"), -1))
d4b <- d4 %>% distinct(filename)
d4c <- d4 %>% distinct(transect_id)

dbDisconnect(con)

################################################################################################################
### Data wrangling: 

### Data wrangling session get started at..... "2_Random_wrangling.R"






