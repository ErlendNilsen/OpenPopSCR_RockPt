

library(abind)
library(tidyverse)
library(reshape2)

########################################################################################
########################################################################################
## 1: Filtering the relevant data. In this case, only data from the spring sampling

event_table_spring <- event_table %>% mutate(season=verbatimEventDate) %>%
                    separate(season, c("year", "season"), sep="_") %>%
                    filter(season=="spring")

temp <- event_table_spring %>% filter (eventRemarks=="Secondary_session") %>% 
        dplyr::select(eventID, eventDate, locationID)

occurrence_table_spring <- inner_join(occurrence_table, temp)

#######################################################################################
#######################################################################################
## Create the y-array, holding the observation data
## This is a 3D array with dimension n_birds x n_dertectors x n_sessions

## Extracting the relevant columns; 
temp <- event_table_spring %>% dplyr::select(eventID, locationID, parentEventID)

##  same - for the primary sessions
temp2 <- event_table_spring %>% filter(eventRemarks=="Primary_session") %>% 
        dplyr::select(eventID, locationID, verbatimEventDate)  

## Creating tibble  from event- and occurence tables, 
## with n_detections, organismName, locationID and eventID
temp_rpt <- occurrence_table_spring %>% filter(scientificName=="Lagopus muta") %>%
            left_join(., temp) %>%
            group_by(parentEventID, organismName) %>%
            summarise(n_detections=n()) %>%
            ungroup() %>%
            mutate(eventID=parentEventID) %>%
            left_join(., temp2) %>%
            separate(locationID, c("cell", "Cell_nr"), sep="_") %>%
            mutate(cell_numerical=as.numeric(Cell_nr))

#### Expanding the tibble to include also detectors where there are no observations 
#### - also those that are not sampled; 

cell_numerical <- tibble(cell_numerical=seq(1,400))

temp_rpt2 <- full_join(cell_numerical, temp_rpt) %>%
              filter(!is.na(cell_numerical)) %>%
              mutate(verbatimEventDate=factor(verbatimEventDate, levels=c("2014_spring",
                                                              "2015_winter", "2015_spring", 
                                                              "2016_winter", "2016_spring", 
                                                              "2017_winter", "2017_spring")))

###### Creating y.alive; 
CH2 <- reshape2::acast(temp_rpt2, organismName~cell_numerical~verbatimEventDate, 
                       value.var="n_detections", fun.aggregate=sum)
###### 

y.alive <- CH2[,,1:4]

## Adding augmented birds, and adding to the y-array; 

ind.augmented <- 1.2
y.augmented <- array(0, dim=c((dim(y.alive)[1]*ind.augmented), dim(y.alive)[2], dim(y.alive)[3]))

y <- abind::abind(y.alive, y.augmented, along=1)

######################################################################
### Creating z-array

z.array <- y 
z.array[z.array ==0] <-  NA
z.array[z.array >0] <- 1

################################################################################
### Creating lowerHabCoords & upperHabcoords

## The bounding box around the study area
bounding_box <- tibble(x=c(408000, 413000), y=c(7148000, 7153000)) 

lowerHabCoords <- tibble(x=rep(seq(min(bounding_box$x), max(bounding_box$x)-250, length.out = 20), 20), 
                         y=rep(seq(max(bounding_box$y)-250, min(bounding_box$y), length.out = 20), each=20))


upperHabCoords <- tibble(x=rep(seq(min(bounding_box$x)+250, max(bounding_box$x), length.out = 20), 20), 
                         y=rep(seq(max(bounding_box$y), min(bounding_box$y)+250, length.out = 20), each=20))

################################################################################
#### detCovs; we use track logs for each detector each session

temp_detCovs <- event_table_spring %>% filter(eventRemarks=="Primary_session")%>%
  mutate(verbatimEventDate=factor(verbatimEventDate, levels=c("2014_spring",
                                                              "2015_winter", "2015_spring", 
                                                              "2016_winter", "2016_spring", 
                                                              "2017_winter", "2017_spring"))) %>%
  separate(locationID, c("cell", "Cell_nr"), sep="_") %>%
  mutate(cell_numerical=as.numeric(Cell_nr))
  

detCov <- reshape2::acast(temp_detCovs, cell_numerical~verbatimEventDate, 
                       value.var="samplingEffort", fun.aggregate=sum)


################################################################################
################################################################################
#### detector.xy

temp_xy <- event_table_spring %>% filter(eventRemarks=="Primary_session" & verbatimEventDate=="2014_spring") 

temp_xy_sp <- SpatialPoints(coords=cbind(temp_xy$decimalLongitude,temp_xy$decimalLatitude), 
                            proj4string=crs_temp_longlat)

temp_xy_sp2 <- sp::spTransform(detectors_small_spldf, crs_temp_utm33) 

detector.xy <- tibble(locationID=temp_xy$locationID, x=temp_xy_sp2@coords[,1], 
                      y=temp_xy_sp2@coords[,2]) %>%
  separate(locationID, c("cell", "Cell_nr"), sep="_") %>%
  mutate(cell_numerical=as.numeric(Cell_nr)) %>%
  arrange(cell_numerical)


################################################################################
################################################################################
### habitatIDDet

habitatIDDet <- matrix(seq(1,400), ncol=20, byrow=T)

################################################################################
################################################################################
### Making the sparseY array
sparseY <- getSparseY(y=y)
yDets <- sparseY$yDets
nbDetections <- sparseY$nbDetections










    
            