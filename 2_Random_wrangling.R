

library(raster)
library(sp)
library(sp)
library(uuid)

## Defining some coordinate systems

crs_temp_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs_temp_utm33 <- CRS("+proj=utm +zone=33 +datum=WGS84")

################################################################################
################################################################################
### Creating space states (i.e. the sampling grid)
### The coordinates for the orignal grids are found in tibble "d2"

### Using the original raster grid; 
d2_temp <- d2 %>% filter(id!=100)

## Re-creating the sampling grid (Space state)

Space_state_large <- raster(ncol=5, nrow=5, xmn=min(d2_temp$min_x), 
                                        xmx=max(d2_temp$max_x), 
                                        ymn=min(d2_temp$min_y), 
                                        ymx=max(d2_temp$max_y), 
                                        crs=crs_temp_utm33)

## Re-creating original cell numbers (which is not consecutive)
values(Space_state_large) <- c(seq(21,25), seq(1:20))
plot(Space_state_large)

###################################################################
###################################################################
### A small (250 X 250 m) grid; 

Space_state_small <- raster(ncol=5*4, nrow=5*4, xmn=min(d2_temp$min_x), 
                            xmx=max(d2_temp$max_x), 
                            ymn=min(d2_temp$min_y), 
                            ymx=max(d2_temp$max_y), 
                            crs=crs_temp_utm33)

## Re-creating original cell numbers (which is not consequtive)
values(Space_state_small) <- 1:ncell(Space_state_small)
plot(Space_state_small)

###################################################################
###### Getting the center point of the grid cell = detectors

detectors_large <- tibble(id=seq(1:ncell(Space_state_large)), cell_nr= c(seq(21,25), seq(1:20)), 
                     x=xyFromCell(Space_state_large,1:ncell(Space_state_large))[,1],
                     y=xyFromCell(Space_state_large,1:ncell(Space_state_large))[,2])

detectors_small <- tibble(id=seq(1:ncell(Space_state_small)), cell_nr= seq(1:ncell(Space_state_small)), 
                          x=xyFromCell(Space_state_small,1:ncell(Space_state_small))[,1],
                          y=xyFromCell(Space_state_small,1:ncell(Space_state_small))[,2])

####################################################################
###### Creating a spatial points object; 

detectors_small_spldf <- SpatialPointsDataFrame(coords=cbind(detectors_small$x,detectors_small$y), proj4string=crs_temp_utm33, 
                                                data=tibble(cell_nr=detectors_small$cell_nr), match.ID=TRUE)

detectors_small_spldf_latlong <- sp::spTransform(detectors_small_spldf, crs_temp_longlat) 
detectors_small_spldf_latlong@data$locationID <- str_c("cellNR_", detectors_small_spldf_latlong@data$cell_nr)

##### And - creating a tibble with lat-long and locationID

detectors_small_LatLong <- tibble(locationID=detectors_small_spldf_latlong@data$locationID, 
                                  decimalLatitude=detectors_small_spldf_latlong@coords[,2], 
                                  decimalLongitude=detectors_small_spldf_latlong@coords[,1])

####################################################################
#### Plotting the detectors: 

pp <- ggplot(detectors_small, aes(x=x, y=y)) +
      geom_point(colour="orange", pch=16) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(colour="none") +
  xlab("UTM easting") +
  ylab("UTM northing")

pp

###################################################################
###################################################################
###### 

### Starting to prepare occurrence table; 

occurrences <- d %>% mutate(eventDate=lubridate::date(date)) %>%
                  mutate(Year=year(eventDate), Month=month(eventDate),
                  Day=day(eventDate)) %>%
                  filter(!is.na(ind_id)) %>%
                  mutate(Season=if_else(Month<5, "winter", "spring"), 
                         Session=paste(Year, Season, sep="_")) %>%
                  mutate(organismName=ind_id, 
                         scientificName=if_else(species_dna=="Fjellrype", "Lagopus muta", 
                                                "Lagopus lagopus")) %>%
                  filter(!is.na(x)) %>%
                  mutate(taxonID=if_else(scientificName=="Lagopus muta", "https://www.biodiversity.no/4070", "https://www.biodiversity.no/4066"), 
                  kingdom="Animalia", phylum="Chordata", family="Phasianidae", genus="Lagopus", speciesEpithet=if_else(scientificName=="Lagopus muta", "muta", "lagopus"), 
                  lifeStage="unknown", occurrenceID=uuid, recordNumber=sample_id, catalogNumber="Genlab_ID", recordedByID="https://orcid.org/0000-0002-5119-8331", 
                  individualCount=1, occurrenceStatus="Present", organismID="Create a unique ID", basisOfRecord="MaterialSample") 
                  
### Overaly trap grid (Space_state_small) to obtain trap ID for each "capture"
## Creating sp_SpatialPointsDataFrame
temp_points <-  SpatialPointsDataFrame(coords=cbind(occurrences$x,occurrences$y), proj4string=crs_temp_utm33, 
                                       data=as.data.frame(occurrences$uuid), match.ID=TRUE)
## Extract "trapID"
occurrences <- occurrences %>% mutate(locationID_small=extract(Space_state_small, temp_points), 
                                           locationID_large=extract(Space_state_large, temp_points))

## Append decimalLatitude og decimalLongitude to "occurrences"

temp_points2 <- spTransform(temp_points, crs_temp_longlat)
occurrences <- occurrences %>% mutate(decimalLatitude=temp_points2@coords[,2], 
                                      decimalLongitude=temp_points2@coords[,1])

###########################################################################
###########################################################################
### Plotting - rock ptarmigan only

d_r <- occurrences %>% filter(scientificName=="Lagopus muta") %>%
        mutate(Session2=factor(Session, levels=c("2014_winter", "2014_spring",
                                                 "2015_winter", "2015_spring", 
                                                 "2016_winter", "2016_spring", 
                                                 "2017_winter", "2017_spring")))

p2 <- ggplot(data=d_r, aes(x=x, y=y, colour=Session)) +
      geom_point(size=1, alpha=0.7) +
     # geom_point(data=detectors_small, aes(x=x, y=y)) +
      facet_wrap(~Session2) +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(colour="none") +
  xlab("UTM easting") +
  ylab("UTM northing")

p2

##############################################################
##############################################################
## Number of samples pr. individual

temp_n <- d_r %>% group_by(organismName) %>%
          summarise(n_samples=n()) %>%
          arrange(., desc(n_samples))


#############################################################
#############################################################
## A simple plot to look at locations for each bird

temp3 <- temp_n %>% filter(n_samples>6) %>%
          left_join(., d_r)


p2 <- ggplot(data=temp3, aes(x=x, y=y, colour=Session)) +
  geom_point(size=2, alpha=0.5) +
  # geom_point(data=detectors_small, aes(x=x, y=y)) +
  facet_wrap(~organismName) +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(colour="none") +
  xlab("UTM easting") +
  ylab("UTM northing") 

p2

         
#############################################################################################
#############################################################################################
#############################################################################################
##### Creating track log - with one segment for each visit in a cell
##### on track-log might have several visits in one cell
         
track_log_points <- SpatialPoints(coords=cbind(d4$x,d4$y), proj4string=crs_temp_longlat)  
track_log_points <- sp::spTransform(track_log_points, crs_temp_utm33)         

##### Extract grid ID values for each track point, and add to data set; 

track_log <- d4 %>% mutate(locationID_small=extract(Space_state_small, track_log_points), 
                           locationID_large=extract(Space_state_large, track_log_points)) %>%
                    filter(!is.na(locationID_small)) %>%
                    mutate(track_box_day_id=paste(transect_id, locationID_small, sep="_")) %>%
                    arrange(transect_id, id) 


########################################
##### Then we need to identify unique tracks. I think this will do the trick.

track_box_day_id <- track_log %>% dplyr::select(track_box_day_id)

track_seg_point_id_new <- transect_id_new <- numeric()
track_seg_point_id_new[1] <- transect_id_new[1] <- 1                                                                   

for (i in 2:dim(track_box_day_id)[1]){
  track_seg_point_id_new[i] <- ifelse(track_box_day_id[i,1]==track_box_day_id[i-1,1], track_seg_point_id_new[i-1]+1, 1)
  transect_id_new[i] <- ifelse(track_box_day_id[i,1]==track_box_day_id[i-1,1], transect_id_new[i-1], transect_id_new[i-1]+1)
}


track_log <- bind_cols(list(track_log, track_seg_point_id_new=track_seg_point_id_new, transect_id_new=transect_id_new))

##########################################
### Calculating length of the new track log segments pr grid cell; 

## This is a placeholder for the line lengths
Line_lengths <- tibble(transect_id_new=seq(1,length(unique(track_log$transect_id_new))), line_length=numeric(length=length(unique(track_log$transect_id_new))))

### Filtering - using only segments with >1 point 
Track_lines <- track_log %>% group_by(transect_id_new) %>%
                summarize(n_points=n()) %>%
                filter(n_points>1) %>%
                dplyr::select(transect_id_new)

#### Looping over each track log segment - this takes a bit of time
#### Could/should be made faster, but it works.  

for(i in 1:dim(Track_lines)[1]){
temp1 <- track_log %>% filter(transect_id_new==Track_lines$transect_id_new[i])
temp_coord <- cbind(temp1$x, temp1$y)
line_temp <- spLines(temp_coord, crs=crs_temp_longlat)
line_temp2 <- spTransform(line_temp, crs_temp_utm33)
Line_lengths[Track_lines$transect_id_new[i],2] <- rgeos::gLength(line_temp2)
              }
                  
################################################################################
#### Adding this to a new data set - with one row for each track segment + date + grid cell id
#### This will be the secondary sessions in the event table. 

 effort_secondary <- track_log %>% dplyr::select(transect_id_new, datum, locationID_small, track_seg_point_id_new) %>%
          group_by(transect_id_new) %>%
          slice(which.min(track_seg_point_id_new)) %>%
          ungroup() %>%
           left_join(., Line_lengths) %>%
           group_by(datum, locationID_small) %>%
          summarize(samplingEffort=sum(line_length)) %>%
          ungroup() %>%
          mutate(eventDate=lubridate::date(datum)) %>%
          mutate(Year=year(eventDate), Month=month(eventDate),
          Day=day(eventDate)) %>%
          mutate(Season=if_else(Month<5, "winter", "spring"), 
          verbatimEventDate=paste(Year, Season, sep="_")) %>%
          dplyr::select(-datum) %>%
          mutate(samplingEffort=if_else(samplingEffort==0, 1, samplingEffort)) %>%
          mutate(locationID=str_c("cellNR_", locationID_small), coordinateUncertaintyInMeters=125, eventRemarks="Secondary_session") %>%
          left_join(., detectors_small_LatLong)  

### Generating unique IDs for each secondary session - and attach to 
temp_uuid1 <- tibble(eventID=UUIDgenerate(n=dim(effort_secondary)[1])) 
effort_secondary <- bind_cols(temp_uuid1, effort_secondary)

### Adding eventID to occurrence table - joining by locationID_small and eventDate; 
test <- effort_secondary %>% dplyr::select(eventID, eventDate, locationID_small) 
occurrences <- left_join(occurrences, test) 

##############################################################################################
### summarizing effort - example by primary session (can be summarized for other time spans)
effort_primary <- effort_secondary %>% group_by(locationID, verbatimEventDate) %>%
                  dplyr::summarise(samplingEffort=sum(samplingEffort), 
                  eventDate_start=min(eventDate), eventDate_stop=max(eventDate)) 
                  #mutate(eventDate=str_c(eventDate_start, eventDate_stop, sep="/"))  
                  ## use separate(eventDate, c("start", "stop"), sep="/") to cast back

#### Completing the secondary sessions - including cells that were not surveyed and setting samplingEffort=0

temp <- tibble(verbatimEventDate=rep(unique(effort_primary$verbatimEventDate), each=400), locationID_small=rep(seq(1,400), 7)) %>%
        mutate(locationID=str_c("cellNR_", locationID_small))

#### joining with effort_primary
effort_primary <- effort_primary %>% right_join(., temp) %>%
                  mutate(samplingEffort=replace_na(samplingEffort, 0), 
                         coordinateUncertaintyInMeters=125, eventRemarks="Primary_session") %>%
                  left_join(., detectors_small_LatLong)

### Generating unique IDs for each secondary session - and attach to 
temp2_uuid <- tibble(eventID=UUIDgenerate(n=dim(effort_primary)[1])) 
effort_primary <- bind_cols(temp2_uuid , effort_primary)

### Extracting the events + locationID and verbatimEventDate, 
### as this will be the parentEventID

temp <- effort_primary %>% dplyr::select(eventID, locationID, verbatimEventDate) %>%
        rename(parentEventID=eventID)
effort_secondary <- effort_secondary  %>% left_join(., temp)

############################################################################################
############################################################################################
#### Adding primary and secondary sessions together in one event table;

primary_short <- effort_primary %>% dplyr::select(eventID, verbatimEventDate, samplingEffort, locationID, decimalLatitude, decimalLongitude, eventRemarks) %>%
                 mutate(coordinateUncertaintyInMeters=125)

secondary_short <- effort_secondary %>% dplyr::select(eventID, parentEventID, eventDate, verbatimEventDate, eventRemarks, samplingEffort, locationID, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters) 
                   

event_table <- bind_rows(secondary_short, primary_short) %>%
                mutate(country="Norway", countryCode="NO", stateProvince="Trøndelag", municipality="Lierne", locality="Lifjellet", geodeticDatum="EPSG:4326")

## To do: Add footprintWKT, footprintSRS

occurrence_table <- occurrences %>% dplyr::select(eventID, eventDate, occurrenceID, basisOfRecord, recordNumber, catalogNumber, organismName, organismID, 
                                                  scientificName, lifeStage, sex, individualCount, 
                                                  taxonID, kingdom, phylum, family, genus, speciesEpithet, 
                                                  occurrenceStatus, recordedByID)

############################################################################################
############################################################################################
### Plotting nice raster maps with effort pr small cell: 

effort_map2 <- event_table %>% filter(eventRemarks=="Primary_session") %>%
              mutate(samplingEffort=if_else(samplingEffort>500, 500, samplingEffort)) %>%
              mutate(verbatimEventDate=factor(verbatimEventDate, levels=c("2014_spring",
                                           "2015_winter", "2015_spring", 
                                           "2016_winter", "2016_spring", 
                                           "2017_winter", "2017_spring")))

## converting to UTM 33 - to create a regularly spaced grid;
##use integer-values to ensure regularity; 

temp_points <- sp::SpatialPoints(coords = cbind(effort_map2$decimalLongitude, effort_map2$decimalLatitude), proj4string=crs_temp_longlat)
temp_points <- sp::spTransform(temp_points, crs_temp_utm33)
temp <- tibble(x=as.integer(temp_points@coords[,1]), y=as.integer(temp_points@coords[,2]))
effort_map2 <- bind_cols(effort_map2, temp)

## Plotting effort by session; 

m2 <- ggplot() +
  xlab("UTM Easting") + 
  ylab("UTM Northing") +
  geom_raster(data=effort_map2, aes(x=x, y=y, fill=samplingEffort), interpolate=TRUE) +
  scale_fill_viridis_c(alpha=0.6) +
  facet_wrap(~verbatimEventDate) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill="Effort")

m2
























