

library(raster)
library(sp)
library(sp)

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
                  filter(!is.na(x))
                  
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

Effort <- track_log %>% dplyr::select(transect_id_new, datum, locationID_small, track_seg_point_id_new) %>%
          group_by(transect_id_new) %>%
          slice(which.min(track_seg_point_id_new)) %>%
          ungroup() %>%
          left_join(., Line_lengths) %>%
          mutate(eventDate=lubridate::date(datum)) %>%
          mutate(Year=year(eventDate), Month=month(eventDate),
          Day=day(eventDate)) %>%
          mutate(Season=if_else(Month<5, "winter", "spring"), 
          Session=paste(Year, Season, sep="_"))%>%
          dplyr::select(-track_seg_point_id_new, -datum)

### summarizing effort - example by primary session (can be summarized for other time spans)
effort_session <- Effort %>% group_by(locationID_small, Session) %>%
                  dplyr::summarise(tl_effort=sum(line_length))


############################################################################################
############################################################################################
### Plotting nice raster maps with effort pr small cell: 

effort_raster <- Space_state_small

temp <- tibble(locationID_small=seq(1,ncell(Space_state_small)))
temp2 <- unique(effort_session$Session)

effort_map <- tibble(locationID_small=numeric(), Session=character(), tl_effort=numeric())

## by session:
for(i in 1:length(temp2)){
effort_temp <- effort_session %>% filter(Session==temp2[i]) %>% 
               right_join(.,temp) %>%
               arrange(locationID_small) %>%
               mutate(tl_effort=replace_na(tl_effort, 0), Session=replace_na(Session, temp2[i]))
effort_map <- bind_rows(effort_map, effort_temp)
                        }
 
effort_map <- effort_map %>% mutate(cell_nr=locationID_small) %>% 
              left_join(., detectors_small) %>%
              mutate(Effort=if_else(tl_effort<500, tl_effort, 500)) %>%
              mutate(Session2=factor(Session, levels=c("2014_winter", "2014_spring",
                                           "2015_winter", "2015_spring", 
                                           "2016_winter", "2016_spring", 
                                           "2017_winter", "2017_spring")))


## Plotting effort by session; 

m2 <- ggplot() +
  xlab("UTM - easting") + 
  ylab("UTM - northing") +
  geom_raster(data=effort_map, aes(x=x, y=y, fill=Effort), interpolate=TRUE) +
  scale_fill_viridis_c(alpha=0.7) +
  facet_wrap(~Session2) +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill="Effort")

m2





















