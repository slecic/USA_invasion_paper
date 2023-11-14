# USA_invasion_paper
# Instructions below contain the code needed to re-create the figures and statistics from the paper.
---
title: "wCer2Modeling"
author: "Sonja Lecic"
date: "26/1/2023"
output: 
   cleanrmd::html_document_clean:
     theme: minicss
     toc: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = F)
```

```{css}
.columns {display: flex;}
h1 {color: purple;}
h2 {color: darkorange;}
```


```{R, echo=T, include=T, message=FALSE, warning=FALSE}
# load the libraries
library(gstat)   # geostatistics
library(sf)      # spatial data handling
library(sp)       # spatial data handling
library(ggplot2) # for ploting
library(ggmap)    # maps
library(dplyr)    # data manipulation
library(tidyr)     # data manipulation
library(tidyverse) # data manipulation
library(tsibble)   # data maipulation
library(ggpubr)    # visualisation
library(geosphere) # geo maps
library(RColorBrewer)  # data visualisation
library(maptools)  # geo maps
library(geodist)   # geo maps
library(rnaturalearth) # river data
library(viridis)   # data visualisation
library(scatterpie) # data visualisation
library(circular) # coordinates
library(Hmisc) # statistics
library(reshape2) # data manipulation
```

First import and save the maps
```{R, echo=T, warning=FALSE, message=FALSE}
#####################
### Map of Europe
#####################

####### make a white and grey map
map <- get_stamenmap( bbox = c(left = -11, bottom = 35, right = 45, top = 61), zoom = 4, maptype = "toner-background")

## Rivers
# download river shape file to plot Danube on the map
#rivers110 <- ne_download(scale = 110, type = 'rivers_lake_centerlines', category = 'physical')
# transform the shape file into a data frame read by ggplot
#rivers_df <- fortify(rivers110)
# and save 
#write.csv(rivers_df, "/Volumes/LaCie/rivers_df.csv", row.names=FALSE)
# import rivers file
rivers_df <- read.csv("/Volumes/LaCie/rivers_df.csv")

# change opacity of basemap
mapatt <- attributes(map)
map_transparent <- matrix(adjustcolor(map, alpha.f = 0.4), nrow = nrow(map))
attributes(map_transparent) <- mapatt
# plot and remove the grid
mygmap <- ggmap(map_transparent) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white')) +
  geom_path(
    data = rivers_df,
    aes(long, lat, fill = NULL),
    color = "#2A788EFF",
    alpha = 0.5, size = 2)
mygmap


##########
#Map USA
##########
####### make a white and grey map of the US
usmap <- get_stamenmap( bbox = c(left = -85, bottom = 30, right = -73, top = 50), zoom = 8, maptype = "toner-background")

# change opacity of basemap
usmapatt <- attributes(usmap)
usmap_transparent <- matrix(adjustcolor(usmap, alpha.f = 0.4), nrow = nrow(usmap))
attributes(usmap_transparent) <- usmapatt
# plot and remove the grid
usgmap <- ggmap(usmap_transparent) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white'))

```

Plot Figure 1 Native range
```{R, echo=T, warning=FALSE}
##################
## Pie charts Fig. 1
##################
# import the data filerter for duplicated longitude/latitude locations
all_dedup <- read.csv("/Volumes/LaCie/all_dedup.csv")
head(all_dedup) # your data set with coordinates and infection frequncies

# subset for prior and data from this study
all_dedup$uninfected <- all_dedup$total - all_dedup$infected
past <- subset(all_dedup, publication != "Sonja 2022")
present <- subset(all_dedup, publication == "Sonja 2022")
## Europe
past$radiusW=past$total/12
mycol2 = c("#818689", "#E60001") # red and grey
##  wCer2 infection status
pie.europe <- mygmap +
  geom_point(data = past,
             aes(x = lon, y = lat),
             colour = "black", size =  0.7) +
  geom_scatterpie(aes(lon, lat, r = 0.25), data = present,
                  cols = c("uninfected", "infected"), 
                  color="black",
                  linetype = 1,
                  #size = 3,
                  alpha = 0.95, 
                  #coord_fixed(), 
                  legend_name = "legend") +
  coord_fixed() +
  geom_scatterpie_legend(seq(1, ceiling(max(past$total/12)), length = 1), 
                         x=-9, y=59, labeller = function(x) x * 12) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  scale_color_manual(values = mycol2,
                     aesthetics = c("colour", "fill"))

pie.europe
#ggsave("/Volumes/LaCie/pie.past.png", plot = pie.europe,  width = 15, height = 10, dpi = 300, units = "in", device='png')
#ggsave("/Volumes/LaCie/pie.past.svg", plot = pie.europe,  width = 15, height = 10, dpi = 300, units = "in", device='svg')
```

Plot Figure 1 introduced range
```{R, echo=T, warning=FALSE}
popUsa <- read.csv("/Volumes/LaCie/popUsa.csv")
# plot the pie chart
popUsa$radius=popUsa$total/40
mycol3 = c("#818689")
pie.usa <- usgmap +
  geom_point(data = popUsa,
             aes(x = lon, y = lat),
             colour = "black", size =  2.5) +
  geom_scatterpie(aes(lon, lat, r = radius), data = popUsa,
                  cols = c("wCer2positive", "wCer2negative"), 
                  color="black", 
                  alpha = 0.95,
                  #size = 1,
                  legend_name = "legend") +
  coord_fixed() +
  geom_scatterpie_legend(seq(1, ceiling(max(popUsa$total/40)), length = 1), 
                         x=-83.5, y=48.5, labeller = function(x) x * 40) +
  xlab("Longitude (°E)") + ylab("Latitude (°S)") +
  scale_color_manual(values = mycol3,
                     aesthetics = c("colour", "fill"))

pie.usa
```

Check Spatial and temporal variation in freq
```{R, echo=T, warning=FALSE}
# load the whole data set
all <- read.csv("/Volumes/LaCie/all.csv")
# extract france
franceall <- subset(all, state == "France")
# split france by year, clumping binomial confidence intervals for locations from the same year 
confData <- tibble()
splitData <- split(franceall, franceall$year)

for(i in 1:length(splitData)){
  confData[i,1] <- splitData[[i]]$state[1] #id
  confData[i,2] <- splitData[[i]]$location[1] #pop
  confData[i,3] <- splitData[[i]]$publication[1] # lat
  confData[i,4] <- splitData[[i]]$year[1] #lon
  confData[i,5] <- splitData[[i]]$lat[1] # year
  confData[i,6] <- splitData[[i]]$lon[1] # infection
  confData[i,7] <- splitData[[i]]$host[1] # number of infected individuals
  confData[i,8] <- splitData[[i]]$inf_status[1] # number of uninfected individuals
  confData[i,9] <- splitData[[i]]$total[1]
  confData[i,10] <- splitData[[i]]$infected[1]
  confData[i,11] <- splitData[[i]]$inf_perc[1]
  confData[i,12] <- splitData[[i]]$time[1]
  confData[i,13] <- splitData[[i]]$dist[1]
  
  #Estimate infection freq. and 95% binomial confidence intervals
  confData[i,14] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[2]
  confData[i,15] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[3]
  
} 
rm(i)
rm(splitData)

colnames(confData) <- c("state", "location", "publication", "year", 
                        "lat", "lon", "host", "inf_status", "total", "infected", "inf_perc", "time", "dist", "infFreq_lowerCI", "infFreq_upperCI")

# plot
p1 <- ggplot(confData, aes(x=factor(year, level = year), y=inf_perc)) + 
  geom_errorbar(aes(ymin=infFreq_lowerCI, 
                    ymax=infFreq_upperCI), 
                width=.2, size=1.5, color=c('grey40', "grey40", "grey40", "grey40")) +
  #geom_point(size=4) +
  geom_point(size=4, color=c('grey40', "grey40", "grey40", "#ff7700")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1),
        axis.title.x = element_text(vjust= 0.1),
        axis.text.x = element_text(angle = 0, vjust= 0.5),
        axis.text = element_text(size = 12),
        text = element_text(size = 18),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection frequency", x = "Sampling year") +
  ggtitle("France")

#### Fishers exact test
# France 2021 vs 1999
france1 <- matrix(c(confData$infected[4], confData$total[4]-confData$infected[4],
                    confData$infected[1], confData$total[1]-confData$infected[1]),
                  ncol=2)
fisher.test(france1)
# France 2021 vs 2000
france2 <- matrix(c(confData$infected[4], confData$total[4]-confData$infected[4],
                    confData$infected[2], confData$total[2]-confData$infected[2]),
                  ncol=2)
fisher.test(france2)
# France 2021 vs 20001
france3 <- matrix(c(confData$infected[4], confData$total[4]-confData$infected[4],
                    confData$infected[3], confData$total[3]-confData$infected[3]),
                  ncol=2)
fisher.test(france3)

# extract poland
polandall <- subset(all, state == "Poland")

confData <- tibble()
splitData <- split(polandall, polandall$year)

for(i in 1:length(splitData)){
  confData[i,1] <- splitData[[i]]$state[1] #id
  confData[i,2] <- splitData[[i]]$location[1] #pop
  confData[i,3] <- splitData[[i]]$publication[1] # lat
  confData[i,4] <- splitData[[i]]$year[1] #lon
  confData[i,5] <- splitData[[i]]$lat[1] # year
  confData[i,6] <- splitData[[i]]$lon[1] # infection
  confData[i,7] <- splitData[[i]]$host[1] # number of infected individuals
  confData[i,8] <- splitData[[i]]$inf_status[1] # number of uninfected individuals
  confData[i,9] <- splitData[[i]]$total[1]
  confData[i,10] <- splitData[[i]]$infected[1]
  confData[i,11] <- splitData[[i]]$inf_perc[1]
  confData[i,12] <- splitData[[i]]$time[1]
  confData[i,13] <- splitData[[i]]$dist[1]
  
  #Estimate infection freq. and 95% binomial confidence intervals
  confData[i,14] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[2]
  confData[i,15] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[3]
  
} 
rm(i)
rm(splitData)

colnames(confData) <- c("state", "location", "publication", "year", 
                        "lat", "lon", "host", "inf_status", "total", "infected", "inf_perc", "time", "dist", "infFreq_lowerCI", "infFreq_upperCI")

p2 <- ggplot(confData, aes(x=factor(year, level = year), y=inf_perc)) + 
  geom_errorbar(aes(ymin=infFreq_lowerCI, 
                    ymax=infFreq_upperCI), 
                width=.2, size=1.5, color=c('grey40', "grey40", "grey40")) +
  geom_point(size=4) +
  geom_point(size=4, color=c('grey40', "grey40", "#ff7700")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1),
        axis.title.x = element_text(vjust= 0.1),
        axis.text.x = element_text(angle = 0, vjust= 0.5),
        axis.text = element_text(size = 12),
        text = element_text(size = 18),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  geom_signif(y_position = c(1.02), xmin = c("2000"), xmax = c("2019"),
              annotation=c("*"), tip_length=0.02, color = "grey40", size =0.6, textsize=9) +
  labs(y = "", x = "Sampling year") +
  ggtitle("Poland")

# Poland 2019 vs 2000
poland1 <- matrix(c(confData$infected[3], confData$total[3]-confData$infected[3],
                    confData$infected[1], confData$total[1]-confData$infected[1]),
                  ncol=2)
fisher.test(poland1)
# Poland 2019 vs 2007
poland2 <- matrix(c(confData$infected[3], confData$total[3]-confData$infected[3],
                    confData$infected[2], confData$total[2]-confData$infected[2]),
                  ncol=2)
fisher.test(poland2)


# extract austria
austriaall <- subset(all, state == "Austria")

confData <- tibble()
splitData <- split(austriaall, austriaall$year)

for(i in 1:length(splitData)){
  confData[i,1] <- splitData[[i]]$state[1] #id
  confData[i,2] <- splitData[[i]]$location[1] #pop
  confData[i,3] <- splitData[[i]]$publication[1] # lat
  confData[i,4] <- splitData[[i]]$year[1] #lon
  confData[i,5] <- splitData[[i]]$lat[1] # year
  confData[i,6] <- splitData[[i]]$lon[1] # infection
  confData[i,7] <- splitData[[i]]$host[1] # number of infected individuals
  confData[i,8] <- splitData[[i]]$inf_status[1] # number of uninfected individuals
  confData[i,9] <- splitData[[i]]$total[1]
  confData[i,10] <- splitData[[i]]$infected[1]
  confData[i,11] <- splitData[[i]]$inf_perc[1]
  confData[i,12] <- splitData[[i]]$time[1]
  confData[i,13] <- splitData[[i]]$dist[1]
  
  #Estimate infection freq. and 95% binomial confidence intervals
  confData[i,14] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[2]
  confData[i,15] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[3]
  
} 
rm(i)
rm(splitData)

colnames(confData) <- c("state", "location", "publication", "year", 
                        "lat", "lon", "host", "inf_status", "total", "infected", "inf_perc", "time", "dist", "infFreq_lowerCI", "infFreq_upperCI")

p3 <- ggplot(confData, aes(x=factor(year, level = year), y=inf_perc)) + 
  geom_errorbar(aes(ymin=infFreq_lowerCI, 
                    ymax=infFreq_upperCI), 
                width=.2, size=1.5, color=c('grey40', "grey40", "grey40", "grey40", "grey40", "grey40")) +
  geom_point(size=4, color=c('grey40', "grey40", "grey40", "grey40", "grey40", "#ff7700")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1),
        axis.title.x = element_text(vjust= 0.1),
        axis.text.x = element_text(angle = 0, vjust= 0.5),
        axis.text = element_text(size = 12),
        text = element_text(size = 18),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "", x = "Sampling year") +
  ggtitle("Austria")

# Austria 2021 vs 1998
austria1 <- matrix(c(confData$infected[6], confData$total[6]-confData$infected[6],
                    confData$infected[1], confData$total[1]-confData$infected[1]),
                  ncol=2)
fisher.test(austria1)
# Austria 2021 vs 1999
austria2 <- matrix(c(confData$infected[6], confData$total[6]-confData$infected[6],
                     confData$infected[2], confData$total[2]-confData$infected[2]),
                   ncol=2)
fisher.test(austria2)

# Austria 2021 vs 2000
austria3 <- matrix(c(confData$infected[6], confData$total[6]-confData$infected[6],
                     confData$infected[3], confData$total[3]-confData$infected[3]),
                   ncol=2)
fisher.test(austria3)
# Austria 2021 vs 2008
austria4 <- matrix(c(confData$infected[6], confData$total[6]-confData$infected[6],
                     confData$infected[4], confData$total[4]-confData$infected[4]),
                   ncol=2)
fisher.test(austria4)
# Austria 2021 vs 2015
austria5 <- matrix(c(confData$infected[6], confData$total[6]-confData$infected[6],
                     confData$infected[5], confData$total[5]-confData$infected[5]),
                   ncol=2)
fisher.test(austria5)


# extract greece
greeceall <- subset(all, state == "Greece")

confData <- tibble()
splitData <- split(greeceall, greeceall$year)

for(i in 1:length(splitData)){
  confData[i,1] <- splitData[[i]]$state[1] #id
  confData[i,2] <- splitData[[i]]$location[1] #pop
  confData[i,3] <- splitData[[i]]$publication[1] # lat
  confData[i,4] <- splitData[[i]]$year[1] #lon
  confData[i,5] <- splitData[[i]]$lat[1] # year
  confData[i,6] <- splitData[[i]]$lon[1] # infection
  confData[i,7] <- splitData[[i]]$host[1] # number of infected individuals
  confData[i,8] <- splitData[[i]]$inf_status[1] # number of uninfected individuals
  confData[i,9] <- splitData[[i]]$total[1]
  confData[i,10] <- splitData[[i]]$infected[1]
  confData[i,11] <- splitData[[i]]$inf_perc[1]
  confData[i,12] <- splitData[[i]]$time[1]
  confData[i,13] <- splitData[[i]]$dist[1]
  
  #Estimate infection freq. and 95% binomial confidence intervals
  confData[i,14] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[2]
  confData[i,15] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[3]
  
} 
rm(i)
rm(splitData)

colnames(confData) <- c("state", "location", "publication", "year", 
                        "lat", "lon", "host", "inf_status", "total", "infected", "inf_perc", "time", "dist", "infFreq_lowerCI", "infFreq_upperCI")


p4 <- ggplot(confData, aes(x=factor(year, level = year), y=inf_perc)) + 
  geom_errorbar(aes(ymin=infFreq_lowerCI, 
                    ymax=infFreq_upperCI), 
                width=.2, size=1.5, color=c('grey40', "grey40", "grey40")) +
  geom_point(size=4, color=c('grey40', "grey40", "#ff7700")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1),
        axis.title.x = element_text(vjust= 0.1),
        axis.text.x = element_text(angle = 0, vjust= 0.5),
        axis.text = element_text(size = 12),
        text = element_text(size = 18),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "", x = "Sampling year") +
  ggtitle("Greece")

# Greece 2019 vs 1999
greece1 <- matrix(c(confData$infected[3], confData$total[3]-confData$infected[3],
                   confData$infected[1], confData$total[1]-confData$infected[1]),
                 ncol=2)
fisher.test(greece1)

# Greece 2019 vs 2000
greece2 <- matrix(c(confData$infected[3], confData$total[3]-confData$infected[3],
                    confData$infected[2], confData$total[2]-confData$infected[2]),
                  ncol=2)
fisher.test(greece2)
# combine the plots
compare <- ggarrange(p1, p3, p4, p2, nrow = 1)
compare

```

Create a spatial grid out of the data and coordinates
```{R, echo=T, warning=FALSE}
meuse <- st_as_sf(x = all_dedup,
                  coords = c("lon", "lat"),
                  crs = "+proj=longlat") #+datum=WGS84 +ellps=WGS84 +towgs84=0,0,0

bbox <- st_bbox(meuse)

meuse_grid <- meuse %>% 
  st_bbox() %>%     # determines bounding box coordinates from meuse
  st_as_sfc() %>%   # creates sfc object from bounding box
  st_make_grid(     # create grid 0.1 x 0.1 pixel size
    cellsize = c(0.1, 0.1), 
    what = "centers") %>%
  st_as_sf(crs=st_crs(meuse)) # convert to sf object
### Convert meuse samples to SpatialPointsDataFrame
meuse_sp <- as(meuse, "Spatial")

### Convert meuse grid to SpatialPixelsDataFrame, the raster/grid equivalent in the sp world
meuse_grid_sp <- as(as(meuse_grid, "Spatial"), "SpatialPixels")
```

First check for anisotropy with a directional map
```{R, echo=T, warning=FALSE}

# Create a "gstat" object 
TheGStat <- gstat(id="Infection_frequency", formula=inf_perc ~ 1, data=meuse)

TheVariogram=variogram(TheGStat, map=TRUE, cutoff=2000, width=60)
## check whether there is a north-south, etc. trend with a variogram map
plot(TheVariogram, threshold=10)
varmap <- plot(TheVariogram, threshold=10)
png(file="/Volumes/Lacie/varmap.png",
width=8, height=8, units="in", res=300)
varmap
dev.off()

```

Ok, directional map look good, no sign of directionality. Now fit a variogram model to this data. If there is no directionality, constant variance will alwas result in the presence of a sill (the curve will have a plateau).
```{R, echo=T, warning=FALSE}

#Create directional empirical variograms at 0, 45, 90, 135 degrees from north (y-axis)
meuse.aniso <- gstat::variogram(inf_perc ~ 1, meuse, cressie = T, cutoff = 2000, width = 60, alpha = c(0, 45, 90, 135))
# plot the empirical variograms
plot(meuse.aniso, xlab = "Distance (km)", ylab = "Semivariance")
# Now choose a model - I choose the exponential model here
TheModel=vgm(model='Exp' , anis=c(0, 0.5))
# Fit the model to the variogram
FittedModel <- fit.variogram(meuse.aniso, model=TheModel)
## plot results
plot(meuse.aniso, model=FittedModel, as.table=TRUE, xlab = "Distance (km)", ylab = "Semivariance")
varaniso <- plot(meuse.aniso, model=FittedModel, as.table=TRUE, xlab = "Distance (km)", ylab = "Semivariance")
#png(file="/Volumes/Lacie/varaniso.png",
#width=10, height=5, units="in", res=300)
#varaniso
#dev.off()
```
Looks good. All four directions have a distinct plateau.

Since we do not have directionality in our data, the exponential variogram model is appropriate to fit here.
```{R, echo=T, warning=FALSE}
meuse.z <- gstat::variogram(inf_perc ~ 1, cressie = T, meuse, cutoff = 2000, width = 60)
#meuse.v <- gstat::variogram(inf_perc ~ 1, cressie = T, meuse, cutoff = 6500, width = 60)
meuse.wav <- vgm(psill=0.1, model = c("Exp"), range = 100, nugget=0.05, covariance=T) 

meuse.vfit <- fit.variogram(meuse.z, meuse.wav)
meuse.vfit # give the theorethcal sill, the nugget and the range
#plot(meuse.z, meuse.vfit)

## create a fucntion to plot a fitted semivariogram with ggplot
plot_variogram <- function(v, vz, m) {
  preds = variogramLine(m, maxdist = max(v$dist))
  ggplot() + 
  theme_classic() +
    geom_point(data = vz, aes(x = dist, y = gamma, size=np)) +
    geom_line(data = preds, aes(x = dist, y = gamma)) +
    geom_vline(xintercept = 193.7, linetype="dotted", 
               color = "black", size=0.8) +
    theme(axis.title.y = element_text(vjust= 1),
          axis.title.x = element_text(vjust= 0.1),
          axis.text.x = element_text(angle = 0, vjust= 0.5),
          axis.text = element_text(size = 12),
          text = element_text(size = 18),
          axis.title = element_text(face = "bold")) +
    labs(x = "Distance (km)", y = "Semivariance", color = "sh")
}
variogram <- plot_variogram(v=meuse.z, vz=meuse.z, m=meuse.vfit)
variogram
#ggsave("/Volumes/LaCie/variogram_cut2000_w60.png", plot = variogram, width = 5, height = 5, dpi = 300, units = "in", device='png')
```

Now, using spatial interpolation of the infection frequencies at observed locations, predict infection frequency at unobserved locations using ordinary kriging. 
```{R, echo=T, warning=FALSE}
###### make kriging function to test different parameters and models #########
mykrige <- function(data, formula, model, inits1, inits2) {
  
  # create and sf object
  meuse <- st_as_sf(x = data,
                    coords = c("lon", "lat"),
                    crs = "+proj=longlat") #+datum=WGS84 +ellps=WGS84 +towgs84=0,0,0
  
  bbox <- st_bbox(meuse)
  
  # create a grid
  meuse_grid <- meuse %>% 
    st_bbox() %>%     # determines bounding box coordinates from meuse
    st_as_sfc() %>%   # creates sfc object from bounding box
    st_make_grid(     # create grid 0.1 x 0.1 pixel size
      cellsize = c(0.1, 0.1), 
      what = "centers") %>%
    st_as_sf(crs=st_crs(meuse)) # convert to sf object
  
  ### Convert meuse samples to SpatialPointsDataFrame
  meuse_sp <- as(meuse, "Spatial")
  
  ### Convert meuse grid to SpatialPixelsDataFrame, the raster/grid equivalent in the sp world
  meuse_grid_sp <- as(as(meuse_grid, "Spatial"), "SpatialPixels")
  
  meuse.v <- gstat::variogram(formula, meuse, cressie = T, cutoff = inits1[1], width = inits1[2])
  
  meuse.wav <- vgm(psill=inits2[1], model = model, range = inits2[2], nugget=inits2[3], covariance = T) 
  
  meuse.vfit <- fit.variogram(meuse.v, meuse.wav)
  meuse.vfit
  #plot(meuse.vfit)
  
  ### ordinary kriging
  prediction <- krige(formula = formula, meuse_sp, meuse_grid_sp, meuse.vfit)
  prediction
  
}
# call the function with chosent set of parameters based on the theoretical prediction
mk <- mykrige(data=all_dedup, formula=inf_perc ~ 1, model="Exp", inits1=c(1000, 50), inits2 = c(0.1, 100, 0.05))

# Coerce raster to dataframe, including coordinates
mk.df <- mk %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame;

# plot kriging predicted contour lines on the map with Danube river
krigplot <- mygmap +
  geom_contour(data=mk.df, aes(x=coords.x1, y=coords.x2, z = var1.pred, colour = ..level..), 
               binwidth = 0.025,
               size = 1.5) +
  theme(legend.title = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, 'in')) +
  scale_color_viridis_c("Prediction", option = "plasma")

krigplot
#ggsave("/Volumes/LaCie/krigplot.GLMresid_cut2000_w60.png", plot = krigplot, width = 30, height = 30, dpi = 300, units = "in", device='png')
#ggsave("/Volumes/LaCie/krigplot.MarkdownProba.png", plot = krigplot, width = 30, height = 30, dpi = 300, units = "in", device='png')
```

Mathematical modeling of the intergenrations change in infection frequency of the wCer2 in the introduced range
```{R, echo=T, warning=FALSE}

# load the us population
popUsa <- read.csv("/Volumes/LaCie/popUsa.csv")
# calculate confidence intervals
for(i in 1:nrow(popUsa)){
  
  #Estimate infection freq. and 95% binomial confidence intervals
  popUsa[i,10] <- binconf(x=as.numeric(popUsa[i,6]), n=as.numeric(popUsa[i,6]) + as.numeric(popUsa[i,7]))[2]
  popUsa[i,11] <- binconf(x=as.numeric(popUsa[i,6]), n=as.numeric(popUsa[i,6]) + as.numeric(popUsa[i,7]))[3]
}
colnames(popUsa) <- c("state", "location", "lat", "lon", 
                   "total", "inf", "uninf", "HT2", "HT1", "infFreq_lowerCI", "infFreq_upperCI")
popUsa$time<- 10
popUsa

## write a function to simulate different invasion scenarios
simInvasion <- function(relF, mu, sh, pw, time) {
  inf <- pw
  for (i in time) {
    pwt <- (relF*(1-mu)*inf[i])/(1-sh*inf[i]*(1-inf[i])-mu*sh*(inf[i])^2)
    inf <- c(inf, pwt)
    
  }
  print(inf)
  params <- data_frame(relFit = rep(relF, length(time)),
                       matfid = rep(mu, length(time)),
                       ci = rep(sh, length(time)),
                       t = time)
  tot <- cbind(params, inf)
  #print(tot)
}
sim1 <- simInvasion(relF=1, mu=0, sh=0.98, pw=0.01, time=c(0:25))
sim2 <- simInvasion(relF=1.02, mu=0, sh=0.98, pw=0.01, time=c(0:25))
sim3 <- simInvasion(relF=1.06, mu=0, sh=0.98, pw=0.01, time=c(0:25))
sim4 <- simInvasion(relF=1.1, mu=0, sh=0.98, pw=0.01, time=c(0:25))
simall <- rbind(sim1, sim2, sim3, sim4)
simall$variable <- c(rep("sim1", nrow(sim1)), rep("sim2", nrow(sim1)), rep("sim3", nrow(sim3)),
                     rep("sim4", nrow(sim4))
                     )
#simall[which(simall$inf > 1),]$inf <- 1
head(simall)

gsim <- ggplot(data=simall, aes(x=t, y=inf 
                              ,col=variable
                              )) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values=c('red','blue', 'orange', 'forestgreen'), labels = c("relF = 0.9", "relF= 1", "relF = 1.06", "relF = 1.2"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0, 25)) +
  labs(y = "Infection frequency", x = "Generations", color = "relF")
  #ggtitle(paste("sh =", sh))

gusa <- ggplot(data=popUsa, aes(x=time , y=inf)) +
  geom_point(color="grey40", size=2.5) +
  geom_errorbar(aes(ymin=infFreq_lowerCI, 
                    ymax=infFreq_upperCI), 
                width=0.8, size=1, color="grey40") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0, 25)) +
  theme_classic()

g_pw_0.01 <- gusa + geom_line(data=simall, aes(x=t, y=inf, col=variable), size=1) +
  theme_classic() +
  scale_color_manual(values=c('red','blue', 'orange', 'forestgreen'), 
  #                   labels = c("F = 1", "F= 1.02", "F = 1.06", "F = 1.1"), 
                     guide = guide_legend(reverse=TRUE)) +
  labs(y = "Infection frequency", x = "Generations", color = "relF") +
  theme(axis.title.y = element_text(vjust= 1),
        axis.title.x = element_text(vjust= 0.1),
        axis.text.x = element_text(angle = 0, vjust= 0.5),
        axis.text = element_text(size = 13),
        text = element_text(size = 18),
        #legend.text = element_text(size = 11),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position="none")
  #ggtitle(paste("sh =", sh))
g_pw_0.01


```

Mathematical modeling of the wCer2 frequency dynamics and equilibirum in the native range
```{R, echo=T, warning=FALSE}

for(i in 1:nrow(all)){
  
  #Estimate infection freq. and 95% binomial confidence intervals
  all[i,14] <- binconf(x=as.numeric(all[i,10]), n=as.numeric(all[i,10] + (all[i, 9] - all[i,10])))[2]
  all[i,15] <- binconf(x=as.numeric(all[i,10]), n=as.numeric(all[i,10] + (all[i, 9] - all[i,10])))[3]
}
colnames(all) <- c("state", "location", "publication", "year", 
                   "lat", "lon", "host", "inf_status", "total", "infected", "inf_perc", "time", "dist", "infFreq_lowerCI", "infFreq_upperCI")


all$all <- "all"
# extract bakovic
bak <- subset(all, publication == "Bakovic et al. 2018")
# extract bakovic czeck republic
bakcz<-subset(bak, state != "Hungary")
# extract bakovic hungary
bakhun <- subset(bak, state == "Hungary")

# calculate binomial confidence intervals
confData <- tibble()
splitData <- split(bakhun, bakhun$all)

for(i in 1:length(splitData)){
  confData[i,1] <- splitData[[i]]$state[1] #id
  confData[i,2] <- splitData[[i]]$location[1] #pop
  confData[i,3] <- splitData[[i]]$publication[1] # lat
  confData[i,4] <- splitData[[i]]$year[1] #lon
  confData[i,5] <- splitData[[i]]$lat[1] # year
  confData[i,6] <- splitData[[i]]$lon[1] # infection
  confData[i,7] <- splitData[[i]]$host[1] # number of infected individuals
  confData[i,8] <- splitData[[i]]$inf_status[1] # number of uninfected individuals
  confData[i,9] <- splitData[[i]]$total[1]
  confData[i,10] <- splitData[[i]]$infected[1]
  confData[i,11] <- splitData[[i]]$inf_perc[1]
  confData[i,12] <- splitData[[i]]$time[1]
  confData[i,13] <- splitData[[i]]$dist[1]
  confData[i,14] <- splitData[[i]]$all[1]
  
  #Estimate infection freq. and 95% binomial confidence intervals
  confData[i,15] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[2]
  confData[i,16] <- binconf(x=as.numeric(confData[i,10]), n=as.numeric(confData[i,10] + (confData[i, 9] - confData[i,10])))[3]
  
} 
rm(i)
rm(splitData)

colnames(confData) <- c("state", "location", "publication", "year", 
                        "lat", "lon", "host", "inf_status", "total", "infected", "inf_perc", "time", "dist", "all", "infFreq_lowerCI", "infFreq_upperCI")



relF<- seq(0.8, 1.1, by=0.01) # range of fitness effects
# set different values for maternal transmission

#mu=0.0    # perfect maternal transmission
#mu=0.01   # imperfect maternal transmission
#mu=0.03   # imperfect maternal transmission
mu=0.06    # imperfect maternal transmission

# set different levels of CI and calculate unstable and stable equilibria
sh=0.45 # CI
#phatu1 <- (1-relF)/sh # unstable equilibrium when F(1-mu)<1
phatu1 <- (sh + 1 - relF - sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) # unstable equilibrium when F(1-mu)>1
phats1 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) # stable equilibrium when F(1-mu)>1

sh=0.75 # CI
#phatu2 <- (1-relF)/sh
phatu2 <- (sh + 1 - relF - sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 
phats2 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

sh=0.98 # CI
#phatu3 <- (1-relF)/sh
phatu3 <- (sh + 1 - relF - sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 
phats3 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

yak_low_data <- as.data.frame(cbind(relF, phatu1
                                    , phatu2, phatu3, phats1, phats2, phats3))
yak_low_data <- melt(yak_low_data, id.vars="relF")
yak_low_data$lineffect <- c(rep("unstable", 3*length(relF)), rep("stable", 3*length(relF)))

# plot for different mus
mu006 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  #geom_ribbon(aes(ymin=confData$infFreq_lowerCI, ymax=confData$infFreq_upperCI), fill="grey") +
  geom_hline(yintercept=confData$inf_perc, color = "grey50", size=0.8) +
  geom_hline(yintercept=confData$infFreq_lowerCI, linetype="dotted", color = "grey50", size=0.8) +
  geom_line(data=subset(yak_low_data, lineffect=="unstable"), linetype=2, size=1) +
  geom_line(data=subset(yak_low_data, lineffect=="stable"), linetype=1, size=1) +
  theme_classic() +
  theme(axis.title.y = element_text(vjust= 1),
        axis.title.x = element_text(vjust= 0.1),
        axis.text.x = element_text(angle = 0, vjust= 0.5),
        axis.text = element_text(size = 13),
        text = element_text(size = 18),
        #legend.text = element_text(size = 11),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
        ,legend.position="none"
        ) +
  scale_color_manual(values=c('lightblue3','plum4', 'darkorange2', 'lightblue3','plum4', 'darkorange2'), labels = c("sh = 0.45", "sh = 0.75", "sh = 0.98", "sh = 0.45", "sh = 0.75", "sh = 0.98"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(y = "Equlibrium frequency", x = "Fitness", color = "sh")
  #ggtitle(paste("mu =", mu))
mu006

```
