
#########################
### load required libraries
########################
library(gstat)   # geostatistics
library(sf)      # spatial data handling
library(sp)       # spatial data handling
library(ggplot2)  # plotting
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


######################################
# Native Range analysis and ploting 
######################################


#####################
### Map of Europe
#####################

####### make a white and grey map
map <- get_stadiamap( bbox = c(left = -11, bottom = 35, right = 45, top = 61), zoom = 4, maptype = "stamen_toner_background")
ggmap(map) + 
  theme_void() + 
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
  )

# download river shape file to plot Danube on the map
rivers110 <- ne_download(scale = 110, type = 'rivers_lake_centerlines', category = 'physical')
# transform the shape file into a data frame read by ggplot
rivers_df <- fortify(rivers110)

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


#### make a green terrain map
#mymap <- get_stamenmap( bbox = c(left = -13, bottom = 35, right = 46, top = 63), zoom=6, maptype = "terrain")
#myterrain <- ggmap(map) +
#  theme_void() + 
#  theme(
#    plot.title = element_text(colour = "orange"), 
#    panel.border = element_rect(colour = "grey", fill=NA, size=2)
#  )

# plot map with the river
#myriver <- myterrain +
#  labs(x = "Longitude", y = "Latitude") +
#  geom_path(
#    data = rivers_df,
#    aes(long, lat, fill = NULL),
#    color = "#2A788EFF",
#    alpha = 0.9)


###############################
#### load the data
##############################
all <- rbind(Riegler, Schuler, Bakovic, Arthofer, USApaper)
all$inf_perc <- all$infected/all$total # infection frequency
all$time <- 2022 - all$year
all$time <- max(all$time) - all$time
#all <- subset(all, total >= 5) # filter out anything with sample size less than 5
all_dedup <- all[!duplicated(all$location),] # duplicated locations removed
write.csv(all_dedup, "/Volumes/LaCie/all_dedup.csv", row.names=FALSE)

###### make kriging function and test different parameters and models #########
mykrige <- function(data, formula, model, inits1, inits2) {
  
  meuse <- st_as_sf(x = data,
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
  #zn.idw <- idw((1-inf_perc) ~ 1, locations=meuse_sp, newdata=meuse_grid_sp, idp = 2)
  #spplot(zn.idw["var1.pred"])
  #points(meuse_sp, col="white", cex=0.5)
  
  meuse.v <- gstat::variogram(formula, meuse, cressie = T, cutoff = inits1[1], width = inits1[2])
  #plot(meuse.v, ylab=bquote(gamma), xlab=c("h (separation distance in m)"))
  
  meuse.wav <- vgm(psill=inits2[1], model = model, range = inits2[2], nugget=inits2[3], covariance = T) 
  
  meuse.vfit <- fit.variogram(meuse.v, meuse.wav)
  #plot(meuse.v, meuse.vfit)
  #plot(meuse.v, meuse.vfit, xlim = c(0, 4000), ylim = c(0, 1))
  
  ### ordinary kriging
  prediction <- krige(formula = formula, meuse_sp, meuse_grid_sp, meuse.vfit)
  prediction
  #spplot(prediction['var1.pred'], contour = T)
  
  #im <- as.image.SpatialGridDataFrame(prediction['var1.pred']) 
  #cl <- ContourLines2SLDF(contourLines(im)) 
  #spplot(cl)
}

mk <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Wav", inits1=c(2200, 100), inits2 = c(0.2, 1500, 0.05))
mk2 <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Wav", inits1=c(1800, 100), inits2 = c(0.2, 800, 0.05))
mk3 <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Sph", inits1=c(1500, 100), inits2 = c(0.2, 1000, 0.05))
mk4 <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Exp", inits1=c(2200, 100), inits2 = c(0.2, 1500, 0.05))
mk5 <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Exp", inits1=c(1500, 60), inits2 = c(0.2, 1500, 0.05)) # the one showing the cline following Danube
mk6 <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Exp", inits1=c(1000, 60), inits2 = c(0.2, 1000, 0.05)) # The best model
mk7 <- mykrige(data=all_dedup, formula=(1-inf_perc) ~ 1, model="Wav", inits1=c(1000, 60), inits2 = c(0.2, 1000, 0.05)) #  Mean square normalized error = 1.037328 NOT GOOD!!!

# Coerce raster to dataframe, including coordinates
mk.df <- mk6 %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame;

write.csv(mk.df, "/Volumes/LaCie/krig.pred.csv", row.names = F)

subset(mk.df, var1.pred > 0.9)

# plot kriging predicted contour lines on the map & Danube river
krigplot <- mygmap +
  geom_contour(data=mk.df, aes(x=coords.x1, y=coords.x2, z = var1.pred, colour = ..level..), binwidth = 0.025, size = 1.5) +
  #theme_void() +
  theme(legend.title = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, 'in')) +
  scale_color_viridis_c("Prediction", option = "plasma")
ggsave("/Volumes/LaCie/krigplot.rus.png", plot = krigplot, width = 30, height = 30, dpi = 300, units = "in", device='png')


### model cross-validation using gstat.cv function
meuse.g <- gstat(formula = (1-inf_perc) ~ 1, data = meuse)
meuse.gtv <- gstat(meuse.g, model = vgm(0.2, "Exp", 1000, 0.05), fill.all = TRUE)
meuse.gev <- variogram(meuse.g, cutoff = 1000, width = 60)
meuse.fit = fit.lmc(meuse.gev, meuse.gtv)
out = gstat.cv(meuse.fit, nmax=60, nfold = nrow(meuse), maxdist = 1000) 
summary(out)
mean(out$residual, na.rm = T) # mean error, ideally 0:
mean(out$residual^2, na.rm = T) # MSPE, ideally small
mean(out$zscore^2, na.rm = T) # Mean square normalized error, ideally close to 1
cor(out$observed, out$observed - out$residual, use = "complete.obs") # correlation observed and predicted, ideally 1
cor(out$observed - out$residual, out$residual, use = "complete.obs") # correlation predicted and residual, ideally 0


##################
## Pie charts
##################

## Europe
pops$radiusW=pops$total/12
mycol2 = c("#E60001", "#818689") # red and grey
##  wCer2 infection status
pie.europe <- mygmap +
  #geom_point(data = all_dedup,
   #          aes(x = lon, y = lat),
   #          colour = "grey20", size =  2.3) +
  geom_scatterpie(aes(lon, lat, r = radiusW), data = pops,
                  cols = c("wCer2positive", "wCer2negative"), 
                  color="black",
                  linetype = 1,
                  #size = 3,
                  alpha = 0.95, 
                  #coord_fixed(), 
                  legend_name = "legend") +
  coord_fixed() +
  geom_scatterpie_legend(seq(1, ceiling(max(pops$total/12)), length = 2), 
                         x=-9, y=59, labeller = function(x) x * 12) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  scale_color_manual(values = mycol2,
                     aesthetics = c("colour", "fill"))

#ggsave("/Volumes/LaCie/pie.europe.png", plot = pie.europe,  width = 30, height = 30, dpi = 300, units = "in", device='png')
ggsave("/Volumes/LaCie/pie.europe.rus.png", plot = pie.europe,  width = 10, height = 10, dpi = 300, units = "in", device='png')
ggsave("/Volumes/LaCie/pie.europe.rus.svg", plot = pie.europe,  width = 15, height = 10, dpi = 300, units = "in", device='svg')


###########
# USA
###########

####### make a white and grey map of the US
usmap <- get_stadiamap( bbox = c(left = -85, bottom = 30, right = -73, top = 50), zoom = 8, maptype = "stamen_toner_background")
ggmap(usmap) + 
  theme_void() + 
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
)

# change opacity of basemap
usmapatt <- attributes(usmap)
usmap_transparent <- matrix(adjustcolor(usmap, alpha.f = 0.4), nrow = nrow(usmap))
attributes(usmap_transparent) <- usmapatt
# plot and remove the grid
usgmap <- ggmap(usmap_transparent) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white'))

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

ggsave("/Volumes/LaCie/pie.usa.png", plot = pie.usa,  width = 5, height = 10, dpi = 300, units = "in", device='png')
ggsave("/Volumes/LaCie/pie.usa.svg", plot = pie.usa,  width = 5, height = 10, dpi = 300, units = "in", device='svg')


####################
# Supplements
####################

# plot Fig. S1 krigign prediction variance as contour lines and overlay observed locations
krigplot.var <- mygmap +
  geom_contour(data=mk.df, aes(x=coords.x1, y=coords.x2, z = var1.var, colour = ..level..), binwidth = 0.005, size = 2) +
  #theme_void() +
  geom_point(data=all_dedup, aes(x=lon, y=lat), color = "grey20", size = 3) +
  theme(legend.title = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, 'in')) +
  scale_color_viridis_c("Prediction variances", option = "viridis")
ggsave("/Volumes/LaCie/krigplot.var.rus.png", plot = krigplot.var, width = 30, height = 30, dpi = 300, units = "in", device='png')

# plot kriging predicted contour lines on the map & Danube river
krigplot <- mygmap +
  geom_contour(data=mk.df, aes(x=coords.x1, y=coords.x2, z = var1.pred, colour = ..level..), binwidth = 0.025, size = 1.5) +
  geom_point(data=all_dedup, aes(x=lon, y=lat), color = "grey20", size = 6) +
  theme(legend.title = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, 'in')) +
  scale_color_viridis_c("Probability", option = "plasma")
ggsave("/Volumes/LaCie/krigplot.rus.locs.png", plot = krigplot, width = 30, height = 30, dpi = 300, units = "in", device='png')

# plot Fig. S2 observed locations only
pointplot <- mygmap +
  geom_point(data=all_dedup, aes(x=lon, y=lat), color = "grey20", size = 6)
ggsave("/Volumes/LaCie/pointplot2.png", plot = pointplot, width = 30, height = 30, dpi = 300, units = "in", device='png')
ggsave("/Volumes/LaCie/pointplot2.svg", plot = pointplot, width = 30, height = 30, dpi = 300, units = "in", device='svg')


##################
## Pie charts as per Hannes's request to have pies & past location in one figure
##################
all_dedup$uninfected <- all_dedup$total - all_dedup$infected
past <- subset(all_dedup, publication != "Sonja 2022")
present <- subset(all_dedup, publication == "Sonja 2022")
## Europe
past$radiusW=past$total/42
mycol2 = c("#818689", "#E60001") # red and grey
##  wCer2 infection status
pie.europe <- mygmap +
  geom_point(data = past,
             aes(x = lon, y = lat),
             colour = "black", size =  0.7) +
  geom_scatterpie(aes(lon, lat, r = 0.25), data = past,
                  cols = c("uninfected", "infected"), 
                  color="black",
                  linetype = 1,
                  #size = 3,
                  alpha = 0.95, 
                  #coord_fixed(), 
                  legend_name = "legend") +
  coord_fixed() +
  #geom_scatterpie_legend(seq(1, ceiling(max(past$total/42)), length = 4), 
  #                       x=-9, y=59, labeller = function(x) x * 42) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  scale_color_manual(values = mycol2,
                     aesthetics = c("colour", "fill"))

#ggsave("/Volumes/LaCie/pie.europe.png", plot = pie.europe,  width = 30, height = 30, dpi = 300, units = "in", device='png')
ggsave("/Volumes/LaCie/pie.past.png", plot = pie.europe,  width = 15, height = 10, dpi = 300, units = "in", device='png')
ggsave("/Volumes/LaCie/pie.past.svg", plot = pie.europe,  width = 15, height = 10, dpi = 300, units = "in", device='svg')
#


