###
### Illustration of BlockR Method on user-specified county/state
###

#start clean
rm(list = ls())

#load libraries
library(sp)
library(sf)
library(ggplot2)
library(broom)
library(ggmap)
library(raster)
library(tidyverse)
library(gridExtra)
library(tigris)
library(vec2dtransf) 
library(rgeoda)
library(terra)
library(pracma)
library(spdep)
library(igraph)
library(smoothr)

#load blockR functions
source("00_bg_functions.R")
source("00_rgeos_functions.R")
source("00_maptools_functions.R")
source("00_surveillance_functions.R")
source("00_border_change.R")

###
### Input of parameters by user
###
state_run <- "MA" #type the chosen state here
county_run <- "Middlesex" #type the chosen county here
num_x_blocks_run <- 10 #specify number of blocks in the the x-direction
angle <- 100 #specify rotation angle in degrees
num_clust <- 9 #set the number of clusters to used for the areal units

###
### Execution of BlockR methods
###

#download census tracts for county of interest
#can vary to a different type of areal unit, if desired
spoly_tracts <- tracts(state = state_run, county = county_run) %>% as("Spatial")

#scale and shift polygon
scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
bbox(scaled_tracts)
scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
bbox(scaled_tracts) #verify the bounding box is within (-1,-1)x(1,1)

#cluster the tracts according to provided number
clust_tracts <- clust_units(scaled_tracts, num_clust = num_clust)
plot(clust_tracts)

#define center of the current shape
shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% 
  as("Spatial") %>% as.data.frame()

#center the shape at the origin
clust_tracts <- elide.polygonsdf(clust_tracts, 
                                 shift = c(-shp_centroid$coords.x1[1],
                                           -shp_centroid$coords.x2[1]))

#rotate the tracts
rot_tracts <- elide.polygonsdf(clust_tracts, rotate = angle, center = c(0,0))
proj4string(rot_tracts) <- proj4string(spoly_tracts)

#blockify tracts
block_tracts <- blockR(rot_tracts, num_x_blocks = num_x_blocks_run) 

#round all blocks
block_tracts <- round_poly_coords(block_tracts, digits = 5)

#reassign projection
proj4string(block_tracts) <- proj4string(rot_tracts)

#add border change process
#this border change process takes time
block_tracts <- border_change(block_tracts, showoutput = T, quiet = F) 

#scale and shift polygon
block_tracts <- elide.polygonsdf(block_tracts, scale = 2)
block_tracts <- elide.polygonsdf(block_tracts, shift = c(-1,-1))

#define center of the current shape
block_centroid <- st_centroid(st_union(st_as_sf(block_tracts))) %>% 
  as("Spatial") %>% as.data.frame()

#center the shape at the origin
block_tracts <- elide.polygonsdf(block_tracts, 
                                 shift = c(-block_centroid$coords.x1[1],
                                           -block_centroid$coords.x2[1]))

#project the shapes
proj4string(block_tracts) <- proj4string(spoly_tracts)

#randomize GEOID for color for public release
num_geoid <- length(unique(block_tracts$GEOID))
cols = rainbow(num_geoid)

#new GEOIDs, for public release
public_geoids <- paste("GEO", 1:num_geoid, sep = "") %>% sample()
names(public_geoids) <- unique(block_tracts@data$GEOID)

#recode the GEOIDs
block_tracts@data <- block_tracts@data %>% 
  mutate(public_GEOID = recode(GEOID, !!!public_geoids))

#deleting original GEOID for public release
block_tracts@data <- block_tracts@data %>% select(public_GEOID)

#convert format of blocks for 
block_sf <- st_as_sf(block_tracts)

#image of blockified shape with randomized colors
block_plot <- ggplot() + 
  geom_sf(aes(geometry = geometry, fill = public_GEOID), 
          data = block_sf, col = "grey",
          inherit.aes = F, linewidth =0.7) + 
  theme_void() + theme(legend.position = "none") + 
  scale_fill_manual(values=cols)

#clean up environment
rm(list=setdiff(ls(), c("block_plot", "block_sf", "block_tracts")))

block_plot

###
### The block_sf or block_tracts object can be used as the anonymized shape
###