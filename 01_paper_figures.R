###
### Reproduction of figures for BlockR paper
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
library(vec2dtransf) #for transformations
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

#download census tracts for Middlesex
spoly_tracts <- tracts(state = "MA", county = "Middlesex") %>% as("Spatial")

#clust units before scaling, for plotting
clust_tracts_plot <- clust_units(spoly_tracts, num_clust = 10)

#scale and shift polygon
bbox(spoly_tracts)
plot(spoly_tracts)
scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
bbox(scaled_tracts)
scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
bbox(scaled_tracts)
plot(scaled_tracts)

#cluster the tracts
clust_tracts <- clust_units(scaled_tracts, num_clust = 10)

#define center of the current shape
shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% 
  as("Spatial") %>% as.data.frame()

#center the shape at the origin
clust_tracts <- elide.polygonsdf(clust_tracts, 
                                 shift = c(-shp_centroid$coords.x1[1],
                                           -shp_centroid$coords.x2[1]))

#plot clustered shapes
plot(clust_tracts)

#rotate the tracts
rotate_angle <-110 #fix this value for now
rot_tracts <- elide.polygonsdf(clust_tracts, rotate = rotate_angle, 
                               center = c(0,0))
bbox(rot_tracts)

#plot to see result
plot(rot_tracts)
plot(clust_tracts, add = T) #see them both overlayed

#blockify tracts
block_tracts <- blockR(rot_tracts, num_x_blocks = 13)
block_tracts_orig <- block_tracts

#round all blocks
block_tracts <- round_poly_coords(block_tracts, digits = 5)

#reassign projection
proj4string(block_tracts) <- proj4string(spoly_tracts)

#add border change process
#this process takes time
block_tracts <- border_change(block_tracts, showoutput = T, quiet = F)

#plot to see result
plot(block_tracts)

#rotate the first shape 0 degrees
spoly_tracts_rot <- elide.polygonsdf(scaled_tracts, rotate = 0, center = c(0,0))

#assign projections
proj4string(spoly_tracts_rot) <- proj4string(spoly_tracts)
proj4string(block_tracts) <- proj4string(spoly_tracts)

#make outlines of both shapes
outlines <- make_outlines(spoly_tracts_rot, block_tracts)

#store in separate objects for each shape
shp1_coords <- outlines[[1]]
shp2_coords <- outlines[[2]]

#plot the coordinates from both shapes
plot(shp1_coords, type = "l")
plot(shp2_coords, type = "l")

#distance between two shapes
hausdorff_dist(shp1_coords, shp2_coords)

#format shapes for plotting
block_sf_orig <- st_as_sf(block_tracts_orig)

block_tracts_plot <- ggplot() + 
  geom_sf(aes(geometry = geometry, fill = GEOID), 
          data = block_sf_orig, col = "grey",
          inherit.aes = F, linewidth =0.7) + 
  theme_void() + theme(legend.position = "none") 


#show image
block_tracts_plot

#produce plot for paper: Figure 4c
cairo_pdf(file = paste(getwd(), "/output/images/figure4c", sep =""),
          width = 4.5,
          height = 5)
block_tracts_plot
dev.off()

#for mapping on original scale
mid_map <- ggmap(get_stadiamap(c(left = bbox(clust_tracts_plot)[1,1],
                                 bottom = bbox(clust_tracts_plot)[2,1], 
                                  right = bbox(clust_tracts_plot)[1,2] ,
                                 top = bbox(clust_tracts_plot)[2,2]), zoom = 10)) + 
  theme_void()

mid_map 

#format shapes for plotting
clust_block <- st_as_sf(clust_tracts_plot)

mid_clust <- mid_map + 
  geom_sf(aes(geometry = geometry, fill = GEOID), 
          data = clust_block, col = "grey",
          inherit.aes = F, linewidth =0.7) + 
  theme_void() + theme(legend.position = "none") 

mid_clust

#produce figure for paper: Figure 4b
cairo_pdf(file = paste(getwd(), "/output/images/figure4b.pdf", sep =""),
          width = 7, 
          height = 5) 
mid_clust
dev.off()

#format shapes for plotting
mid_orig_shp <- st_as_sf(spoly_tracts)

mid_orig <- mid_map + 
  geom_sf(aes(geometry = geometry), 
          data = mid_orig_shp, col = "grey",  fill = "darkred",
          inherit.aes = F, linewidth =0.3, alpha = 0.4) + 
  theme_void() + theme(legend.position = "none") 

mid_orig

#produce plot for paper: Figure 4a
cairo_pdf(file = paste(getwd(), "/output/images/figure4a.pdf", sep =""),
          width = 7, 
          height = 5) 
mid_orig
dev.off()


###
### Figure 6-7 series: on Middlesex
###

#one example of a blockifed shape
block_sf2 <- st_as_sf(block_tracts)

p2_no_title <- ggplot() + 
  geom_sf(aes(geometry = geometry, fill = GEOID), 
          data = block_sf2, col = "grey",
          inherit.aes = F, linewidth =0.7) + 
  theme_void() + theme(legend.position = "none") 
p2_no_title

#before border change process
block_sf2 <- st_as_sf(block_tracts_orig)
p2_orig <- ggplot() + 
  geom_sf(aes(geometry = geometry, fill = GEOID), 
          data = block_sf2, col = "grey",
          inherit.aes = F, linewidth =0.7) + 
  theme_void() + theme(legend.position = "none")

p2_orig

###
### Figure 6 subfigures
###
#plot for paper
cairo_pdf(file = paste(getwd(), "/output/images/figure6a.pdf", sep =""),
          width = 5,
          height = 5)
p2_no_title #seed 6, angle = 110 (same image)
dev.off()

cairo_pdf(file = paste(getwd(), "/output/images/figure6b.pdf", sep =""),
          width = 5,
          height = 5)
p2_no_title #seed = 11, angle = 110 (same image)
dev.off()

###
### Figure 7 subfigures
###
cairo_pdf(file = paste(getwd(), "/output/images/figure7a.pdf", sep =""),
          width = 5,
          height = 5)
p2_orig #seed = 7, angle = 50
dev.off()

cairo_pdf(file = paste(getwd(), "/output/images/figure7b.pdf", sep =""),
          width = 5,
          height = 5)
p2_no_title #seed = 7, angle = 50
dev.off()


cairo_pdf(file = paste(getwd(), "/output/images/figure7c.pdf", sep =""),
          width = 5,
          height = 5)
p2_no_title #seed = 11, angle = 50
dev.off()


##############################################################################
### Maps for Durham
##############################################################################

#start clean
rm(list = ls())

#load blockify functions
source("00_bg_functions.R")
source("00_rgeos_functions.R")
source("00_maptools_functions.R")
source("00_surveillance_functions.R")
source("00_border_change.R")

#download census tracts for Durham
spoly_tracts <- tracts(state = "NC", county = "Durham") %>% as("Spatial")

#clust units before scaling, for plotting
clust_tracts_plot <- clust_units(spoly_tracts, num_clust = 10)

#scale and shift polygons
bbox(spoly_tracts)
plot(spoly_tracts)
scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
bbox(scaled_tracts)
scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
bbox(scaled_tracts)
plot(scaled_tracts)

#cluster the tracts
clust_tracts <- clust_units(scaled_tracts, num_clust = 10)

#plot to see result
plot(clust_tracts)

#rotate the tracts
rotate_angle <-50 #fix this value, for now

#define center of the current shape
shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% 
  as("Spatial") %>% as.data.frame()

#center the shape at the origin
clust_tracts <- elide.polygonsdf(clust_tracts, 
                                 shift = c(-shp_centroid$coords.x1[1],
                                           -shp_centroid$coords.x2[1]))

#rotate shapes
rot_tracts <- elide.polygonsdf(clust_tracts, rotate = rotate_angle,
                               center = c(0,0))
bbox(rot_tracts)

#plot to see result
plot(rot_tracts)
plot(clust_tracts, add = T) #if you want to see both, overlayed

#blockify tracts
block_tracts <- blockR(rot_tracts, num_x_blocks = 10)

#round all blocks
block_tracts <- round_poly_coords(block_tracts, digits = 5)

#reassign projection
proj4string(block_tracts) <- proj4string(spoly_tracts)

#add border change process
block_tracts <- border_change(block_tracts, showoutput = T, quiet = F)

#plot to see result
plot(block_tracts)

#rotate the first shape 0 degrees
spoly_tracts_rot <- elide.polygonsdf(scaled_tracts, rotate = 0, center = c(0,0))

proj4string(spoly_tracts_rot) <- proj4string(spoly_tracts)
proj4string(block_tracts) <- proj4string(spoly_tracts)

#make outlines of both shapes
outlines <- make_outlines(spoly_tracts_rot, block_tracts)

#store in separate objects for each shape
shp1_coords <- outlines[[1]]
shp2_coords <- outlines[[2]]

#plot the coordinates from both shapes
plot(shp1_coords, type = "l")
plot(shp2_coords, type = "l")

#distance between two shapes
hausdorff_dist(shp1_coords, shp2_coords)

#for mapping on original scale
durh_map <- ggmap(get_stadiamap(c(left = bbox(clust_tracts_plot)[1,1]-0.03, 
                                  bottom = bbox(clust_tracts_plot)[2,1]-0.02, 
                                  right = bbox(clust_tracts_plot)[1,2]+0.01 ,
                                  top = bbox(clust_tracts_plot)[2,2])+0.01, 
                                zoom = 10)) + 
  theme_void()

# format shapes for plotting
block_sf_orig <- st_as_sf(spoly_tracts)

durh_orig <- durh_map + 
  geom_sf(aes(geometry = geometry),
          alpha=0.25, color = "grey", fill = "black",
          data = block_sf_orig,
          inherit.aes = F, linewidth =0.5) + 
  theme_void() + theme(legend.position = "none") 

durh_orig

#produce figure for paper: Figure 1a
cairo_pdf(file = paste(getwd(), "/output/images/figure1a.pdf", sep =""),
          width = 7, 
          height = 5) 
durh_orig
dev.off()

# format shapes for plotting
block_tracts_sf <- st_as_sf(block_tracts)

block_tracts_plot <- ggplot() + 
  geom_sf(aes(geometry = geometry, fill = GEOID),
          color = "grey", 
          data = block_tracts_sf,
          inherit.aes = F, linewidth =0.5) + 
  theme_void() + theme(legend.position = "none") 

#produce figure for paper: Figure 1b
cairo_pdf(file = paste(getwd(), "/output/images/figure1b.pdf", sep =""),
          width = 5,
          height = 5)
block_tracts_plot
dev.off()
