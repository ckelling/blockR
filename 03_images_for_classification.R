###
### BlockR code - Code to produce counties for image classification
###

set.seed(3)

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
library(surveillance) #to see if a polygon is on the border
library(igraph)
library(smoothr)

#load blockR functions
source("00_bg_functions.R")
source("00_rgeos_functions.R")
source("00_maptools_functions.R")
source("00_surveillance_functions.R")
source("00_border_change.R")

library(usdata)
data("county_complete")
county_data <- county_complete[, c("name", "state")]

#random 50 counties
r_counties <- sample(1:nrow(county_data), size = 50, replace = FALSE)

county_data <- county_data[r_counties,]
county_data$name <- str_replace(county_data$name, " County", "")

###
###  Not random color
###
for(j in 1:nrow(county_data)){ 
  #overall parameters
  state_run <- county_data[j,2]
  county_run <- county_data[j,1]
  num_x_blocks_run <- 10
  num_angles <- 100 

  #download census tracts for county of interest
  spoly_tracts <- tracts(state = state_run, county = county_run) %>% as("Spatial")
  
  #identify angles that will be used
  rotate_deg <- runif(n = num_angles, 0,360)
  
  #scale and shift polygon
  scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
  bbox(scaled_tracts)
  scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
  bbox(scaled_tracts)
  
  #set number of clusters as being minimum of 9 or number of tracts
  num_clust <- min(9, length(spoly_tracts))
  clust_tracts <- clust_units(scaled_tracts, num_clust = num_clust)
  shape_list <- list(length=length(num_angles))
  
  plot(clust_tracts)
  
  for(i in 1:num_angles){
    if(length(clust_tracts) < 7){break}
    
    #if not starting over, load the old shape list
    if(i != 1){
      load(file = paste("data/image_classification/", county_run,"_",state_run, "_shapes.Rdata", sep = ""))
    }
    
    #extract rotation angle
    angle_i <- rotate_deg[i]
    
    #define center of the current shape
    shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% as("Spatial") %>% as.data.frame()
    #center the shape at the origin
    clust_tracts <- elide.polygonsdf(clust_tracts, 
                                     shift = c(-shp_centroid$coords.x1[1],
                                               -shp_centroid$coords.x2[1]))
    
    #rotate the tracts
    rot_tracts <- elide.polygonsdf(clust_tracts, rotate = angle_i, center = c(0,0))
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
    block_centroid <- st_centroid(st_union(st_as_sf(block_tracts))) %>% as("Spatial") %>% as.data.frame()
   
     #center the shape at the origin
    block_tracts <- elide.polygonsdf(block_tracts, 
                                     shift = c(-block_centroid$coords.x1[1],
                                               -block_centroid$coords.x2[1]))
    
    #rotate the first shape 0 degrees
    spoly_tracts_rot <- elide.polygonsdf(clust_tracts, rotate = 0, center = c(0,0))
    proj4string(spoly_tracts_rot) <- proj4string(spoly_tracts)
    proj4string(block_tracts) <- proj4string(spoly_tracts)
  
    
    #save shape in list
    shape_list[[i]] <- block_tracts
    save(shape_list, file = paste("data/image_classification/", county_run,"_",state_run, "_shapes.Rdata", sep = ""))
    print(paste(county_run, "*******"))
    
    #format blocked tracts as sf objects
    block_sf <- st_as_sf(block_tracts)
    
    #not random color
    num_geoid <- length(unique(block_sf$GEOID))
    cols = rainbow(num_geoid)
    
    block_plot <- ggplot() + 
      geom_sf(aes(geometry = geometry, fill = GEOID), 
              data = block_sf, col = "grey",
              inherit.aes = F, linewidth =0.7) + 
      theme_void() + theme(legend.position = "none") + 
      scale_fill_manual(values=cols)
    
    #save images with non-random images
    ggsave(block_plot, filename = paste("data/image_classification/not_random/", 
                                        county_run,state_run,"_", i, ".png", 
                                        sep = ""), 
           device = "png", width = 2, height = 2)
    
  }
}


####
#### remake pictures with random color
####
rm(list = ls())
library(tidyverse)
library(sf)

#check the file names (counties that are done)
done_counties <- list.files(path = "data/image_classification/")

#format completed counties
done_counties <- done_counties %>% str_replace("_[:digit:]*.png", "") %>%
  unique() %>%
  str_subset(".Rdata")


for(j in 1:length(done_counties)){
  #load this county
  file_j <- done_counties[j]
  
  load(file= paste("data/image_classification/", file_j, sep = ""))
  
  for(i in 1:length(shape_list)){
    #format blocked tracts as sf objects
    block_sf <- st_as_sf(shape_list[[i]])
    county_run <- str_split(file_j, pattern = "_")[[1]][1]
    state_run <- str_split(file_j, pattern = "_")[[1]][2]
    
    #randomize color
    num_geoid <- length(unique(block_sf$GEOID))
    cols = rainbow(num_geoid)[sample(1:num_geoid,num_geoid)]
    
    #create plot
    block_plot <- ggplot() + 
      geom_sf(aes(geometry = geometry, fill = GEOID), 
              data = block_sf, col = "grey",
              inherit.aes = F, linewidth =0.7) + 
      theme_void() + theme(legend.position = "none") + 
      scale_fill_manual(values=cols)
    
    ggsave(block_plot, filename = paste("data/image_classification/random/", 
                                        county_run,state_run,"_", i, ".png", sep = ""), 
           device = "png", width = 2, height = 2)
  }
  
}


