###
### BlockR Method
###  Code to reproduce paper analysis that involves more than one rotation
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
library(surveillance)
library(igraph)
library(smoothr)
library(reshape2)


#load blockR functions
source("00_bg_functions.R")
source("00_rgeos_functions.R")
source("00_maptools_functions.R")
source("00_surveillance_functions.R")
source("00_border_change.R")


###
### Figures 9-11, part b: rotate original shape in 5 degree increments and 
###                     calculate H distance
###
rm(list = ls())

#overall parameters
state_run <- "NC"
county_run <- "Durham"
num_x_blocks_run <- 10
num_angles <- 73 #set as 73 for five degree increments

#number of x blocks
#durham, NC: 10
#middlesex, MA: 13
#warren, OH: 10

#load blockR functions
source("00_bg_functions.R")
source("00_rgeos_functions.R")
source("00_maptools_functions.R")
source("00_surveillance_functions.R")
source("00_border_change.R")

#download census tracts for chosen state and county
spoly_tracts <- tracts(state = state_run, county = county_run) %>% as("Spatial")

h_dist <- rep(NA, num_angles)
rotate_deg <- seq(0,360, length.out=num_angles)


#scale and shift polygon
scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
bbox(scaled_tracts)
scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
bbox(scaled_tracts)


#make blank plot to see rotation of the shape
plot(1, type = "n", xlab = "",
     ylab = "", xlim = c(-2,2),
     ylim = c(-2,2))
for(i in 1:length(h_dist)){
  #extract rotation angle
  angle_i <- rotate_deg[i]
  
  #define center of the current shape
  shp_centroid <- st_centroid(st_union(st_as_sf(scaled_tracts))) %>% 
    as("Spatial") %>% as.data.frame()
  
  #center the shape at the origin
  scaled_tracts <- elide.polygonsdf(scaled_tracts, 
                                   shift = c(-shp_centroid$coords.x1[1],
                                             -shp_centroid$coords.x2[1]))
  
  #rotate around the origin
  rot_tracts <- elide.polygonsdf(scaled_tracts, rotate = angle_i, 
                                 center = c(0,0))
  proj4string(rot_tracts) <- proj4string(spoly_tracts)
  
  #rotate the first shape 0 degrees
  spoly_tracts_rot <- elide.polygonsdf(scaled_tracts, rotate = 0, 
                                       center = c(0,0))
  proj4string(spoly_tracts_rot) <- proj4string(spoly_tracts)

  
  #make outlines of both shapes
  outline1 <- outline_one_shp(spoly_tracts_rot)
  outline2 <- outline_one_shp(rot_tracts)
  
  shp1_coords <- outline1 %>% as.matrix()
  shp2_coords <- outline2 %>% as.matrix()
  
  #plot the coordinates from both shapes
  points(shp1_coords, type = "l")
  points(shp2_coords, type = "l")
  
  #distance between two shapes
  h_dist[i] <- hausdorff_dist(shp1_coords, shp2_coords)
  
}

#combine data for use in later plot
df_plot1 <- cbind(rotate_deg, h_dist) %>% as.data.frame()

###
### Completing Figures 9-11, part b: Rotated distance between shape and rotated 
###                     and blockified shape
###

#download census tracts for chosen state and county
spoly_tracts <- tracts(state = state_run, county = county_run) %>% as("Spatial")

#for storage of distances, create angles for rotation
h_dist <- rep(NA, num_angles)
rotate_deg <- seq(0,360, length.out=num_angles)

#scale and shift polygon
scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
bbox(scaled_tracts)
scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
bbox(scaled_tracts)

#10 clusters for all state/counties
clust_tracts <- clust_units(scaled_tracts, num_clust = 10)

#make blank plot to see rotation of the shape
plot(1, type = "n", xlab = "",
     ylab = "", xlim = c(-2,2),
     ylim = c(-2,2))
for(i in 1:length(h_dist)){
  #extract rotation angle
  angle_i <- rotate_deg[i]
  
  #define center of the current shape
  shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% 
    as("Spatial") %>% as.data.frame()
  #center the shape at the origin
  clust_tracts <- elide.polygonsdf(clust_tracts, 
                                   shift = c(-shp_centroid$coords.x1[1],
                                             -shp_centroid$coords.x2[1]))
  
  #rotate the tracts
  rot_tracts <- elide.polygonsdf(clust_tracts, rotate = angle_i, 
                                 center = c(0,0))
  proj4string(rot_tracts) <- proj4string(spoly_tracts)
  
  #blockify tracts
  block_tracts <- blockR(rot_tracts, num_x_blocks = num_x_blocks_run) 
  
  #round all blocks
  block_tracts <- round_poly_coords(block_tracts, digits = 5)
  
  #reassign projection
  proj4string(block_tracts) <- proj4string(rot_tracts)
  
  #add border change process, this process takes some time
  block_tracts <- border_change(block_tracts, showoutput = F, quiet = F)
  
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

  #rotate the first shape 0 degrees
  spoly_tracts_rot <- elide.polygonsdf(clust_tracts, rotate = 0, 
                                       center = c(0,0))
  proj4string(spoly_tracts_rot) <- proj4string(spoly_tracts)
  proj4string(block_tracts) <- proj4string(spoly_tracts)
  
  #make outlines of both shapes
  outline1 <- outline_one_shp(spoly_tracts_rot)
  outline2 <- outline_one_shp(block_tracts)
  
  #format outline
  shp1_coords <- outline1 %>% as.matrix()

  
  if(class(outline2) == "SpatialPolygons"){
    outline2 <- round_poly_coords_sub(outline2, digits =4)
    outline2 <- unionSpatialPolygons_pck(outline2, IDs = rep(1, length(outline2)))
    
    #if there are still more than one polygon defining the border, round to 3 digits
    if(length(outline2@polygons[[1]]@Polygons) > 1){
      outline2 <- round_poly_coords_sub(outline2, digits =3)
      outline2 <- unionSpatialPolygons_pck(outline2, IDs = rep(1, length(outline2)))
      shp2_coords <- outline2@polygons[[1]]@Polygons[[1]]@coords
    }else{
      #otherwise accept the outline
      shp2_coords <- outline2@polygons[[1]]@Polygons[[1]]@coords
    }
  }else{
    shp2_coords <- outline2 %>% as.matrix()
  }
  
  #distance between two shapes
  h_dist[i] <- hausdorff_dist(shp1_coords, shp2_coords)
  
  #plot the coordinates from both shapes
  points(shp1_coords, type = "l")
  points(shp2_coords, type = "l")

  #print progress
  print(i)
  
}

#combine data
df_plot2 <- cbind(rotate_deg, h_dist) %>% as.data.frame()

#combine the plots
df_plot2$shape <- "Blockified Shape"
df_plot1$shape <- "Original Shape"

#combine and plot
full_df <- bind_rows(df_plot1, df_plot2)
plot_b <- ggplot(data=full_df, aes(x=rotate_deg, y=h_dist, col = shape))+ geom_point()+geom_line() +
  labs(x = "Rotation Degree", y = "Hausdorff Distance", 
       title = paste("Distance between Original Shape and Rotated Shape,\n",
                     county_run, " County", sep = ""),
       col = "Rotated Shape\nType")+theme_light() +
  theme(axis.text=element_text(size=13, angle=0, vjust=0.3),
        axis.title =element_text(size=16),
        axis.text.y=element_text(size=16),
        legend.text=element_text(size=11),
        legend.title = element_text(size=13),
        plot.title=element_text(size=17))
plot_b


#store data for this specific county for use for selberg bound
save(full_df, county_run, file = paste(getwd(), "/data/plot_data/hdist_plotb_", 
                                       county_run, ".Rdata", sep = ""))


#create part b of figures 9-11 for current county
cairo_pdf(file = paste(getwd(), "/output/",county_run,"_figureb.pdf",
                       sep =""),
          width = 7, 
          height = 5) 
plot_b
dev.off()


###
### Figures 9-11, part c: Heatmaps of many combinations of rotation angles 
###

#download census tracts for chosen state and county
spoly_tracts <- tracts(state = state_run, county = county_run) %>% as("Spatial")

#number of angles to test (fewer here)
num_angles <- 18

#create storage and rotation objects
h_dist_mat <- matrix(nrow=num_angles,ncol=num_angles)
rotate_deg <- seq(0,360, length.out=num_angles)

#scale and shift polygon
scaled_tracts <- elide.polygonsdf(spoly_tracts, scale = 2)
bbox(scaled_tracts)
scaled_tracts <- elide.polygonsdf(scaled_tracts, shift = c(-1,-1))
bbox(scaled_tracts)

#rotate the first shape 0 degrees
spoly_tracts_rot <- elide.polygonsdf(scaled_tracts, rotate = 0, center = c(0,0))
proj4string(spoly_tracts_rot) <- proj4string(spoly_tracts)

#cluster units
clust_tracts <- clust_units(scaled_tracts, num_clust = 10)

for(i in 1:num_angles){
  #extract rotation angle
  angle_i <- rotate_deg[i]
  
  #define center of the current shape
  shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% 
    as("Spatial") %>% as.data.frame()
  
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
  
  #add border change process, this process takes some time
  block_tracts <- border_change(block_tracts, showoutput = F, quiet = F)
  
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
  
  proj4string(block_tracts) <- proj4string(spoly_tracts)
  
  #make outline of this blockified shape
  outline1 <- outline_one_shp(block_tracts)
  
  #format outline of this shape
  shp1_coords <- outline1 %>% as.matrix()
  
  for(j in 1:nrow(h_dist_mat)){
    #extract rotation angle
    angle_j <- rotate_deg[j]
    
    #define center of the current shape
    shp_centroid <- st_centroid(st_union(st_as_sf(clust_tracts))) %>% 
      as("Spatial") %>% as.data.frame()
    #center the shape at the origin
    clust_tracts <- elide.polygonsdf(clust_tracts, 
                                     shift = c(-shp_centroid$coords.x1[1],
                                               -shp_centroid$coords.x2[1]))
    
    #rotate the tracts
    rot_tracts <- elide.polygonsdf(clust_tracts, rotate = angle_j, center = c(0,0))
    proj4string(rot_tracts) <- proj4string(spoly_tracts)
    
    #blockify tracts
    block_tracts <- blockR(rot_tracts, num_x_blocks = num_x_blocks_run)
    
    #round all blocks
    block_tracts <- round_poly_coords(block_tracts, digits = 5)
    
    #reassign projection
    proj4string(block_tracts) <- proj4string(rot_tracts)
    
    #add border change process, this process takes some time
    block_tracts <- border_change(block_tracts, showoutput = F, quiet = F)
    
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
    
    proj4string(block_tracts) <- proj4string(spoly_tracts)
    
    #make outline of this blockified shape
    outline2 <- outline_one_shp(block_tracts)
    
    #format outline of this shape
    shp2_coords <- outline2 %>% as.matrix()
    
    #distance between two shapes
    h_dist_mat[i,j] <- hausdorff_dist(shp1_coords, shp2_coords)
  
  }
  print(i)
  
}

#format data for plotting
h_dist_mat <- as.data.frame(h_dist_mat)
longData<-melt(h_dist_mat)
colnames(longData)[2] <- "h_dist"
longData$r1 <- rep(1:num_angles, num_angles)*(360/num_angles)
longData$r2 <- longData$r1[order(longData$r1)]

#make part c plots for Figures 9-11
plot_c <- ggplot(longData, aes(x = r1, y = r2)) + 
  geom_raster(aes(fill=h_dist)) + 
  scale_fill_distiller(palette = "Spectral") +
  labs(x="r", y="r'", title= 
         paste("Hausdorff Distance between Rotated Blockified \nShapes, ",
                                   county_run, " County", sep = ""),
       fill = "Hausdorff \nDistance") +theme_light() +
  theme(axis.text=element_text(size=13, angle=0, vjust=0.3),
        axis.title =element_text(size=16),
        axis.text.y=element_text(size=16),
        legend.text=element_text(size=11),
        legend.title = element_text(size=13),
        plot.title=element_text(size=17))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_continuous(expand = expansion(mult = c(0, 0)))

#preview of plots for part c figures
plot_c

cairo_pdf(file = paste(getwd(), "/output/images/",county_run,"_figurec.pdf", sep =""),
          width = 7, 
          height = 5) 
plot_c
dev.off()


###
### Selberg plots
###
county_int <- "Warren"
load(file = paste(getwd(), "/data/plot_data/hdist_plotb_", county_int, ".Rdata",
                  sep = ""))
df_plot2_warren <- full_df
df_plot2_warren$county_int <- county_int

county_int <- "Durham"
load(file = paste(getwd(), "/data/plot_data/hdist_plotb_", county_int, ".Rdata",
                  sep = ""))
df_plot2_durham <- full_df
df_plot2_durham$county_int <- county_int

county_int <- "Middlesex"
load(file = paste(getwd(), "/data/plot_data/hdist_plotb_", county_int, ".Rdata",
                  sep = ""))
df_plot2_midd <- full_df
df_plot2_midd$county_int <- county_int


full_selb_eval_data <- rbind(full_selb_eval_fun(df_plot2=df_plot2_warren,
                                                county_int="Warren", num_test=150),
                             full_selb_eval_fun(df_plot2=df_plot2_durham, 
                                                county_int="Durham", num_test=150),
                             full_selb_eval_fun(df_plot2=df_plot2_midd, 
                                                county_int="Middlesex", num_test=150))
full_selb_eval_data$delta <- as.numeric(full_selb_eval_data$delta)
full_selb_eval_data$bound <- as.numeric(full_selb_eval_data$bound)

p4 <- ggplot(data = full_selb_eval_data) +
  geom_line(aes(x=delta, y=bound, col = county_int, linetype=bound_type), linewidth = 1.5) +
  labs(x=expression(delta), y = "Bound", linetype = "Bound Type", col = "County")+
  scale_colour_brewer(palette = "Set2") + theme_light()

p4

#Reproduce Figure 12 in the paper
cairo_pdf(file = paste(getwd(), "/output/images/figure12.pdf", sep =""),
          width = 7, 
          height = 5) 
p4
dev.off()

