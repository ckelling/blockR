###
### BlockR Helper Functions
###


#make outline of just one shape
outline_one_shp <- function(shp2, return_shp = F){
  
  #extract window of tract/shape data
  shp2_window <- st_transform(st_union(st_as_sf(shp2)), 
                              crs = proj4string(shp2)) %>% as("Spatial") %>%
    suppressMessages()

  new_shp <- unionSpatialPolygons_pck(shp2, rep(1,length(shp2))) %>%
    suppressMessages()
  
  if(round(area(new_shp)/area(shp2_window), 7) == 1){
    shp2_coords <- new_shp@polygons[[1]]@Polygons[[1]]@coords
    
    #test area
    #reform polygon
    new_shp_poly <- st_polygon(list(shp2_coords)) %>% as("Spatial")
    
    if(suppressWarnings(round(area(new_shp_poly)/area(new_shp),7) == 1)){
      #only complete if the shape is representative
      shp2_coords <- as.data.frame(shp2_coords)
      return(shp2_coords)
      break
    }
  }
  
  
  #need to convert all of the polygons within the window to a spatialpolygons object
  first_poly <- new_shp@polygons[[1]]@Polygons[[1]]
  
  #first polygon in the correct format
  first_poly_sp <- SpatialPolygons(list(Polygons(list(first_poly),1)))
  
  if(length(new_shp@polygons[[1]]@Polygons) > 1){
    for(i in 2:length(new_shp@polygons[[1]]@Polygons)){
      next_poly <- new_shp@polygons[[1]]@Polygons[[i]]
      
      next_poly_sp <- SpatialPolygons(list(Polygons(list(next_poly),1)))
      
      first_poly_sp <- terra::union(first_poly_sp, next_poly_sp)
    }

    shp2_window <- unionSpatialPolygons_pck(first_poly_sp,
                                            IDs = rep(1, length(first_poly_sp)))
  }else{
    #first_poly_sp is invalid is because of tails, not an issue though because
    #the tail is incorporated into the window
    shp2_window <- unionSpatialPolygons_pck(first_poly_sp,
                                            IDs = rep(1, length(first_poly_sp)))
  }
  
  #keep unrounded version
  shp2_window_unr <- shp2_window
  
  proj4string(first_poly_sp) <- proj4string(shp2)

    if(length(shp2_window@polygons[[1]]@Polygons)>1){
      # if the window still has multiple polygons, 
      #     round the coordinates and re-merge
      shp2_window <- round_poly_coords_sub(shp2_window, digits = 4)
      
      #recombine
      shp2_window <- unionSpatialPolygons_pck(shp2_window, 
                               rep(1, length(shp2_window)))
      
      #store coordinates
      shp2_coords <- shp2_window@polygons[[1]]@Polygons[[1]]@coords %>% 
        as.data.frame()
      
      if(length(shp2_window@polygons[[1]]@Polygons)>1){
        # if the window still has multiple polygons, round the coordinates 
        #        and re-merge
        shp2_window <- round_poly_coords_sub(shp2_window, digits = 3)
        
        #recombine
        shp2_window <- unionSpatialPolygons_pck(shp2_window, 
                                                rep(1, length(shp2_window)))
        
        #store coordinates
        shp2_coords <- shp2_window@polygons[[1]]@Polygons[[1]]@coords %>% 
          as.data.frame()
      }
    }else{
         shp2_coords <- shp2_window@polygons[[1]]@Polygons[[1]]@coords %>% 
           as.data.frame()
    }
  
  if(return_shp == T){
    return(shp2_window_unr)
  }else{
    return(shp2_coords)
  }
}

#find bordering cells of rasterized shape
border_cells_fn <- function(block_tracts){
  #store in separate objects for each shape
  shp2_coords <- outline_one_shp(block_tracts)
  
  #check for mismatch, convert to 
  if(suppressWarnings(area(as(st_polygon(list(as.matrix(shp2_coords))), "Spatial")) < 
     0.99* sum(area(block_tracts)))){
    shp2_coords <- outline_one_shp(block_tracts, return_shp = T)
  }
  
  #make an outline
  if(class(shp2_coords) == "data.frame"){
    shp2_coords <- as.matrix(shp2_coords)
    border <- st_as_sf(as(st_polygon(list(shp2_coords)), "Spatial"))
  }else{
    border <- st_as_sf(shp2_coords)
  }
  st_bbox(border)
  
  #convert to sf object
  block_sf <- st_as_sf(block_tracts)
  
  #see if polygons are on the border
  on_border_ind <- polyAtBorder_cust(as(block_sf, "Spatial"),
                                     W=as(border, "Spatial"), method = "sf")

  
  #subset to just those on the border
  border_shapes <- block_tracts[which(on_border_ind),]
  
  return(border_shapes)
}


#function for the selberg bound
selb_bound <- function(delta_mean, var){
  round(delta_mean^2/(delta_mean^2 + var),3)
}

#function for the sample prob to compare to selberg bound
samp_prob <- function(one_minus_delta_mean_vec, results_df){
  length(which(results_df$h_dist > (one_minus_delta_mean_vec)))/nrow(results_df)
}


full_selb_eval_fun <- function(df_plot2, county_int, num_test){
  #plot over multiple values of delta
  delta_vec <- seq(0,0.9,length=num_test)
  one_minus_delta_mean_vec <- (1-delta_vec)*mean(df_plot2$h_dist) #(1-delta)*mu_A
  delta_mean_vec <- (delta_vec)*mean(df_plot2$h_dist) #(delta)*mu_A
  selb_vec <- rep(NA, num_test)
  
  #calculate selberg bound
  for(i in 1:num_test){
    selb_vec[i] <- selb_bound(delta_mean = delta_mean_vec[i], 
                              var = var(df_plot2$h_dist))
  }

  #calculate sample probability
  samp_vec <- rep(NA, num_test)
  for(i in 1:num_test){
    samp_vec[i] <- samp_prob(one_minus_delta_mean_vec= one_minus_delta_mean_vec[i], 
                             results_df= df_plot2)
  }

  #store data
  full_eval_df <- cbind(c(delta_vec, delta_vec),
                        c(selb_vec, samp_vec),
                        c(rep("Selberg Bound", num_test),
                          rep("Sample Prop", num_test))) %>% as.data.frame()
  colnames(full_eval_df) <- c("delta", "bound", "bound_type")
  full_eval_df$county_int <- rep(county_int, nrow(full_eval_df))
  
  return(full_eval_df)
}

#function to make a new spatial polygon
make_poly_new_coord <- function(grid, new_shape_coords, shps){
  #format so that it can be transformed into Spatial points 
  new_shape_coords <- new_shape_coords %>% as.data.frame()
  colnames(new_shape_coords) <- c("lon", "lat")
  coordinates(new_shape_coords) <- ~lon + lat
  
  # grid is a data.frame. To change it to a spatial data set we have to
  new_shape_coords <- SpatialPoints(new_shape_coords, proj4string = 
                                      CRS(proj4string(shps)))
  
  #append to grid
  grid2 <- union(grid, new_shape_coords) 
  
  #convert to polygons
  sp_pix <- SpatialPixelsDataFrame(points = grid2@coords, data = 
                                     as.data.frame(1:nrow(grid2@coords))) %>%
    as("SpatialPolygonsDataFrame")
  
  proj4string(sp_pix) <- proj4string(shps)
  
  ###
  ### Overlay onto the block groups
  ###
  st_pix <- st_as_sf(sp_pix)
  st_cbg <- st_as_sf(shps)
  st_pix_over <- st_join(st_pix, st_cbg, largest = TRUE)
  sp_pix_over <- as(st_pix_over, "Spatial")
  
  return(sp_pix_over)
}


#function to connect islands
connect_islands <- function(grid, sp_pix_over, shps){
  #detect where the island is located
  prob_ind <- which(colSums(nb2mat(poly2nb(sp_pix_over), zero.policy = T))==0)
  
  #connect the island by adding a new shape: try adding in all four corners 
  #and if it does not work, try again.
  if(length(prob_ind) == 1){
    #try to reconnect the island
    #new set of coordinates to try
    grid_x_res <- bbox(sp_pix_over[prob_ind,])[1,2] - bbox(sp_pix_over[prob_ind,])[1,1]
    grid_y_res <- bbox(sp_pix_over[prob_ind,])[2,2] - bbox(sp_pix_over[prob_ind,])[2,1]
    
    #try a new shape
    #add to the right
    new_shape_coords <- coordinates(grid[prob_ind,])
    new_shape_coords[1] <- new_shape_coords[1] + grid_x_res

    sp_pix_over2 <- make_poly_new_coord(grid, new_shape_coords, shps)
    
    #if this works, you're done. If it doesn't try adding shapes on other parts
    if(length(which(colSums(nb2mat(poly2nb(sp_pix_over2), zero.policy = T))==0)) != 0){
      new_shape_coords <- coordinates(grid[prob_ind,])
      new_shape_coords[2] <- new_shape_coords[2] + grid_y_res
      sp_pix_over2 <- make_poly_new_coord(grid, new_shape_coords, shps)
      if(length(which(colSums(nb2mat(poly2nb(sp_pix_over2), zero.policy = T))==0)) != 0){
        new_shape_coords <- coordinates(grid[prob_ind,])
        new_shape_coords[2] <- new_shape_coords[2] - grid_y_res
        sp_pix_over2 <- make_poly_new_coord(grid, new_shape_coords, shps)
        if(length(which(colSums(nb2mat(poly2nb(sp_pix_over2), zero.policy = T))==0)) != 0){
          new_shape_coords <- coordinates(grid[prob_ind,])
          new_shape_coords[1] <- new_shape_coords[1] + grid_x_res
          sp_pix_over2 <- make_poly_new_coord(grid, new_shape_coords, shps)
          if(length(which(colSums(nb2mat(poly2nb(sp_pix_over2), zero.policy = T))==0)) != 0){
            print("Attempt at connecting islands did not work.")
            plot(sp_pix_over)
            stop()
          }
        }
      }
    }
    
    return(sp_pix_over2)
  }else{
    print("There are more than two islands, try a different resolution?")
  }
  
}


#function to cluster tracts
clust_units <- function(shps, num_clust, census_var = "ALAND"){
  num_clust <- num_clust
  rook_w <- rook_weights(st_as_sf(shps))
  data <- shps@data[,c(census_var)] %>% as.data.frame()
  shp_clusters <- rgeoda::skater(k = num_clust, 
                                  w = rook_w, 
                                  df = data)
  
  #cluster the tracts
  shps <- as(shps, "SpatialPolygons")
  shp_tracts_union <- unionSpatialPolygons_pck(shps,
                                                IDs = shp_clusters$Clusters)
  shps <- as(shp_tracts_union, "SpatialPolygonsDataFrame")
  shps@data$GEOID <- paste("GEO", 1:num_clust, sep = "")
  
  return(shps)
}


#functions to rotate tracts
#Source: https://rpubs.com/geospacedman/rotatespatial
rotateProj = function(spobj, rot_angle) {
  # get bounding box as spatial points object
  boxpts = SpatialPoints(t(bbox(spobj)), proj4string = CRS(proj4string(spobj)))
  # convert to lat-long
  boxLL = bbox(spTransform(boxpts, CRS("+init=epsg:4326")))
  # find the centre
  llc = apply(boxLL, 1, mean) %>% as.numeric()
  # construct the proj4 string
  prj = paste0("+proj=omerc +lat_0=", llc[2], " +lonc=", llc[1], " +alpha=",
               rot_angle,
               " +gamma=0.0 +k=1.000000 +x_0=0.000 +y_0=0.000 +ellps=WGS84 +units=m ")
  # return as a CRS:
  CRS(prj)
  
  # transform
  shp_rotate = spTransform(spobj, CRS(prj))
}



#function to blockify tracts
blockR <- function(shps, num_x_blocks){
  suppressMessages(sf_use_s2(FALSE))
  
  #Create a lattice (grid) over the surface
  grid_num_x <- num_x_blocks
  
  #need to calculate number of y blocks to make them square 
  #so that they have the same width and height (approx)
  bbox_width <- bbox(shps)[1,2] - bbox(shps)[1,1]
  bbox_height <- bbox(shps)[2,2] - bbox(shps)[2,1]
  
  block_width <- bbox_width/grid_num_x
  grid_num_y <- round(bbox_height/block_width)
  
  if(grid_num_y < grid_num_x){
    grid_num_y <- grid_num_x
    
    block_height <- bbox_height/grid_num_y
    grid_num_x <- round(bbox_width/block_height)
  }else{
    grid_num_x <- grid_num_x + 1
  }
  
  #need to add a little bit of wiggle room here so the points right off the 
  #    border don't get deleted
  grid <- expand.grid(lon = seq(bbox(shps)[1,1]+0.01,bbox(shps)[1,2]+0.01,
                                length.out = grid_num_x),
                      lat = seq(bbox(shps)[2,1]+0.01,bbox(shps)[2,2]+0.01,
                                length.out = grid_num_y))
  
  #need to be exactly square grids, 7/21 change to digits = 5
  block_width_grid <- round((range(grid[,1])[2]-range(grid[,1])[1])/(grid_num_x-1),
                            digits = 5)
  
  e <- as(raster::extent(round(bbox(shps), digits = 5)), "SpatialPolygons") %>% 
    st_as_sf()
  
  grd_lrg <- st_make_grid(e, cellsize = c(block_width_grid, block_width_grid)) %>%
    st_as_sf() %>%
    as("Spatial")
  grd_lrg@data <- as.data.frame(1:length(grd_lrg))
  
  sp_pix <- grd_lrg
  
  #remove near empty pixels
  pixel_areas <- suppressWarnings(area(sp_pix))
  
  #empty pixels
  empty_pix <- which(pixel_areas < range(pixel_areas)[2]/2)
  
  if(length(empty_pix) > 0){
    save_data <- sp_pix@data
    sp_pix <- sp_pix[-empty_pix,]
    
    #overwrite data to match dimension
    sp_pix@data <- save_data
  }
  
  proj4string(sp_pix) <- proj4string(shps)
  
  ###
  ### Overlay onto the block groups
  ###
  st_pix <- st_as_sf(sp_pix)
  st_cbg <- st_as_sf(shps)
  st_pix_over <- suppressMessages(suppressWarnings(st_join(st_pix, st_cbg, 
                                                           largest = TRUE)))
  sp_pix_over <- as(st_pix_over, "Spatial")
  
  #remove those outside of the shape
  if(length(which(is.na(sp_pix_over@data$GEOID)))>0){
    sp_pix_over <- sp_pix_over[-which(is.na(sp_pix_over@data$GEOID)),]
  }
  
  #are all of them covered?
  length(shps) == length(unique(sp_pix_over$GEOID))
  
  
  
  #test to see if there are islands
  if(length(which(colSums(nb2mat(suppressMessages(poly2nb(sp_pix_over)), 
                                 zero.policy = T))==0))>0){
    print("There are islands in the blockified shape. Attempting to reconnect.")
    
    sp_pix_over <- connect_islands(grid, sp_pix_over, shps)
  }
  
  return(sp_pix_over)
  
}


#source of functions: https://spatialanalysis.github.io/handsonspatialdata/contiguity-based-spatial-weights.html
st_rook <- function(a, b = a) st_relate(a, b, pattern = "F***1****")
st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")

make_outlines <- function(shp1, shp2){
  
  #extract window of tract/shape data
  #window of both shapes
  shp1_window <- st_transform(st_union(st_as_sf(shp1)), 
                              crs = proj4string(shp1)) %>% as("Spatial")
  shp2_window <- st_transform(st_union(st_as_sf(shp2)),
                              crs = proj4string(shp2)) %>% as("Spatial")
  
  #extract points from both shape outlines
  shp1_coords <- shp1_window@polygons[[1]]@Polygons[[1]]@coords
  
  #need to convert all of the polygons within the window to a spatialpolygons object
  first_poly <- shp2_window@polygons[[1]]@Polygons[[1]]
  first_poly_sp <- SpatialPolygons(list(Polygons(list(first_poly),1)))
  if(length(shp2_window@polygons[[1]]@Polygons) > 1){
    for(i in 2:length(shp2_window@polygons[[1]]@Polygons)){
      next_poly <- shp2_window@polygons[[1]]@Polygons[[i]]
      next_poly_sp <- SpatialPolygons(list(Polygons(list(next_poly),1)))
      
      first_poly_sp <- terra::union(first_poly_sp, next_poly_sp)
    }
    first_poly_sp <- gBuffer_block(first_poly_sp, byid = T, width = 0)
    shp2_window <- unionSpatialPolygons_pck(first_poly_sp, 
                                            IDs = rep(1, length(first_poly_sp)))
  }else{
    #first_poly_sp is invalid is because of tails, not an issue though because
    #the tail is incorporated into the window
    shp2_window <- unionSpatialPolygons_pck(first_poly_sp, 
                                            IDs = rep(1, length(first_poly_sp)))
  }
  
  #plot(first_poly_sp)
  proj4string(first_poly_sp) <- proj4string(shp2)
  
  #less than 10 points defining the window indicates not a likely reasonable boundary
  if(nrow(shp2_window@polygons[[1]]@Polygons[[1]]@coords) < 10){ 
    shp2_coords <- NULL
    #if it just takes the tail as the only polygon, bind the rest
    for(k in 1:length(shp2_window@polygons[[1]]@Polygons)){
      shp2_coords <- rbind(shp2_coords,
                           shp2_window@polygons[[1]]@Polygons[[k]]@coords)
    }
  }else{
    shp2_coords <- shp2_window@polygons[[1]]@Polygons[[1]]@coords
  }

  return(list(shp1_coords, shp2_coords))
  
}


round_poly_coords <- function(poly, digits = 5){
  #need to convert all of the polygons within the window to a spatialpolygons object
  first_poly <- poly@polygons[[1]]@Polygons[[1]]
  
  #join from (rounded) coordinates, instead of from polys
  first_poly_coords <- first_poly@coords %>% round(digits = digits)
  new_poly <- as(st_polygon(list(first_poly_coords)), "Spatial")

  for(i in 2:length(poly)){
    next_poly <- poly@polygons[[i]]@Polygons[[1]]

    #round coordinates of next poly and bind from coordinates, not polygon
    next_poly_coords <- next_poly@coords %>% round(digits = digits)
    next_poly_from_coords <- as(st_polygon(list(next_poly_coords)), "Spatial")

    #new_poly_keep <- new_poly #delete when fixed
    new_poly <- terra::union(new_poly, next_poly_from_coords)
  }
  
  new_poly <- as(new_poly, "SpatialPolygonsDataFrame")
  new_poly@data <- poly@data
  
  return(new_poly)
}



round_poly_coords_sub <- function(poly, digits = 5){
  #need to convert all of the polygons within the window to a spatialpolygons object
  first_poly <- poly@polygons[[1]]@Polygons[[1]]
  #join from (rounded) coordinates, instead of from polys
  first_poly_coords <- first_poly@coords %>% round(digits = digits)
  new_poly <- as(st_polygon(list(first_poly_coords)), "Spatial")
  
  poly_length <- length(poly@polygons[[1]]@Polygons)
  
  for(i in 2:poly_length){
    next_poly <- poly@polygons[[1]]@Polygons[[i]]
    
    #round coordinates of next poly and bind from coordinates, not polygon
    next_poly_coords <- next_poly@coords %>% round(digits = digits)
    next_poly_from_coords <- as(st_polygon(list(next_poly_coords)), "Spatial")
    
    new_poly <- terra::union(new_poly, next_poly_from_coords)
  }
  
  proj4string(new_poly) <- proj4string(poly)
  
  return(new_poly)
}
