###
### Function to redo border
###

border_change <- function(block_tracts){
  #assign ids to sample from
  block_tracts@data$id <- 1:length(block_tracts)
  block_tracts@data <- block_tracts@data %>% select(-X1.length.grd_lrg.)
  
  #find blocks that are on the border
  border_shapes <- border_cells_fn(block_tracts)
  plot(border_shapes)
  
  #initialize indices that we've tried
  tried_ind <- NULL
  eligible_id <- border_shapes@data$id
  
  #keep a count of successful moves and a threshold of number of desired 
  #        successes
  success_count <- 0
  success_threshold <- 0.5*length(border_shapes)
  k <- length(tried_ind)
  
  #keep track of the shapes that have been moved
  mv_id <- NULL
  
  orig_border2 <- border_shapes
  orig_border <- border_shapes

  block_tracts_orig <- block_tracts
  
  while((length(eligible_id) > 0) & (success_count < success_threshold)){
    border_shapes <- orig_border
    i <- sample(eligible_id, size = 1)
    tried_ind <- c(tried_ind, i) #keep track of indices that have been tried
    k <- length(tried_ind)
    
    #print what unit id we are on
    print(i)
    
    #reform block tracts without this block
    block_tracts_wo <- block_tracts[-which(block_tracts@data$id == i),]
    
    #try to move outer shape so that it preserves adjacency matrix and 
    #   doesn't overlap or form an island
    
    #pick one of the outer shapes
    mv_shp <- block_tracts[which(block_tracts@data$id ==i),]
    #extract coordinates
    mv_shp_coords <- mv_shp@polygons[[1]]@Polygons[[1]]@coords
    
    #remove this shape from existing polygon
    border_shapes <- border_shapes[-which(border_shapes@data$id ==i),]
    
    #calculate height and width
    width <- range(mv_shp_coords[,1])[2] - range(mv_shp_coords[,1])[1]
    height <- range(mv_shp_coords[,2])[2] - range(mv_shp_coords[,2])[1]
    
    #round width and height
    width <- round(width, digits = 5)
    height <- round(height, digits = 5)
    
    #round coordinates
    mv_shp_coords <- round(mv_shp_coords, digits = 5)
    
    #initialize new coordinates
    mv_shp_coords_new <- mv_shp_coords
    
    #move up
    mv_shp_coords_up <- mv_shp_coords
    mv_shp_coords_up[,2] <- mv_shp_coords_up[,2] + height
    #move down
    mv_shp_coords_down <- mv_shp_coords
    mv_shp_coords_down[,2] <- mv_shp_coords_down[,2] - height
    #move left
    mv_shp_coords_left <- mv_shp_coords
    mv_shp_coords_left[,1] <- mv_shp_coords_left[,1] - width
    #move right
    mv_shp_coords_right <- mv_shp_coords
    mv_shp_coords_right[,1] <- mv_shp_coords_right[,1] + width
    #move left up diagonal
    mv_shp_coords_lud <- mv_shp_coords
    mv_shp_coords_lud[,2] <- mv_shp_coords_lud[,2] + height
    mv_shp_coords_lud[,1] <- mv_shp_coords_lud[,1] - width
    #move left down diagonal
    mv_shp_coords_ldd <- mv_shp_coords
    mv_shp_coords_ldd[,2] <- mv_shp_coords_ldd[,2] - height
    mv_shp_coords_ldd[,1] <- mv_shp_coords_ldd[,1] - width
    #move right up diagonal
    mv_shp_coords_rud <- mv_shp_coords
    mv_shp_coords_rud[,2] <- mv_shp_coords_rud[,2] + height
    mv_shp_coords_rud[,1] <- mv_shp_coords_rud[,1] + width
    #move right down diagonal
    mv_shp_coords_rdd <- mv_shp_coords
    mv_shp_coords_rdd[,2] <- mv_shp_coords_rdd[,2] - height
    mv_shp_coords_rdd[,1] <- mv_shp_coords_rdd[,1] + width
    
    #test shape of interests
    all_poss_shapes <- list(mv_shp_coords_down, mv_shp_coords_up, 
                            mv_shp_coords_left, mv_shp_coords_right,
                            mv_shp_coords_lud, mv_shp_coords_ldd,
                            mv_shp_coords_rud, mv_shp_coords_rdd)
    
    #randomly reorder the shapes
    all_poss_shapes <- sample(all_poss_shapes, replace = FALSE)
    
    for(j in 1:length(all_poss_shapes)){
      #extract test shape
      test_shape <- all_poss_shapes[[j]]
      #make into polygon
      new_pol = st_polygon(list(test_shape)) %>% as("Spatial")
      #fix projection
      proj4string(new_pol) <- proj4string(border_shapes)
      
      #calculate overlap area
      overlap_area <- sum(st_area(st_intersection(st_as_sf(new_pol), 
                                                  st_as_sf(block_tracts_wo)))) %>% 
        as.numeric()
      #calculate to see if there is line overlap
      overlap_line <- sum(st_length(st_intersection(st_as_sf(new_pol), 
                                                    st_as_sf(block_tracts_wo)))) %>% 
        as.numeric() %>% round()
      
      
      #test to see if it is an island (overlap area ==0) OR
      # if it is covering another shape (overlap area > 0.1)
      # if it is neither, break the for loop and use this shape
      if(overlap_line != 0 & overlap_area < 0.1){
        
        #need to also test neighborhood matrix
        new_pol <- new_pol %>% as("SpatialPolygonsDataFrame")
        new_pol@data$GEOID <- mv_shp@data$GEOID
        new_pol@data$id <- mv_shp@data$id
        
        #check neighborhood matrix
        all_shapes <- terra::union(new_pol, block_tracts_wo)
        
        #before calculating all neighborhood matrices, re-round coordinates, added 7/21
        all_shapes <- round_poly_coords(all_shapes, digits = 5)
        
        #reassign projection
        proj4string(all_shapes) <- proj4string(block_tracts_wo)
        
        #store shape
        all_shapes_neighb <- all_shapes
        
        if(length(which(is.na(all_shapes_neighb@data$id)))>0){
          #delete repetitive shapes
          all_shapes_neighb <- all_shapes_neighb[-which(is.na(all_shapes_neighb@data$id)),]
          #remove unnecessary columns
          if(ncol(all_shapes_neighb@data) > 3){all_shapes_neighb@data <- 
            all_shapes_neighb@data[,1:3]}
        }
        
        poly_by_geoid_new <- unionSpatialPolygons_pck(all_shapes_neighb,
                                                      all_shapes_neighb@data$GEOID)
        poly_mat_new <- st_queen(st_as_sf(poly_by_geoid_new)) %>% 
          as.matrix() %>% round()
        
        poly_by_geoid_orig <- unionSpatialPolygons_pck(block_tracts_orig, 
                                                       block_tracts_orig@data$GEOID)
        poly_mat_orig <- st_queen(st_as_sf(poly_by_geoid_orig)) %>% 
          as.matrix() %>% round()
      
        #test to see if the shape is disconnected now
        geoid_poly <- all_shapes_neighb[which(all_shapes_neighb@data$GEOID == 
                                                new_pol$GEOID),]
        geoid_old <- block_tracts_orig[which(block_tracts_orig@data$GEOID ==
                                               new_pol$GEOID),]
        
        if(length(geoid_poly)>1){
          geoid_nb <- poly2nb(geoid_poly, queen = F)
          card_nb <- sum(card(geoid_nb))
          if(card_nb > 0){
            #test for island using st, not nb2mat (nb2mat occasionally misidentifies rook as queen)
            rook_st_nums_geoid <- st_rook(st_as_sf(geoid_poly)) %>%
              as.matrix() %>% round()

            #record if this creates islands
            shp_island_bin <- length(which(rowSums(rook_st_nums_geoid) == 0)) == 0
          }else{
            #in this case, there are only islands
            shp_island_bin <- F
          }
        }else{
          #if there is only one shape, it is okay for this test
          shp_island_bin <- T
        }
        
        #only need to perform this test if test for islands is satisfied
        if(shp_island_bin == T & length(geoid_poly) > 1){
          #test to see if there are multiple components to the shape
          poly_mat_geoid_q <- st_queen(st_as_sf(geoid_poly)) %>% 
            as.matrix() %>% round()
          
          g_geoid <- graph_from_adjacency_matrix(round(poly_mat_geoid_q), 
                                                 mode="undirected", weighted=NULL)
          component_num_geoid <- components(g_geoid)$no
        }else if(length(geoid_poly) == 1){
          component_num_geoid <- 1
        }else{
          #otherwise, make the test fail, just like the island test failed
          component_num_geoid <- 2
        }
        
        
        #test to see if the neighborhood matrices are identical
        identical_test <- identical(poly_mat_new, poly_mat_orig)
        
        #test for island using st, not nb2mat (nb2mat occasionally misidentifies 
        #      rook as queen)
        rook_st_nums <- st_rook(st_as_sf(all_shapes_neighb)) %>% 
          as.matrix() %>% round()
        num_island <- length(which(rowSums(rook_st_nums) == 0))
        
        #check to make sure that the number of *rook* islands has not
        #   increased (not that it is 0)
        rook_st_orig <- st_rook(st_as_sf(block_tracts_orig)) %>%
          as.matrix() %>% round()
        num_island_orig <- length(which(rowSums(rook_st_orig) == 0))
        
        #adjacency matrix from sf, not poly2nb
        cell_nb_mat_q <- st_rook(st_as_sf(all_shapes_neighb)) %>% 
          as.matrix() %>% round()
        cell_nb_mat_orig <- st_rook(st_as_sf(block_tracts_orig)) %>% 
          as.matrix() %>% round()

        #make igraph on cell nb mat and test to see if there are 
        #     separate components
        g <- graph_from_adjacency_matrix(round(cell_nb_mat_q), 
                                         mode="undirected", weighted=NULL)
        component_num <- components(g)$no
        g_orig <- graph_from_adjacency_matrix(round(cell_nb_mat_orig), 
                                              mode="undirected", weighted=NULL)
        component_num_orig <- components(g_orig)$no
        
        #are the number of rook connected components the same as the original?
        rook_comp_test <- (component_num <= component_num_orig)
        
        
        #check to see if this new shape creates holes/lakes
        all_shapes_fill <- round_poly_coords(all_shapes, digits = 3)
        all_shp_area <- sum(area(all_shapes_fill))
        
        #combine all of the polygons, fill the holes, then calculate area
        #threshold set so that all holes are filled
        fill_poly <- fill_holes(unionSpatialPolygons_pck(all_shapes_fill, 
                                                         rep(1, length(all_shapes_fill))), 
                                threshold= sum(area(all_shapes_fill)))
        fill_hole_area <- area(fill_poly)
        
        #how many polygons are there? If there are more than one, can indicate a lake
        #solution: reject the case that this polygon has more than one polygon
        num_union_poly <- length(fill_poly@polygons[[1]]@Polygons)
        
        #need to keep shapes that just have a tail
        row_rook <- st_rook(st_as_sf(all_shapes_fill)) %>% 
          as.matrix() %>% round() %>% rowSums()
        row_queen <- st_queen(st_as_sf(all_shapes_fill)) %>% 
          as.matrix() %>% round() %>% rowSums()
        if(length(which(row_rook == 0)) > 0){
          row_queen_tail <- row_queen[which(row_rook == 0)]
          row_rook_tail <- row_rook[which(row_rook == 0)]
          just_a_tail <- 0
          for(i in 1:length(row_queen_tail)){
            if(row_rook_tail[i] == 0 & row_queen_tail[i] != 0){
              #testing to see if there is a tail to the shape
              just_a_tail <- just_a_tail + 0
            }
          }
          #if all of the rook disconnected shapes are tails, it is satisfied
          if(just_a_tail == 0){
            num_union_poly <- 1
          }
        }
        
        #check ratio of areas, want this to be less than 0.1%
        no_hole_ratio <- fill_hole_area/all_shp_area
        
        #if the overlap tests are satisfied, check the neighborhood matrices and
        #    look for islands
        if(identical_test & num_island  <= num_island_orig & shp_island_bin & 
           no_hole_ratio < 1.001 & rook_comp_test & #component_num ==1 & 
           component_num_geoid == 1 & num_union_poly == 1){
          #if all conditions are satisfied, consider it a success
          success_count <- success_count + 1
          #if there is a match, make it part of border shapes
          orig_border <- border_cells_fn(all_shapes)
          border_shapes <- border_cells_fn(all_shapes)
          store_prev_it <- block_tracts
          block_tracts <- all_shapes
          print(paste("Success count is ", success_count, " out of ", 
                      success_threshold, ". (unit ", i, "), shape", j, sep = ""))
          
          #store ids for those that have moved
          mv_id <- c(mv_id, i)
          
          #update eligible ids, if cells are no longer on the border
          no_longer_border <- which(!(eligible_id %in% border_shapes$id))
          if(length(no_longer_border) > 0) eligible_id <- 
            eligible_id[-no_longer_border]
          
          #newly on border, and haven't already tried
          new_border_ind <- which(!(border_shapes$id %in% tried_ind))
          if(length(new_border_ind) > 0) eligible_id <- 
            c(eligible_id, border_shapes$id[new_border_ind]) %>% unique()
          
          #take out ids that are NAs
          if(length(which(is.na(eligible_id))) > 0){
            eligible_id <- eligible_id[-which(is.na(eligible_id))]
          }
          
          #plot new and old shapes
          plot(st_as_sf(block_tracts_orig)["GEOID"])
          plot(st_as_sf(block_tracts)["GEOID"])
          break
        }
      }
      if(j == length(all_poss_shapes)){
        print(paste("No match found for try ", k, " (unit ", i,
                    ").", sep = ""))
        
        mv_shp_df <- as.data.frame(mv_shp_coords)
        no_match <- ggplot()+ geom_sf(data = st_as_sf(block_tracts), aes(fill = GEOID)) + 
          geom_point(data=mv_shp_df, aes(x=x,y=y), col = "red", size=6, inherit.aes = T)+theme_void()
        print(no_match)
      }
    }
    
    #subset to only eligible indices
    if(length(tried_ind)>0){eligible_id <- eligible_id[-which(eligible_id %in% tried_ind)]}
    
    #take out NA geoid
    na_geoid <- length(which(is.na(block_tracts@data$GEOID)))
    if(na_geoid>0){block_tracts <- block_tracts[-which(is.na(block_tracts@data$GEOID)),]}
    
    print(paste("Try #:", length(tried_ind)))
  }
  
  return(block_tracts)
}


