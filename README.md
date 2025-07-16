# BlockR: An Areal Spatial Anonymization and Visualization Tool

Citation to full paper: Haensch, Anna, and Claire Kelling. "BlockR: An Areal Spatial Anonymization and Visualization Tool." *Transactions in GIS* 29.4 (2025): e70070.

This repository contains the code necessary to implement the BlockR method
to anonymize areal data. 

First, we describe the helper files which includes functions to implement this method that we created. 

- **00_bg_functions.R**: This file contains the vast majority of the functions that we've written/adapted to support the BlockR method. This includes functions called `blockR()`, `border_cells_fn()`, `clust_units()`, `connect_islands()`, `full_selb_eval_fun()`, `make_outlines()`, `make_outlines()`, `make_poly_new_coord()`, `outline_one_shp()`, `rotateProj()`, `round_poly_coords()`, `round_poly_coords_sub()`, `samp_prob()`, `selb_bound()`, `st_queen()`, and `st_rook()`. We will add descriptions of these functions soon.
- **00_border_change.R**: We created a function to randomly permute the border of the rotated and blockified shape, as described in our paper. This function is called `border_change()`.

Next, we describe the helper files that were adapted from other packages, including some R packages that have been deprecated. 

- **00_maptools_functions.R**: Adaptation of `maptools` R package functions because some functions have been deprecated (including the `elide()` family of functions, `rotateCoords()`, `scaleCoords()`, and `unionSpatialPolygons()` adapted to our naming of `unionSpatialPolygons_pck()`).
- **00_surveillance_functions.R**: Adaptation of surveillance R package function to detect shapes that are at the border of a polygon, (including `polyAtBorder()` adapted to our naming of `polyAtBorder_cust()`).
- **00_rgeos_functions.R**: Adaptation of `rgeos` R package functions as some functions have been deprecated (including `gUnaryUnion()`, `gBuffer_block()`, `rgeosStatus()`).


Finally, we describe the files to implement the method, as is shown in the paper.

- **01_paper_figures.R**: Reproduces the analysis in the paper shown in Figures 1, 4, 6, and 7.
- **02_blockR_reproduce.R**: Reproduces the analysis in the paper shown in Figures 9, 10, 11.
- **03_images_for_classification.R**: Creating images for the image classification analysis, summarized in Table 1. These include both random and non-random coloring of the blockified units.
- **04_blockR_illustration.R**: The previous files are used to reproduce the analysis in our paper. *If you are looking to utilize the BlockR method on a new place (e.g. county), we would suggest adapting this script for the applicable geographic area.* 


Code to implement our classification analysis (Table 1) is in a Google Colab notebook, which can be accessed by [this link]( https://drive.google.com/drive/folders/1ULYRBd1mSrV_u6g7fgZdKSAqTWsTnnA2?usp=sharing).


We welcome suggestions for improvements to our approach. If you utilize our method in your work, we would love to know. Please reach out to us at ckelling at carleton dot edu.
