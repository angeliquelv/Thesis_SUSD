# In this script Geotiff files are stacked, cropped, low quality pixels are filtered out
# Input: Geotiff files of ndvi van quailty flag for each date and three tiles
# Output: stacked data

rm(list = ls())  # clear workspace

# Install packages and load libraries -------------------------------------

# install packages if required
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}
if(!require(lubridate)){install.packages("lubridate")}

# load libraries
library(rgdal)
library(raster) # package for raster manipulation
library(lubridate) # package to handle dates more easily

# Set working directory ---------------------------------------------------

workdirNDVI <- "G:/Thesis_data/NDVI"
workdirNDVILandsat5_7 <- "G:/Thesis_data/NDVI/Landsat5_7"
workdirQFLAGLandsat5_7 <- "G:/Thesis_data/NDVI/Landsat5_7/QFLAG"
workdirNDVILandsat8 <- "G:/Thesis_data/NDVI/Landsat8"
workdirQFLAGLandsat8 <- "G:/Thesis_data/NDVI/Landsat8/QFLAG"
workdir_GIS_Ounila = "G:/Thesis_data/Shapefiles"

# Open reprojected shapefile ----------------------------------------------

setwd(workdir_GIS_Ounila)

Ounila_catchment_reprojected <- readOGR(dsn=workdir_GIS_Ounila,layer="Ounila_catchment_reprojected")

# Unzip and untar tarballs ------------------------------------------------

# For Landsat 5 and 7 data
setwd(workdirNDVILandsat5_7)

tgz_files <- list.files(path = workdirNDVILandsat5_7, pattern = '.tar.gz$') # list of tarballs in data directory

for(i in 1:length(tgz_files)) {
  # untar tarballs
  untar(tgz_files[i])
}

# For Landsat 8 data
setwd(workdirNDVILandsat8)

tgz_files <- list.files(path = workdirNDVILandsat8, pattern = '.tar.gz$') # list of tarballs in data directory

for(i in 1:length(tgz_files)) {
  # untar tarballs
  untar(tgz_files[i])
}

# Now move all quality flag files to another directory workdirQFLAG (I did this using Windows Explorer)

# Create rasterstack of Landsat 5 and 7 NDVI data -------------------------

setwd(workdirNDVILandsat5_7)

# create a list of all the NDVI files in the working directory
NDVI_tif_files <- list.files(path = workdirNDVILandsat5_7, pattern = '.tif$') # list of GeoTIFF in data directory
splitted_names<-strsplit(NDVI_tif_files,"_") # split names by _
# for loop that isolates the dates from the filenames
dates_list <- vector()
for (i in 1:length(NDVI_tif_files)) {
  dates_list[i]<-splitted_names[[i]][4]
}

# save the original file names
saveRDS(NDVI_tif_files,file="orignal_names.Rda")

# create list of all the unique dates
dates_unique <- unique(dates_list)
# sort in chronological order
dates_sorted <- sort(dates_unique)

# The names of the files are structured as follows:
# LE07_L1GT_201038_20001025_20170209_01_T2_sr_ndvi
# LXSS_LLLL_PPPRRR_YYYYMMDD_yyyymmdd_CX_TX_prod_band.ext
# L Landsat
# X Sensor ("E" = ETM+; "T" = TM)
# SS Satellite ("07" = Landsat 7; "05" = Landsat 5; "04" = Landsat 4)
# LLLL Processing correction level ("L1TP" = Precision Terrain; "L1GT" =
#                                      Systematic Terrain; "L1GS" = Systematic)
# PPP Path
# RRR Row
# YYYY Year of acquisition
# MM Month of acquisition
# DD Day of acquisition
# yyyy Year of L1 processing
# mm Month of L1 processing
# dd Day of L1 processing
# CX Collection number ("01")
# TX Collection category ("RT" = Real Time; "T1" = Tier 1; "T2" = Tier 2)
# prod Product, such as "toa", "bt", or "sr"
# band Band, such as "band<1-7>", "qa," or spectral index.
# ext File format extension, such as "tif", "tfw", "xml", "hdf", "hdr", "nc", or "img"

# Rename the files, by removing the date of processing (to avoid mixing up in the for loop)
new_names_tif_files <- vector()
for (i in 1:length(NDVI_tif_files)) {
  new_names_tif_files[i] <- paste(splitted_names[[i]][1],splitted_names[[i]][2],splitted_names[[i]][3],splitted_names[[i]][4],splitted_names[[i]][6],splitted_names[[i]][7],splitted_names[[i]][8],splitted_names[[i]][9],sep="_")
}

file.rename(NDVI_tif_files,new_names_tif_files)

# scenes from four different path/row combinations were retrieved (201038, 201039, 202038, 202039)
# Not for each date, all these four scenes were available
# However I do want to have rasters of the same extent, for the operations that follow (cropping, masking out low quality pixels, etc)
# This extent encloses the 4 scenes combined:
extent_scenes<-extent(420000,900000,3200000,3800000)

# # To avoid making the rasterstack unneccesarily large, the data will immediately be cropped to the rough extent of the Ounila watershed:
# Ounila_catchment_rough_extent<-spTransform(Ounila_catchment_reprojected,"+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# Ounila_catchment_rough_extent
# OUTPUT:
# class       : SpatialPolygonsDataFrame
# features    : 1
# extent      : 652374.4, 693084.4, 3434945, 3473015  (xmin, xmax, ymin, ymax)
# crs         : +proj=utm +zone=29 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0
# variables   : 2
# names       :  DN,      area
# value       : 100, 730863000
# rough_extent <- extent (652374.4,)
rough_extent<-extent(652000,694000,3434000,3474000)

# This for loop checks for each date whether there are 1, 2, 3 or 4 files
# The files are then individually converted to a raster format
# Fill value and saturate value (-9999 and 20000) are set to NA, before the tiles are mosaicked together
# This has to be done first because when the tiles are mosaicked a mean is taken on the areas that overlap
# When a pixel in the overlapping areas has value 20000 or -9999 this can result in unrealistic values
# The individual files are then mosaicked together
# They are then expanded to the full extent of 4 scenes combined by filling the raster up with NA
# Then the rasterlayer is cropped to the rough Ounila watershed extent
# Finally the rasterlayer is added to the raster stack and given a name containing the satellite (Landsat 5 or 7) and the date
# Unnecessary files are removed from hard disk and workspace

NDVI_stack_Landsat_5_7 <- stack()
for (i in 1:length(dates_sorted)) {
  list <- list.files(path=workdirNDVILandsat5_7,pattern=paste0(dates_sorted[i],".+.tif$"))
  if (length(list)==1) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list)
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster_scene1,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent) # a temporary file is created, needed when the stack is written to disk
    # add to rasterstack
    NDVI_stack_Landsat_5_7 <- stack(NDVI_stack_Landsat_5_7,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list,"_")
    names(NDVI_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==2) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list[1])
    NDVI_raster_scene2 <- raster(list[2])
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # mosaic raster files together
    NDVI_raster <- mosaic(NDVI_raster_scene1,NDVI_raster_scene2,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent)
    # add to rasterstack
    NDVI_stack_Landsat_5_7 <- stack(NDVI_stack_Landsat_5_7,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(NDVI_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_scene2_saturate_removed.grd","raster_scene2_saturate_removed.gri",
                  "raster_scene2_fill_removed.grd","raster_scene2_fill_removed.gri",
                  "raster_mosaic.grd","raster_mosaic.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==3) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list[1])
    NDVI_raster_scene2 <- raster(list[2])
    NDVI_raster_scene3 <- raster(list[3])
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # mosaic raster files together
    NDVI_raster <- mosaic(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent)
    # add to rasterstack
    NDVI_stack_Landsat_5_7 <- stack(NDVI_stack_Landsat_5_7,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(NDVI_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_scene2_saturate_removed.grd","raster_scene2_saturate_removed.gri",
                  "raster_scene2_fill_removed.grd","raster_scene2_fill_removed.gri",
                  "raster_scene3_saturate_removed.grd","raster_scene3_saturate_removed.gri",
                  "raster_scene3_fill_removed.grd","raster_scene3_fill_removed.gri",
                  "raster_mosaic.grd","raster_mosaic.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==4) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list[1])
    NDVI_raster_scene2 <- raster(list[2])
    NDVI_raster_scene3 <- raster(list[3])
    NDVI_raster_scene4 <- raster(list[4])
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene4 <- mask(NDVI_raster_scene4, NDVI_raster_scene4, filename="raster_scene4_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene4 <- mask(NDVI_raster_scene4, NDVI_raster_scene4, filename="raster_scene4_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # mosaic raster files together
    NDVI_raster <- mosaic(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,NDVI_raster_scene4,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent)
    # add to rasterstack
    NDVI_stack_Landsat_5_7 <- stack(NDVI_stack_Landsat_5_7,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(NDVI_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,NDVI_raster_scene4,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_scene2_saturate_removed.grd","raster_scene2_saturate_removed.gri",
                  "raster_scene2_fill_removed.grd","raster_scene2_fill_removed.gri",
                  "raster_scene3_saturate_removed.grd","raster_scene3_saturate_removed.gri",
                  "raster_scene3_fill_removed.grd","raster_scene3_fill_removed.gri",
                  "raster_scene4_saturate_removed.grd","raster_scene4_saturate_removed.gri",
                  "raster_scene4_fill_removed.grd","raster_scene4_fill_removed.gri",
                  "raster_mosaic.grd","raster_mosaic.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  # print progress
  print(i)
}

writeRaster(NDVI_stack_Landsat_5_7,"NDVI_stack.grd",overwrite=TRUE)

# Create rasterstack of Landsat 5 and 7 QFLAG data ------------------------

setwd(workdirQFLAGLandsat5_7)

# create a list of all the unique dates
QFLAG_tif_files <- list.files(path = workdirQFLAGLandsat5_7, pattern = '.tif$') # list of GeoTIFF in data directory
splitted_names<-strsplit(QFLAG_tif_files,"_") # split names by _

# save the original file names
saveRDS(QFLAG_tif_files,file="orignal_names.Rda")

# Rename the files, by removing the date of processing (to avoid mixing up in the for loop)
new_names_tif_files <- vector()
for (i in 1:length(QFLAG_tif_files)) {
  new_names_tif_files[i] <- paste(splitted_names[[i]][1],splitted_names[[i]][2],splitted_names[[i]][3],splitted_names[[i]][4],splitted_names[[i]][6],splitted_names[[i]][7],splitted_names[[i]][8],splitted_names[[i]][9],sep="_")
}

file.rename(QFLAG_tif_files,new_names_tif_files)

# taking the mean of overlapping values in the quality flag was probably a mistake. 
# this way that data may have gotten compromised. -- is should have applied the quality flag to each 
# individual frame before mosaicking the raster together.
# however effect is not expected to be large. if the pixel had a clear value in both it will still have that
# however, if one of the two pixels in overlapping areas was cloudy it means that that pixel is now filtered out altogher
# whereas otherwise it could still have a value from the other image.
# this tiny misttake may have resulted in slightly less pixels with data around overlapping areas

# The same operations are performed for loop on the quality flag data
QFLAG_stack_Landsat_5_7 <- stack()
for (i in 1:length(dates_sorted)) {
  list <- list.files(path=workdirQFLAGLandsat5_7,pattern=paste0(dates_sorted[i],".+.tif$"))
  if (length(list)==1) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list)
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster_scene1,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_5_7 <- stack(QFLAG_stack_Landsat_5_7,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list,"_")
    names(QFLAG_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==2) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list[1])
    QFLAG_raster_scene2 <- raster(list[2])
    # mosaic raster files together
    QFLAG_raster <- mosaic(QFLAG_raster_scene1,QFLAG_raster_scene2,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_5_7 <- stack(QFLAG_stack_Landsat_5_7,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(QFLAG_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_mosaic.grd","raster_mosaic.gri","raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==3) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list[1])
    QFLAG_raster_scene2 <- raster(list[2])
    QFLAG_raster_scene3 <- raster(list[3])
    # mosaic raster files together
    QFLAG_raster <- mosaic(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_5_7 <- stack(QFLAG_stack_Landsat_5_7,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(QFLAG_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_mosaic.grd","raster_mosaic.gri","raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==4) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list[1])
    QFLAG_raster_scene2 <- raster(list[2])
    QFLAG_raster_scene3 <- raster(list[3])
    QFLAG_raster_scene4 <- raster(list[4])
    # mosaic raster files together
    QFLAG_raster <- mosaic(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,QFLAG_raster_scene4,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_5_7 <- stack(QFLAG_stack_Landsat_5_7,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(QFLAG_stack_Landsat_5_7[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,QFLAG_raster_scene4,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_mosaic.grd","raster_mosaic.gri","raster_extend.grd","raster_extend.gri"))
  }
  # print progress
  print(i)
}
writeRaster(QFLAG_stack_Landsat_5_7,"QFLAG_stack.grd",overwrite=TRUE)

# Create rasterstack of Landsat 8 NDVI data -------------------------------

setwd(workdirNDVILandsat8)

# create a list of all the NDVI files in the working directory
NDVI_tif_files <- list.files(path = workdirNDVILandsat8, pattern = '.tif$') # list of GeoTIFF in data directory
splitted_names<-strsplit(NDVI_tif_files,"_") # split names by _
# for loop that isolates the dates from the filenames
dates_list <- vector()
for (i in 1:length(NDVI_tif_files)) {
  dates_list[i]<-splitted_names[[i]][4]
}

# save the original file names
saveRDS(NDVI_tif_files,file="original_names.Rda")

# create list of all the unique dates
dates_unique <- unique(dates_list)
# sort in chronological order
dates_sorted <- sort(dates_unique)

# Rename the files, by removing the date of processing (to avoid mixing up in the for loop)
new_names_tif_files <- vector()
for (i in 1:length(NDVI_tif_files)) {
  new_names_tif_files[i] <- paste(splitted_names[[i]][1],splitted_names[[i]][2],splitted_names[[i]][3],splitted_names[[i]][4],splitted_names[[i]][6],splitted_names[[i]][7],splitted_names[[i]][8],splitted_names[[i]][9],sep="_")
}

file.rename(NDVI_tif_files,new_names_tif_files)

# This for loop checks for each date whether there are 1, 2, 3 or 4 files
# The files are then individually converted to a raster format
# The individual files are then mosaicked together
# They are then expanded to the full extent of 4 scenes combined by filling the raster up with NA
# Then the rasterlayer is cropped to the rough Ounila watershed extent
# Finally the rasterlayer is added to the raster stack and given a name containing the date and the satellite (Landsat8)
# Unnecessary files are removed from hard disk and workspace
NDVI_stack_Landsat_8 <- stack()
for (i in 1:length(dates_sorted)) {
  list <- list.files(path=workdirNDVILandsat8,pattern=paste0(dates_sorted[i],".+.tif$"))
  if (length(list)==1) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list)
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster_scene1,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent) # a temporary file is created, needed when the stack is written to disk
    # add to rasterstack
    NDVI_stack_Landsat_8 <- stack(NDVI_stack_Landsat_8,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list,"_")
    names(NDVI_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==2) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list[1])
    NDVI_raster_scene2 <- raster(list[2])
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # mosaic raster files together
    NDVI_raster <- mosaic(NDVI_raster_scene1,NDVI_raster_scene2,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent)
    # add to rasterstack
    NDVI_stack_Landsat_8 <- stack(NDVI_stack_Landsat_8,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(NDVI_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_scene2_saturate_removed.grd","raster_scene2_saturate_removed.gri",
                  "raster_scene2_fill_removed.grd","raster_scene2_fill_removed.gri",
                  "raster_mosaic.grd","raster_mosaic.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==3) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list[1])
    NDVI_raster_scene2 <- raster(list[2])
    NDVI_raster_scene3 <- raster(list[3])
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # mosaic raster files together
    NDVI_raster <- mosaic(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent)
    # add to rasterstack
    NDVI_stack_Landsat_8 <- stack(NDVI_stack_Landsat_8,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(NDVI_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_scene2_saturate_removed.grd","raster_scene2_saturate_removed.gri",
                  "raster_scene2_fill_removed.grd","raster_scene2_fill_removed.gri",
                  "raster_scene3_saturate_removed.grd","raster_scene3_saturate_removed.gri",
                  "raster_scene3_fill_removed.grd","raster_scene3_fill_removed.gri",
                  "raster_mosaic.grd","raster_mosaic.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==4) {
    # convert file to raster
    NDVI_raster_scene1 <- raster(list[1])
    NDVI_raster_scene2 <- raster(list[2])
    NDVI_raster_scene3 <- raster(list[3])
    NDVI_raster_scene4 <- raster(list[4])
    # Set saturated pixels to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    NDVI_raster_scene4 <- mask(NDVI_raster_scene4, NDVI_raster_scene4, filename="raster_scene4_saturate_removed.grd",
                               inverse=FALSE, maskvalue=20000, updatevalue=NA)
    # Set pixels with fill value to NA
    NDVI_raster_scene1 <- mask(NDVI_raster_scene1, NDVI_raster_scene1, filename="raster_scene1_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene2 <- mask(NDVI_raster_scene2, NDVI_raster_scene2, filename="raster_scene2_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene3 <- mask(NDVI_raster_scene3, NDVI_raster_scene3, filename="raster_scene3_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    NDVI_raster_scene4 <- mask(NDVI_raster_scene4, NDVI_raster_scene4, filename="raster_scene4_fill_removed.grd",
                               inverse=FALSE, maskvalue=-9999, updatevalue=NA)
    # mosaic raster files together
    NDVI_raster <- mosaic(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,NDVI_raster_scene4,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    NDVI_raster <- extend(NDVI_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    NDVI_raster <- crop(NDVI_raster,rough_extent)
    # add to rasterstack
    NDVI_stack_Landsat_8 <- stack(NDVI_stack_Landsat_8,NDVI_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(NDVI_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(NDVI_raster_scene1,NDVI_raster_scene2,NDVI_raster_scene3,NDVI_raster_scene4,NDVI_raster) # remove unused file from workspace
    file.remove(c("raster_scene1_saturate_removed.grd","raster_scene1_saturate_removed.gri",
                  "raster_scene1_fill_removed.grd","raster_scene1_fill_removed.gri",
                  "raster_scene2_saturate_removed.grd","raster_scene2_saturate_removed.gri",
                  "raster_scene2_fill_removed.grd","raster_scene2_fill_removed.gri",
                  "raster_scene3_saturate_removed.grd","raster_scene3_saturate_removed.gri",
                  "raster_scene3_fill_removed.grd","raster_scene3_fill_removed.gri",
                  "raster_scene4_saturate_removed.grd","raster_scene4_saturate_removed.gri",
                  "raster_scene4_fill_removed.grd","raster_scene4_fill_removed.gri",
                  "raster_mosaic.grd","raster_mosaic.gri",
                  "raster_extend.grd","raster_extend.gri"))
  }
  # print progress
  print(i)
}

writeRaster(NDVI_stack_Landsat_8,"NDVI_stack.grd",overwrite=TRUE)

# Create rasterstack of Landsat 8 QFLAG data ------------------------------

setwd(workdirQFLAGLandsat8)

# create a list of all the unique dates
QFLAG_tif_files <- list.files(path = workdirQFLAGLandsat8, pattern = '.tif$') # list of GeoTIFF in data directory
splitted_names<-strsplit(QFLAG_tif_files,"_") # split names by _

# save the original file names
saveRDS(QFLAG_tif_files,file="orignal_names.Rda")

# Rename the files, by removing the date of processing (to avoid mixing up in the for loop)
new_names_tif_files <- vector()
for (i in 1:length(QFLAG_tif_files)) {
  new_names_tif_files[i] <- paste(splitted_names[[i]][1],splitted_names[[i]][2],splitted_names[[i]][3],splitted_names[[i]][4],splitted_names[[i]][6],splitted_names[[i]][7],splitted_names[[i]][8],splitted_names[[i]][9],sep="_")
}

file.rename(QFLAG_tif_files,new_names_tif_files)

# The same operations are performed for loop on the quality flag data
QFLAG_stack_Landsat_8 <- stack()
for (i in 1:length(dates_sorted)) {
  list <- list.files(path=workdirQFLAGLandsat8,pattern=paste0(dates_sorted[i],".+.tif$"))
  if (length(list)==1) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list)
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster_scene1,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_8 <- stack(QFLAG_stack_Landsat_8,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list,"_")
    names(QFLAG_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==2) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list[1])
    QFLAG_raster_scene2 <- raster(list[2])
    # mosaic raster files together
    QFLAG_raster <- mosaic(QFLAG_raster_scene1,QFLAG_raster_scene2,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_8 <- stack(QFLAG_stack_Landsat_8,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(QFLAG_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_mosaic.grd","raster_mosaic.gri","raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==3) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list[1])
    QFLAG_raster_scene2 <- raster(list[2])
    QFLAG_raster_scene3 <- raster(list[3])
    # mosaic raster files together
    QFLAG_raster <- mosaic(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_8 <- stack(QFLAG_stack_Landsat_8,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(QFLAG_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_mosaic.grd","raster_mosaic.gri","raster_extend.grd","raster_extend.gri"))
  }
  else if (length(list)==4) {
    # convert file to raster
    QFLAG_raster_scene1 <- raster(list[1])
    QFLAG_raster_scene2 <- raster(list[2])
    QFLAG_raster_scene3 <- raster(list[3])
    QFLAG_raster_scene4 <- raster(list[4])
    # mosaic raster files together
    QFLAG_raster <- mosaic(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,QFLAG_raster_scene4,fun=mean,filename="raster_mosaic.grd")
    # expand file to full extent
    QFLAG_raster <- extend(QFLAG_raster,extent_scenes,value=NA,filename="raster_extend.grd")
    # Crop to rough Ounila extent
    QFLAG_raster <- crop(QFLAG_raster,rough_extent)
    # add to rasterstack
    QFLAG_stack_Landsat_8 <- stack(QFLAG_stack_Landsat_8,QFLAG_raster)
    # Make sure that the layer name still shows the date and whether Landsat 5 or 7 was used
    split <- strsplit(list[1],"_")
    names(QFLAG_stack_Landsat_8[[i]]) <- paste0(split[[1]][1],split[[1]][4])
    # remove unnecessary files from workspace and temporary files
    rm(QFLAG_raster_scene1,QFLAG_raster_scene2,QFLAG_raster_scene3,QFLAG_raster_scene4,QFLAG_raster) # remove unused file from workspace
    file.remove(c("raster_mosaic.grd","raster_mosaic.gri","raster_extend.grd","raster_extend.gri"))
  }
  # print progress
  print(i)
}
writeRaster(QFLAG_stack_Landsat_8,"QFLAG_stack.grd",overwrite=TRUE)

# Remove low quality pixels -----------------------------------------------

setwd(workdirNDVILandsat5_7)

# All values that are not 66 in the quality flag stack are masked out in the NDVI stack
NDVI_clean_Landsat_5_7 <- mask(NDVI_stack_Landsat_5_7, QFLAG_stack_Landsat_5_7, filename="NDVI_clean.grd",overwrite=TRUE,
                   inverse=TRUE, maskvalue=66, updatevalue=NA)

setwd(workdirNDVILandsat8)

# All values that are not 322 in the quality flag stack are masked out in the NDVI stack
NDVI_clean_Landsat_8 <- mask(NDVI_stack_Landsat_8, QFLAG_stack_Landsat_8, filename="NDVI_clean.grd",overwrite=TRUE,
                   inverse=TRUE, maskvalue=322, updatevalue=NA)


# Combine Landsat 5 and 7 with Landsat 8 ----------------------------------

setwd(workdirNDVI)

# Now order and combine the two rasterstacks
NDVI_combined <- stack(NDVI_clean_Landsat_5_7,NDVI_clean_Landsat_8)

dates_names <- substr(names(NDVI_combined),5,12)
ordered_dates <- sort(dates_names)

names_position <- vector()
for (i in 1:length(ordered_dates)) {
names_position[i] <- grep(pattern=ordered_dates[i],names(NDVI_combined))
}

NDVI_clean <- NDVI_combined[[names_position]]

writeRaster(NDVI_clean,"NDVI_clean.grd",overwrite=TRUE)


# Apply scale factor ------------------------------------------------------

NDVI_scaled <- NDVI_clean * 0.0001

# Reproject to WGS84 ------------------------------------------------------

NDVI_reprojected <- projectRaster(NDVI_scaled,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",method="ngb",filename="NDVI_reprojected.grd",overwrite=TRUE)
# Warning messages:
# 1: In showSRID(uprojargs, format = "PROJ", multiline = "NO") :
#   Discarded datum WGS_1984 in CRS definition,
# but +towgs84= values preserved
# 2: In showSRID(uprojargs, format = "PROJ", multiline = "NO") :
#   Discarded datum Unknown based on WGS84 ellipsoid in CRS definition,
# but +towgs84= values preserve
# I've reprojected using proj4string which delivered a warning
# The reason for this is that proj4string is being deprecated
# However projectRaster still relies on the Proj4.string library

# Crop and mask to Ounila watershed extent --------------------------------

NDVI <- crop(NDVI_reprojected,Ounila_catchment_reprojected)
NDVI <- mask(NDVI,Ounila_catchment_reprojected)

writeRaster(NDVI,"NDVI.grd",overwrite=TRUE)

# output of > cat(strwrap(gsub(",", ", ", (comment(NDVI@crs)))), sep="\n")
# BOUNDCRS[ SOURCECRS[ GEOGCRS["unknown", DATUM["Unknown based on WGS84 ellipsoid",
# ELLIPSOID["WGS 84", 6378137, 298.257223563, LENGTHUNIT["metre", 1], ID["EPSG", 7030]]],
# PRIMEM["Greenwich", 0, ANGLEUNIT["degree", 0.0174532925199433], ID["EPSG", 8901]],
# CS[ellipsoidal, 2], AXIS["longitude", east, ORDER[1], ANGLEUNIT["degree",
# 0.0174532925199433, ID["EPSG", 9122]]], AXIS["latitude", north, ORDER[2],
# ANGLEUNIT["degree", 0.0174532925199433, ID["EPSG", 9122]]]]], TARGETCRS[ GEOGCRS["WGS 84",
# DATUM["World Geodetic System 1984", ELLIPSOID["WGS 84", 6378137, 298.257223563,
# LENGTHUNIT["metre", 1]]], PRIMEM["Greenwich", 0, ANGLEUNIT["degree", 0.0174532925199433]],
# CS[ellipsoidal, 2], AXIS["latitude", north, ORDER[1], ANGLEUNIT["degree",
# 0.0174532925199433]], AXIS["longitude", east, ORDER[2], ANGLEUNIT["degree",
# 0.0174532925199433]], ID["EPSG", 4326]]], ABRIDGEDTRANSFORMATION["Transformation from
# unknown to WGS84", METHOD["Position Vector transformation (geog2D domain)", ID["EPSG",
# 9606]], PARAMETER["X-axis translation", 0, ID["EPSG", 8605]], PARAMETER["Y-axis
# translation", 0, ID["EPSG", 8606]], PARAMETER["Z-axis translation", 0, ID["EPSG", 8607]],
# PARAMETER["X-axis rotation", 0, ID["EPSG", 8608]], PARAMETER["Y-axis rotation", 0,
# ID["EPSG", 8609]], PARAMETER["Z-axis rotation", 0, ID["EPSG", 8610]], PARAMETER["Scale
# difference", 1, ID["EPSG", 8611]]]]

raster::wkt(Ounila_catchment_reprojected)
# "GEOGCRS[\"WGS 84\",\n    DATUM[\"World Geodetic System 1984\",\n
#  ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n
#  LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n
#  ANGLEUNIT[\"degree\",0.0174532925199433]],\n
#  CS[ellipsoidal,2],\n
#  AXIS[\"latitude\",north,\n            ORDER[1],\n
#  ANGLEUNIT[\"degree\",0.0174532925199433]],\n
#  AXIS[\"longitude\",east,\n            ORDER[2],\n
#  ANGLEUNIT[\"degree\",0.0174532925199433]],\n
#  ID[\"EPSG\",4326]]"

# The correct way to reproject the data would be to use a wkt argument, however, projecRaster
# does not accept this
# When looking at the data, it does seem that no errors are made in reprojection
# However, when all functions in the rasterpackage have adapted to the new proj upgrade
# it is definetly safer to reproject the data with wkt description

# https://stackoverflow.com/questions/63727886/proj4-to-proj6-upgrade-and-discarded-datum-warnings


# Inspect stack and drop erronous and empty scenes ------------------------

# various scenes have errors, these are obvious from visual inspection as stripes with deviant values
# In many scenes the erronous data have values below zero (f.i. layer 3, 384, 710, 728, 764, 789, 833 etc.),
# since values below zero will be removed anyway,
# these errors are automatically masked out
# In some bands there are erronous data with values above 0 (also recognizable as clear band with unrealistic values)
# 722 has an error that leads to erronous results with values above 0 - these scenes must be removed
# A regular error occurs in 29 Landsat 5 scenes
# 9 31 332 825 1062
# 828 837 841 1004 1008 1010 1011 1015 1070 1080 1084 1087 1091 1094 1097 1099 1104 1107 1156 1163 1166 1179 1186 1288
# Scene 1290 has very different values from the scenes before and after and contains an obvious sensor error,
# so this scene will be removed as well

# remove the scenes mentioned above
NDVI_dropped <- dropLayer(NDVI,c(722,9, 31, 332, 825, 1062,828, 837, 841, 1004, 1008, 1010, 1011, 1015, 1070, 1080,
                                1084, 1087, 1091, 1094, 1097, 1099, 1104, 1107, 1156, 1163, 1166, 1179, 1186, 1288, 1290))

# Now all layers that contain only NA values are removed

# find layers with only NA:
which(is.na(maxValue(NDVI_dropped)==1))
which(is.na(minValue(NDVI_dropped)==-1))
# output:
# [1]   23   28   37   38   69   87   89  105  115  135  139  143  180  182  184  195  221  222  224  228  232
# [22]  252  277  307  315  382  384  393  399  405  409  442  467  474  481  485  496  522  533  578  598  612
# [43]  627  629  651  671  679  684  715  773  786  793  794  798  811  824  855  867  871  915  924  930  942
# [64]  948  953  961  968  970  984  986 1054 1057 1071 1099 1189 1203 1232 1268 1269 1278 1309 1362 1409 1410
# [85] 1452 1470 1471 1485 1499 1550 1551 1557 1644 1659 1664 1732 1754 1797 1830 1831 1870 1928

# remove those layers from the stack
NDVI_empty_scenes_removed <- dropLayer(NDVI_dropped,which(is.na(maxValue(NDVI_dropped)==1)))

# Set all negative values to NA -------------------------------------------

NDVI_no_neg <- reclassify(NDVI_empty_scenes_removed, cbind(-Inf, 0, NA), right=FALSE)

writeRaster(NDVI_no_neg,"NDVI_no_neg.grd",overwrite=TRUE)

# Apply transformation to Landsat 8 data ----------------------------------

# It would have been easier and more neat to do this earlier (before the stacks were combined)
LC08 <- raster::subset(NDVI_no_neg, grep('LC08', names(NDVI_no_neg), value = T))

LT04_LT05_LT07 <- raster::dropLayer(NDVI_no_neg,grep('LC08', names(NDVI_no_neg)))
LC08_transformed <- 0.0029 + 0.9589 * LC08

NDVI_transformed <- NDVI_no_neg

for(i in 1:nlayers(NDVI_no_neg[[grep('LC08', names(NDVI_no_neg))]])) {
  NDVI_transformed[[grep('LC08', names(NDVI_transformed))[[i]]]] <- LC08_transformed[[i]]
} 

writeRaster(NDVI_transformed, "NDVI_final.grd")

LT04_LT05_LT07_transformed_RMA <- 0.0149 + 1.0035 * LT04_LT05_LT07

NDVI_transformed_RMA <- NDVI_no_neg

for(i in 1:nlayers(NDVI_no_neg[[grep('LC08', names(NDVI_no_neg), invert=TRUE)]])) {
  NDVI_transformed_RMA[[grep('LC08', names(NDVI_transformed_RMA), invert=TRUE)[[i]]]] <- LT04_LT05_LT07_transformed_RMA[[i]]
} 

writeRaster(NDVI_transformed_RMA, "NDVI_final_RMA.grd", overwrite=TRUE)

LT04_LT05_LT07_transformed_OLS <- 0.0235 + 0.9723 * LT04_LT05_LT07

NDVI_transformed_OLS <- NDVI_no_neg

for(i in 1:nlayers(NDVI_no_neg[[grep('LC08', names(NDVI_no_neg), invert=TRUE)]])) {
  NDVI_transformed_OLS[[grep('LC08', names(NDVI_transformed_OLS), invert=TRUE)[[i]]]] <- LT04_LT05_LT07_transformed_OLS[[i]]
} 

writeRaster(NDVI_transformed_OLS, "NDVI_final_OLS.grd", overwrite=TRUE)

# The best transformation is: NDVI_final_OLS 

# Cut NDVI stack in smaller blocks ----------------------------------------

NDVI <- brick("NDVI_final_OLS.grd") 

setwd("G:/Thesis_data/NDVI/sub_areas/tinier")

# based on : https://stackoverflow.com/questions/29784829/r-raster-package-split-image-into-multiples
# I want to divide the grid in 8x8=64 blocks
nblocks=26
nrowblock <- ceiling(nrow(NDVI)/nblocks)
ncolblock <- ceiling(ncol(NDVI)/nblocks)
# create a mold to
NDVI_mold<-NDVI[[1]]
mold <- raster::aggregate(NDVI_mold,fact=c(ncolblock,nrowblock),expand=TRUE)
mold[] <- 1:ncell(mold)
mold_polygon <- rasterToPolygons(mold)
names(mold_polygon) <- "mold"

block_list <- list()
for(i in 1:ncell(mold)){
  e1 <- extent(mold_polygon[mold_polygon$mold==i,])
  block_list[[i]] <- crop(NDVI,e1)
}
for(i in 1:length(block_list)){
  writeRaster(block_list[[i]],filename=paste0("NDVI_splitted_",i,".grd"),overwrite=TRUE)  
}

rm(block_list,splitted_NDVI,mold_polgyon,mold,NDVI_mold,e1)

# Test whether the 
countEMPTY <- 0 
countDATA <- 0 
read_list <- list()
for(i in 1:168){ # 64 pieces
  splitted_NDVI <- brick(paste0("NDVI_splitted_",i,".grd"))
  read_list[[i]] <- splitted_NDVI
  # if (which(!is.na(minValue(splitted_NDVI))) print(i)
  print(i)
  print(all(is.na(minValue(splitted_NDVI))))
  if ((all(is.na(minValue(splitted_NDVI))))==TRUE){
    countEMPTY<-1+countEMPTY
  } else  {
    countDATA<-1+countDATA 
  }
}
print(countEMPTY)
print(countDATA)

# Relocate all the empty blocks to another directory

# check which blocks are empty using:which(!is.na(minValue(read_list[[i]])))
# result: 1,2,31,32,40,41,42,48,49,50,55,56,57,58,59,60,63,64
# These tiles are not uploaded to Cartesius, because they contain no data anyway

# Test whether the blocks can be merged correctly -------------------------

# Merge back the tiles to one rasterbrick to see if result is the same as initial NDVI
# This resulted in a lot of problems, probably due to memory limitations in R
# Rows of NA appeared in between rows of blocks
# After much trial and error I found that when first merging the blocks incrementally in columns
# And then merging the columns together, avoids the problem.

#Manually make a list of blocknrs (that actually contains data, in cols (maybe then the NA doesnt occur))
seq1_9 <- c(17,25,33)
seq2_10 <- c(18,26,34)
seq3_3 <- c(11,19,27,35,43,51)
seq4_4 <- c(12,20,28,36,44,52)
seq5_5 <- c(13,21,29,37,45,53,61)
seq6_6 <- c(14,22,30,38,46,54,62)
seq7_7 <- c(15,23,39,47)
seq8_8 <- c(16,24)

# blocks are merged together per column, using only blocks that contain data
NDVI_mosaicked_1<-read_list[[9]]
for (i in seq1_9) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_1 <- merge(NDVI_mosaicked_1,NDVI_add)
}

writeRaster(NDVI_mosaicked_1,"NDVI_mosaicked_1.grd",overwrite=TRUE)


rm(NDVI_mosaicked_1)

NDVI_mosaicked_2<-read_list[[10]]
for (i in seq2_10) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_2 <- merge(NDVI_mosaicked_2,NDVI_add)
}

writeRaster(NDVI_mosaicked_2,"NDVI_mosaicked_2.grd",overwrite=TRUE)

rm(NDVI_mosaicked_2)

NDVI_mosaicked_3<-read_list[[3]]
for (i in seq3_3) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_3 <- merge(NDVI_mosaicked_3,NDVI_add)
}

writeRaster(NDVI_mosaicked_3,"NDVI_mosaicked_3.grd",overwrite=TRUE)

rm(NDVI_mosaicked_3)

NDVI_mosaicked_4<-read_list[[4]]
for (i in seq4_4) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_4 <- merge(NDVI_mosaicked_4,NDVI_add)
}

writeRaster(NDVI_mosaicked_4,"NDVI_mosaicked_4.grd",overwrite=TRUE)

rm(NDVI_mosaicked_4)

NDVI_mosaicked_5<-read_list[[5]]
for (i in seq5_5) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_5 <- merge(NDVI_mosaicked_5,NDVI_add)
}

writeRaster(NDVI_mosaicked_5,"NDVI_mosaicked_5.grd",overwrite=TRUE)

rm(NDVI_mosaicked_5)

NDVI_mosaicked_6<-read_list[[6]]
for (i in seq6_6) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_6 <- merge(NDVI_mosaicked_6,NDVI_add)
}

writeRaster(NDVI_mosaicked_6,"NDVI_mosaicked_6.grd",overwrite=TRUE)

rm(NDVI_mosaicked_6)

NDVI_mosaicked_7<-read_list[[7]]
for (i in seq7_7) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_7 <- merge(NDVI_mosaicked_7,NDVI_add)
}

writeRaster(NDVI_mosaicked_7,"NDVI_mosaicked_7.grd",overwrite=TRUE)

rm(NDVI_mosaicked_7)

NDVI_mosaicked_8<-read_list[[8]]
for (i in seq8_8) {
  NDVI_add <- read_list[[i]]
  NDVI_mosaicked_8 <- merge(NDVI_mosaicked_8,NDVI_add)
}

writeRaster(NDVI_mosaicked_8,"NDVI_mosaicked_8.grd",overwrite=TRUE)

rm(NDVI_mosaicked_8)
   
NDVI_mosaicked_1 <- brick("NDVI_mosaicked_1.grd")
NDVI_mosaicked_2 <- brick("NDVI_mosaicked_2.grd")
NDVI_mosaicked_3 <- brick("NDVI_mosaicked_3.grd")
NDVI_mosaicked_4 <- brick("NDVI_mosaicked_4.grd")
NDVI_mosaicked_5 <- brick("NDVI_mosaicked_5.grd")
NDVI_mosaicked_6 <- brick("NDVI_mosaicked_6.grd")
NDVI_mosaicked_7 <- brick("NDVI_mosaicked_7.grd")
NDVI_mosaicked_8 <- brick("NDVI_mosaicked_8.grd")

# columns are merged together
NDVI_mosaicked <- merge(NDVI_mosaicked_1,NDVI_mosaicked_2,NDVI_mosaicked_3,NDVI_mosaicked_4,NDVI_mosaicked_5,NDVI_mosaicked_6,NDVI_mosaicked_7,NDVI_mosaicked_8)

# check whether the original NDVI is indeed the same as the mosaicked one
names(NDVI_mosaicked)<-names(NDVI)
compareRaster(NDVI,NDVI_mosaicked,values=TRUE)
# [1] TRUE