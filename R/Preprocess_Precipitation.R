rm(list = ls())  # clear workspace

# Install packages and load libraries -------------------------------------

# install packages if required 
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}
if(!require(ncdf4)){install.packages("ncdf4")}
if(!require(lubridate)){install.packages("lubridate")}

# load libraries
library(rgdal)
library(raster) # package for raster manipulation
library(ncdf4) # package for netcdf manipulation 
library(lubridate) # package for netcdf manipulation 

# Set working directory ---------------------------------------------------

workdirdata <- "G:/Thesis_data/Precipitation"
workdir_GIS_Ounila = "G:/Thesis_data/GIS_Ounila"
workdirNDVI <- "G:/Thesis_data/NDVI"

setwd(workdirdata)

# Create rasterstacks -----------------------------------------------------

ncfiles <- list.files(path = workdirdata, pattern = '.nc$') # list of netcdf files

# one nc file per year, 
# containing 1 variable (total precipitation) and 3 dimensions (long, lat, time)
# time has a size of 12 (each month)
# The total precipitation data is reanalysis data aggregated to monthly means of daily total precipitation

# Create a rasterstack and dataframe with dates stored 
date_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(date_df) <- c("Date", "Year", "Month")

tp_stack <- stack()
for(i in 1:length(ncfiles)) {
  # # create a rasterbrick
  tp_brick <- brick(ncfiles[i] , varname = "tp")
  # Add date to dataframe
  date_df_add <- as.data.frame(as.Date(getZ(tp_brick)))
  colnames(date_df_add) <- "Date"
  date_df_add$Year <- lubridate::year(date_df_add$Date)
  date_df_add$Month <- lubridate::month(date_df_add$Date)
  date_df <- rbind(date_df, date_df_add)
  # create raster stack
  tp_stack <- stack(tp_stack, tp_brick)
}

rm(tp_brick,date_df_add)  # remove unused variables from workspace

saveRDS(date_df,file="date_df.Rda")

# Now save the rasterstacks
writeRaster(tp_stack, "tp_stack.grd", overwrite=TRUE)

# Load the rasterstack - if R session was terminated in the meantime
tp_stack <- stack("tp_stack.grd")

# Reproject data to NDVI resolution ---------------------------------------

setwd(workdirNDVI)

NDVI <- stack("NDVI.grd")

setwd(workdirdata)

tp_reprojected <- projectRaster(tp_stack,NDVI,method-"bilinear",filename="tp_reprojected.grd",overwrite=TRUE)

# Mask to extent of Ounila watershed --------------------------------------

setwd(workdir_GIS_Ounila)

Ounila_catchment_reprojected <- readOGR(dsn=workdir_GIS_Ounila,layer="Ounila_catchment_reprojected")

setwd(workdirdata)

tp <- mask(tp_reprojected, Ounila_catchment_reprojected, filename="tp.grd", overwrite=TRUE)


# Convert from m to mm ----------------------------------------------------

tp <- tp * 1000

writeRaster(tp,"tp.grd",overwrite=TRUE)

rm(tp_stack) # remove unused variables from workspace 
