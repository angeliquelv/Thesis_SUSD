# Input files for this scripts are:
# A tif-file with climate classification data from Beck et al. (2018) (Beck_KG_V1_present_0p0083)
# hdf-files with MODIS Land Cover classification data from 2001 and 2018, obtained via https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/

rm(list = ls())  # clear workspace

# Install packages and load libraries  ------------------------------------

# Install requried packages
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}
if(!require(MODIS)){install.packages("MODIS")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(ggmap)){install.packages("ggmap")}
if(!require(gdalUtils)){install.packages("gdalUtils")}
if(!require(broom)){install.packages("broom")}


# load libraries
library(rgdal)
library(raster)
library(MODIS)
library(ggplot2)
library(ggmap)
library(gdalUtils)
library(broom)

# Set working directory ---------------------------------------------------

workdir = "C:/Users/Boris/Documents/MSc Sustainable Development/Thesis - Morrocan Project/R-scripts"
workdir_GIS_Ounila = "G:/Thesis_data/Shapefiles"
workdir_LandUse = "G:/Thesis_data/Land_Use"
workdir_ClimateClassification = "G:/Thesis_data/Climate_Classification"
workdir_Fig = "G:/Thesis_data/Figures"

# Open reprojected shapefile ----------------------------------------------

setwd(workdir_GIS_Ounila)

Ounila_catchment_reprojected <- readOGR(dsn=workdir_GIS_Ounila,layer="Ounila_catchment_reprojected")

# Map of Morocco ----------------------------------------------------------

map_Morocco <- get_map("Morocco",
                       zoom=7,
                       size=c(640,640),
                       scale=2,
                       format="png8",
                       maptype="terrain",
                       language = "en-EN",
                       filename="terrain_Morocco.png",
                       color="color")

setwd(workdir_Fig)

png(filename="Morocco_terrain_map_1.png", width=1200, height=1200, res=300)

Ounila_catchment_fortified <- tidy(Ounila_catchment_reprojected)

ggmap(map_Morocco) +
       ggplot2::geom_polygon(data=Ounila_catchment_fortified,aes(x=long,y=lat,group=group),fill="transparent",col="red",lwd=0.5)

dev.off()

setwd(workdir_Fig)
map_Morocco_2 <- get_map("Morocco",
                         zoom=6,
                         size=c(640,640),
                         scale=2,
                         format="png8",
                         maptype="terrain",
                         language = "en-EN",
                         filename="terrain_Morocco_2.png",
                         color="color")

png(filename="Morocco_terrain_map_2.png", width=1200, height=1200, res=300)

ggmap(map_Morocco_2) +
        ggplot2::geom_polygon(data=Ounila_catchment_fortified,aes(x=long,y=lat,group=group),fill="transparent",col="red",lwd=0.25)
dev.off()

# Satellite image of Ounila valley ----------------------------------------

map_Ounila_3 <- get_googlemap(center=c(-7.151778,31.276726),
                              zoom = 15,
                              size=c(640,640),
                              scale=2,
                              format="png8",
                              maptype="satellite",
                              language = "en-EN",
                              filename="satellite_Ounila_3.png",
                              color="color")

png(filename="Ounila_satellite_map.png", width=1200, height=1200, res=300)

ggmap(map_Ounila_3) +
        ggplot2::geom_polygon(data=Ounila_catchment_fortified,aes(x=long,y=lat,group=group),fill="transparent",col="red",lwd=1) +
        geom_rect(aes(xmin=-7.159,xmax=-7.1565,ymin=31.276,ymax=31.279),alpha=0.05,fill="blue") +
        geom_point(x=-7.151778,y=31.276726, color="red",size=1.5) +
        geom_text(x=-7.151778,y=31.276726,aes(label="Anguelz"),hjust=-0.2,color="white",size=3)

dev.off()

# # Create Land Cover Classification  maps ----------------------------------
#
# # Bounding box of the Ounila valley, where Anguelz is located
# Ounila_extent <- extent(-7.275, -6.981, 31.017, 31.360)
# getTile(Ounila_extent)
# # Output:
# # Slot "tile":
# #         [1] "h17v05"
#
# # Create list with names of hdf-files in working directory
# hdf_files <- list.files(path = workdir, pattern = '.hdf$')
#
# # Now manually download these data tiles from:
# # https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/
# # for the year 2001 and 2018
#
# # check what subsdatasets are in the hdf files
# sds_2001 <-get_subdatasets(hdf_files[1])
# sds_2001
# sds_2018 <-get_subdatasets(hdf_files[2])
# sds_2018
# # I am only interested in the first one: Land Cover Type 1
#
# names <- paste0("LandCover",c("2001","2018"),".tif")
#
# # Now convert these files from hdf to tif format
# for (i in 1:length(hdf_files)){
#         sds <- get_subdatasets(hdf_files[i])
#         gdal_translate(sds[1], dst_dataset=names[i])
# }
#
# names_qc <- paste0("LandCover_qc_",c("2001","2018"),".tif")
#
# # Now convert these files from hdf to tif format
# for (i in 1:length(hdf_files)){
#         sds <- get_subdatasets(hdf_files[i])
#         gdal_translate(sds[12], dst_dataset=names_qc[i])
# }
#
# # Create list with names of tif-files in working directory
# tif_files <- list.files(path = workdir, pattern = '.tif$')
#
# #  Make a rasterstack of the Land Cover data from 2001 and 2018
# LC_stack <- stack(tif_files[4:5])
# QC_stack <- stack(tif_files[2:3])
#
# names(LC_stack) <- c("Land Cover Ounila valley 2001", "Land Cover Ounila valley 2018")
#
# # convert from sinusoidal to longitude latitude
# llcrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
# LC_stack <- projectRaster(LC_stack, crs=llcrs, method = "ngb")
# QC_stack <- projectRaster(QC_stack, crs=llcrs, method = "ngb")
#
# # crop to region of interest (because of sinusoidal, more data is shown than wanted)
# LC_stack <- crop(LC_stack, Ounila_extent)
# QC_stack <- crop(QC_stack, Ounila_extent)
#
# writeRaster(LC_stack, "LC_stack_Ounila.grd", overwrite=TRUE)
# writeRaster(QC_stack, "QC_stack_Ounila.grd", overwrite=TRUE)
#
# # Test whether there are pixels with a quality flag attached to it
# summary(QC_stack)
# # Every pixel has value zero, which means that it the pixels all have been classified as land
# # So there is no need to mask out certain values
#
# # # In case R was terminated in the mean time: Open raster file again
# # LC_stack <- stack("LC_stack_Ounila.grd")
# # QC_stack <- stack("QC_stack_Ounila.grd")
#
# # Define color palette for Land Cover Classifcication
# igbpPalette = c('#05450a', '#086a10', '#54a708', '#78d203', '#009900', '#c6b044', '#dcd159',
#                 '#dade48', '#fbff13', '#b6ff05', '#27ff87', '#c24f44', '#a5a5a5', '#ff6d4c',
#                 '#69fff8', '#f9ffa4', '#1c0dff')
#
#
# # Create Figure containing image of both the Land Cover map of 2000 and 2018
# png(filename="LandCoverOunila.png", width=2400, height=1000, res=300)
#
# par(mfrow = c(1,3), mar=c(5,4,4,2)+0.1, cex.main=1.3)
#
# plot(LC_stack[[1]], legend = FALSE, col = igbpPalette, zlim=c(1,17), main="Land Cover Ounila valley 2001",
#      xlab="Longitude",ylab="Latitude")
#
# points(x=-7.151778,y=31.276726,type="p",pch=19,col="black",cex=0.65)
# text(x=-7.151778,y=31.276726,labels="Anguelz",cex=0.8,adj=-0.2)
#
# plot(LC_stack[[2]], legend = FALSE, col = igbpPalette, zlim=c(1,17), main="Land Cover Ounila valley 2018",
#      xlab="Longitude",ylab="Latitude")
#
# points(x=-7.151778,y=31.276726,type="p",pch=19,col="black",cex=0.65)
# text(x=-7.151778,y=31.276726,labels="Anguelz",cex=0.8,adj=-0.2)
#
# igbpPaletteVis = c('#dcd159',
#                    '#b6ff05',
#                    '#f9ffa4')
#
# legend(x='right', legend =  c("Open Shrubland",
#                               "Grassland",
#                               "Barren"
#                               ),
#        fill = igbpPaletteVis,
#        bty='n', cex=1, xpd=NA,inset=-0.85)
#
# dev.off()
#
# Create detailed Land Use map --------------------------------------------

setwd(workdir_LandUse)

tif_file <- list.files(path = workdir_LandUse, pattern = '.tif$')

Land_Use <- raster(tif_file[1])

Land_Use_Ounila <- crop(Land_Use,Ounila_catchment_reprojected)
Land_Use_Ounila <- mask(Land_Use_Ounila,Ounila_catchment_reprojected)

# Save the cropped and masked Land Use
writeRaster(Land_Use_Ounila,"Land_Use_Ounila.tif",overwrite=TRUE)

# Create Figure of Land Use
png(filename="LandUseOunila.png", width=1890, height=1260, res=300)

par(mfrow = c(1,1), mar=c(5,4,4,12)+0.1)

Land_Use_Palette <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                      rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255),
                      rgb(255,255,255, max=255),rgb(0,70,200, max=255))

plot(Land_Use_Ounila, legend = FALSE, col = Land_Use_Palette, zlim=c(1,10),
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")

axis(side=1,seq(-7.35,-7,0.05),font=1,cex.axis=0.8)
axis(side=2,seq(31.05,31.35,0.05),font=1,cex.axis=0.8)
mtext("Longitude",side=1,line=2.5,font=1,cex=1)
mtext("Latitude",side=2,line=2.5,font=1,cex=1)
mtext("Land Use in the Ounila watershed",side=3,line=0.5,font=2,cex=1.2)

Land_Use_Palette_Vis  <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                           rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                           rgb(255,235,175, max=255),
                           rgb(255,245,215, max=255),rgb(195,20,0, max=255),
                           rgb(0,70,200, max=255))
legend(x='right', legend =  c("Tree cover areas","Shrub cover areas",
                              "Grassland","Cropland",
                              "Lichen Mosses / Sparse Vegetation",
                              "Bare areas", "Built up areas",
                              "Open water"
                              ),
       fill = Land_Use_Palette_Vis,
       bty='n', cex=0.8, xpd=NA,inset=-0.85)

dev.off()

# Create climate classification map of Ounila valley ----------------------

setwd(workdir_ClimateClassification)

# Load climate classification data
KG_Climate <- raster("Beck_KG_V1_present_0p0083.tif")

KG_Climate_Ounila <- crop(KG_Climate,Ounila_catchment_reprojected,snap="out")

writeRaster(KG_Climate_Ounila,"KG_Climate_Ounila.grd",overwrite=TRUE)

# In case R was terminated in the mean time: Open raster file again
KG_Climate_Ounila<- raster("KG_Climate_Ounila.grd")

setwd(workdir_Fig)

# Plot KG Climate Classification
png(filename="ClimateClassificationOunilaWatershed.png", width=1890, height=945, res=300)
par(mfrow = c(1,1), mar=c(4,1,4,16)+0.1)

plot(KG_Climate_Ounila)
plot(Ounila_catchment_reprojected,add=T)

# Define color pallete for climate classification China
KGPalette = c(rgb(245,165,0,max=255), '#ffdc64', '#ffff00', '#c8c800', '#c800c8')

legend(x='right', legend = c("Arid, steppe, hot (Bsh)",
                             "Arid, steppe, cold (Bsk)",
                             "Temperate, dry summer, hot summer (Csa)",
                             "Temperate, dry summer, warm summer (Csb)",
                             "Cold, dry summer, warm summer (Dsb)"),
       fill = KGPalette, bty='n', cex=0.6, xpd=NA,inset=-0.5)

dev.off()
