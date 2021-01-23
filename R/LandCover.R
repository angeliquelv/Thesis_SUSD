# Input files for this scripts are:
# A tif-file with climate classification data from Beck et al. (2018) (Beck_KG_V1_present_0p0083)
# hdf-files with MODIS Land Cover classification data from 2001 and 2018, obtained via https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/
# Detailed Land use map of Africa 20 by 20 m obtained via: http://2016africalandcover20m.esrin.esa.int/download.php
# Google maps and satellite images via Google Cloud Platform
# Output: maps of land cover, climate classification, and satellite image of watershed

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
if(!require(extrafont)){install.packages("extrafont")}


# load libraries
library(rgdal)
library(raster)
library(MODIS)
library(ggplot2)
library(ggmap)
library(gdalUtils)
library(broom)
library(extrafont)

# Set working directory ---------------------------------------------------

workdir = "C:/Users/Boris/Documents/MSc Sustainable Development/Thesis_Morrocan_Project/R_scripts"
workdir_GIS_Ounila = "G:/Thesis_data/Shapefiles"
workdir_LandUse = "G:/Thesis_data/Land_Use"
workdir_ClimateClassification = "G:/Thesis_data/Climate_Classification"
workdir_Fig = "G:/Thesis_data/Figures"

# Open reprojected shapefile ----------------------------------------------

setwd(workdir_GIS_Ounila)

Ounila_catchment_reprojected <- readOGR(dsn=workdir_GIS_Ounila,layer="Ounila_catchment_reprojected")
Ounila_catchment_fortified <- tidy(Ounila_catchment_reprojected)

channels <- readOGR(dsn=workdir_GIS_Ounila,layer="nasadem_channels")
channels_reprojected<-spTransform(channels,"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
channels_reprojected<-crop(channels_reprojected,Ounila_catchment_reprojected)


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

png(filename="Morocco_terrain_map_1.png", family="Calibri",width=1200, height=1200, res=300)


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

map_Ounila_3 <- get_googlemap(center=c(-7.185778,31.210726),
                              zoom = 11,
                              size=c(640,640),
                              scale=2,
                              format="png8",
                              maptype="satellite",
                              language = "en-EN",
                              filename="satellite_Ounila_3.png",
                              color="color")

png(filename="Ounila_satellite_map.png",family="Calibri",width=1200, height=1200, res=300)

ggmap(map_Ounila_3) +
        ggplot2::geom_polygon(data=Ounila_catchment_fortified,aes(x=long,y=lat,group=group),fill="transparent",col="red",lwd=0.8) +
        ggplot2::geom_path(data=channels_reprojected,aes(x=long,y=lat,group=group),col="cyan",lwd=0.3) +
        coord_fixed(ratio=1) +
        geom_rect(aes(xmin=-7.159,xmax=-7.1565,ymin=31.276,ymax=31.279),alpha=0.05,fill="blue") +
        geom_point(x=-7.151778,y=31.276726, color="white",size=1) +
        geom_text(x=-7.151778,y=31.276726,aes(label="Anguelz"),hjust=-0.2,color="white",size=2.5,fontface="bold") +
        geom_point(x=-7.123735,y=31.303808, color="white",size=1) +
        geom_text(x=-7.123735,y=31.303808,aes(label="Tighza"),hjust=-0.2,color="white",size=2.5,fontface="bold") +
        geom_point(x=-7.240218,y=31.288580, color="white",size=1) +
        geom_text(x=-7.240218,y=31.288580,aes(label="Telouet"),hjust=-0.2,color="white",size=2.5,fontface="bold") +
        geom_point(x=-7.132104,y=31.047817, color="white",size=1) +
        geom_text(x=-7.132104,y=31.047817,aes(label="Aït-Ben-Haddou"),hjust=-0.2,color="white",size=2.5,fontface="bold") +
        geom_text(x=-7.150054,y=31.143358,angle=90,aes(label="Asif Ounila"),hjust=-0.2,color="cyan",size=2.5) +
        geom_text(x=-7.250561,y=31.1970677,angle=-40,aes(label="Asif Ounila"),hjust=-0.2,color="cyan",size=2.5) +
        geom_text(x=-7.291051,y=31.255168,angle=-65,aes(label="Asif Mellah"),hjust=-0.2,color="cyan",size=2.5) +
        ggsn::north(x.min = -7.36, x.max = -7.35,
                    y.min = 31.37, y.max = 31.38, scale = 4,
                    symbol=3) + 
        ggsn::scalebar(data=Ounila_catchment_fortified,
                       location="bottomleft", dist=5,
                       dist_unit ="km",
                       transform = TRUE, 
                       model = "WGS84", height = 0.02, 
                       st.bottom=FALSE,st.size=0.8,box.color="transparent",
                       ) +
        geom_text(x=-7.395,y=31.050,aes(label="0"),color="black",size=2) +
        geom_text(x=-7.30,y=31.050,aes(label="10 km"),color="black",size=2) 
        
  
dev.off()

# Create Land Cover Classification  maps ----------------------------------
setwd(workdir)
# Bounding box of the Ounila valley, where Anguelz is located
Ounila_extent <- extent(-7.275, -6.981, 31.017, 31.360)
getTile(Ounila_extent)
# Output:
# Slot "tile":
#         [1] "h17v05"

# Create list with names of hdf-files in working directory
hdf_files <- list.files(path = workdir, pattern = '.hdf$')

# Now manually download these data tiles from:
# https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/
# for the year 2001 and 2018

# check what subsdatasets are in the hdf files
sds_2001 <-get_subdatasets(hdf_files[1])
sds_2001
sds_2018 <-get_subdatasets(hdf_files[2])
sds_2018
# I am only interested in the first one: Land Cover Type 1

names <- paste0("LandCover",c("2001","2018"),".tif")

# Now convert these files from hdf to tif format
for (i in 1:length(hdf_files)){
        sds <- get_subdatasets(hdf_files[i])
        gdal_translate(sds[1], dst_dataset=names[i])
}

names_qc <- paste0("LandCover_qc_",c("2001","2018"),".tif")

# Now convert these files from hdf to tif format
for (i in 1:length(hdf_files)){
        sds <- get_subdatasets(hdf_files[i])
        gdal_translate(sds[12], dst_dataset=names_qc[i])
}

# Create list with names of tif-files in working directory
tif_files <- list.files(path = workdir, pattern = '.tif$')

#  Make a rasterstack of the Land Cover data from 2001 and 2018
LC_stack <- stack(tif_files[4:5])
QC_stack <- stack(tif_files[2:3])

names(LC_stack) <- c("Land Cover Ounila valley 2001", "Land Cover Ounila valley 2018")

# convert from sinusoidal to longitude latitude
llcrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
LC_stack <- projectRaster(LC_stack, crs=llcrs, method = "ngb")
QC_stack <- projectRaster(QC_stack, crs=llcrs, method = "ngb")

# crop to region of interest (because of sinusoidal, more data is shown than wanted)
LC_stack <- crop(LC_stack, Ounila_catchment_reprojected)
QC_stack <- crop(QC_stack, Ounila_catchment_reprojected)

writeRaster(LC_stack, "LC_stack_Ounila.grd", overwrite=TRUE)
writeRaster(QC_stack, "QC_stack_Ounila.grd", overwrite=TRUE)

# Test whether there are pixels with a quality flag attached to it
summary(QC_stack)
# Every pixel has value zero, which means that it the pixels all have been classified as land
# So there is no need to mask out certain values

# # In case R was terminated in the mean time: Open raster file again
# LC_stack <- stack("LC_stack_Ounila.grd")
# QC_stack <- stack("QC_stack_Ounila.grd")

# Define color palette for Land Cover Classifcication
igbpPalette = c('#05450a', '#086a10', '#54a708', '#78d203', '#009900', '#c6b044', '#dcd159',
                '#dade48', '#fbff13', '#b6ff05', '#27ff87', '#c24f44', '#a5a5a5', '#ff6d4c',
                '#69fff8', '#f9ffa4', '#1c0dff')

setwd(workdir_Fig)

# Create Figure containing image of both the Land Cover map of 2000 and 2018
png(filename="LandCoverOunila.png",family="Calibri",width=2031, height=1000, res=300)

par(mfrow = c(1,1), mar=c(3,3,2,13)+0.2)

plot(LC_stack[[1]], legend = FALSE, col = igbpPalette, zlim=c(1,17), asp=1, main="Land Cover Ounila valley 2001",
     xlab="Longitude",ylab="Latitude")
plot(Ounila_catchment_reprojected,add=T)
points(x=-7.151778,y=31.276726,type="p",pch=19,col="black",cex=0.65)
text(x=-7.151778,y=31.276726,labels="Anguelz",cex=0.8,adj=-0.2)

plot(LC_stack[[2]], legend = FALSE, col = igbpPalette, zlim=c(1,17), asp=1, main="Land Cover Ounila valley 2018",
     xlab="Longitude",ylab="Latitude")
plot(Ounila_catchment_reprojected,add=T)
points(x=-7.151778,y=31.276726,type="p",pch=19,col="black",cex=0.65)
text(x=-7.151778,y=31.276726,labels="Anguelz",cex=0.8,adj=-0.2)

igbpPaletteVis = c('#dcd159',
                   '#b6ff05',
                   '#f9ffa4')

legend(x='right', legend =  c("Open Shrubland",
                              "Grassland",
                              "Barren"
                              ),
       fill = igbpPaletteVis,
       bty='n', cex=1, xpd=NA,inset=-0.85)

dev.off()

# Create detailed Land Use map --------------------------------------------

setwd(workdir_LandUse)

tif_file <- list.files(path = workdir_LandUse, pattern = '.tif$')

Land_Use <- raster(tif_file[1])

Land_Use_Ounila <- crop(Land_Use,Ounila_catchment_reprojected)
Land_Use_Ounila <- mask(Land_Use_Ounila,Ounila_catchment_reprojected)

# Save the cropped and masked Land Use
writeRaster(Land_Use_Ounila,"Land_Use_Ounila.tif",overwrite=TRUE)

setwd(workdir_Fig)


sum(Land_Use_Ounila[]==3,na.rm=TRUE)/sum(!is.na(Land_Use_Ounila[]))*100 # grassland
# [1] 34.78632
sum(Land_Use_Ounila[]==4,na.rm=TRUE)/sum(!is.na(Land_Use_Ounila[]))*100 # cropland
# [1] 2.510798
sum(Land_Use_Ounila[]==6,na.rm=TRUE)/sum(!is.na(Land_Use_Ounila[]))*100 # sparse vegetation
# [1] 19.69998
sum(Land_Use_Ounila[]==7,na.rm=TRUE)/sum(!is.na(Land_Use_Ounila[]))*100 # bare
# [1] 42.48259

# Create Figure of Land Use
png(filename="LandUseOunila.png", family="Calibri",width=2031, height=1270, res=300)
par(mfrow = c(1,1), mar=c(3,3,2,9)+0.2)

Land_Use_Palette <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                      rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255),
                      rgb(255,255,255, max=255),rgb(0,70,200, max=255))

plot(Land_Use_Ounila, legend = FALSE, col = Land_Use_Palette, zlim=c(1,10), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("Land Use in the Ounila watershed",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
# text(x=-7.395,y=31.050,labels="0",col="black",cex=0.8) +
# text(x=-7.30,y=31.050,labels="10 km",col="black",cex=0.8) 

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
       bty='n', cex=0.8, xpd=NA,inset=-0.58)

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
png(filename="ClimateClassificationOunilaWatershed.png", family="Calibri",width=2400, height=700, res=300)
par(mfrow = c(1,1), mar=c(0,0,0,30)+0.1)

plot(KG_Climate_Ounila,asp=1)
plot(Ounila_catchment_reprojected,add=T)

# Define color pallete for climate classification China
KGPalette = c(rgb(245,165,0,max=255), '#ffdc64', '#ffff00', '#c8c800', '#c800c8')

legend(x='right', legend = c("Arid, steppe, hot (Bsh)",
                             "Arid, steppe, cold (Bsk)",
                             "Temperate, dry summer, hot summer (Csa)",
                             "Temperate, dry summer, warm summer (Csb)",
                             "Cold, dry summer, warm summer (Dsb)"),
       fill = KGPalette, bty='n', cex=0.7, xpd=NA,inset=-0.005)

dev.off()
