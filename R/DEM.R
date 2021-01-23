# input: NASA Making Earth System Data Records for Use in Research Environments (MeaSUREs) Digital Elevation Model (DEM) dataset,
# shapefile of Ounila watershed
# output: elevation, slope, aspect figures and rastr for the ounila watershed 

rm(list = ls())  # clear workspace

# Install packages and load libraries  ------------------------------------

# Install requried packages
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
if(!require(extrafont)){install.packages("extrafont")}


# load libraries
library(rgdal)
library(raster)
library(RColorBrewer)
library(extrafont)

# Set working directory ---------------------------------------------------

workdir_GIS_Ounila = "G:/Thesis_data/Shapefiles"
workdir_DEM = "G:/Thesis_data/DEM"
workdir_Fig = "G:/Thesis_data/Figures"

# Reproject shapefile Ounila watershed ------------------------------------

setwd(workdir_GIS_Ounila)

Ounila_catchment <- readOGR(dsn=workdir_GIS_Ounila,layer="nasadem_ounila_catchment")

Ounila_catchment_reprojected<-spTransform(Ounila_catchment,"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

writeOGR(Ounila_catchment_reprojected,dsn=workdir_GIS_Ounila, layer="Ounila_catchment_reprojected", driver="ESRI Shapefile")

# Define colour pallette for plots ----------------------------------------
palette <- c(rgb(132,214,0, max=255),rgb(0,171,68, max=255),rgb(0,104,192, max=255),rgb(108,0,163, max=255),rgb(202,0,156, max=255),rgb(225,85,104, max=255),rgb(225,171,71, max=255),rgb(244,250,0, max=255))
colours <- colorRampPalette(colors=palette)

# Plot elevation data -----------------------------------------------------

setwd(workdir_DEM)

DEM_n31w007 <- raster("n31w007.hgt")
DEM_n31w008 <- raster("n31w008.hgt")
DEM <- merge(DEM_n31w007,DEM_n31w008)
DEM_Ounila <- crop(DEM,Ounila_catchment_reprojected)
DEM_Ounila <- mask(DEM_Ounila,Ounila_catchment_reprojected)
writeRaster(DEM_Ounila,"DEM_OUnila.grd",overwrite=TRUE)
writeRaster(DEM_Ounila,"DEM_Ounila.tif")

setwd(workdir_Fig)

png(filename="DEM.png",width=2031,height=1700,res=300,family="Calibri")

plot(DEM_Ounila, col=colours(100), legend=FALSE, axes=TRUE, asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(side=1,seq(-7.35,-7,0.05),font=1,cex.axis=0.8)
axis(side=2,seq(31.05,31.35,0.05),font=1,cex.axis=0.8)
mtext("Longitude",side=1,line=2.5,font=1,cex=1)
mtext("Latitude",side=2,line=2.5,font=1,cex=1)
mtext("Elevation in the Ounila watershed",side=3,line=0.5,font=2,cex=1.2)

plot(Ounila_catchment_reprojected,add=T)

list_labels <- c(minValue(DEM_Ounila),seq(1300,3500,200),maxValue(DEM_Ounila))
plot(DEM_Ounila, legend.only=TRUE, col=colours(100),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=list_labels,
                    labels=list_labels, 
                    cex.axis=0.8,
                    font=1,
                    tcl=-0.4),
     legend.args=list(text='Elevation (m)', side=4, font=1, line=3, cex=1))

dev.off()

# Plot slope data ---------------------------------------------------------

setwd(workdir_DEM)

slope_n31w007 <- raster("n31w007.slope.hgt")
slope_n31w008 <- raster("n31w008.slope.hgt")
slope <- merge(slope_n31w007,slope_n31w008)
# apply scaling factor
slope <- slope * 0.01
slope_Ounila <- crop(slope,Ounila_catchment_reprojected)
slope_Ounila <- mask(slope_Ounila,Ounila_catchment_reprojected)
writeRaster(slope_Ounila,"slope_OUnila.grd")

setwd(workdir_Fig)

png(filename="slope.png",width=1890,height=1890,res=300)

list_labels <- c(seq(minValue(slope_Ounila),maxValue(slope_Ounila),10),70)

plot(slope_Ounila, col=colours(100), legend=FALSE, axes=TRUE, 
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(side=1,seq(-7.35,-7,0.05),font=1,cex.axis=0.8)
axis(side=2,seq(31.05,31.35,0.05),font=1,cex.axis=0.8)
mtext("Longitude",side=1,line=2.5,font=1,cex=1)
mtext("Latitude",side=2,line=2.5,font=1,cex=1)
mtext("Slope in the Ounila watershed",side=3,line=0.5,font=2,cex=1.2)

plot(Ounila_catchment_reprojected,add=T)

plot(slope_Ounila, legend.only=TRUE, col=colours(100),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=list_labels,
                    labels=list_labels, 
                    cex.axis=0.8,
                    font=1,
                    tcl=-0.4),
     legend.args=list(text='Slope (°)', side=4, font=1, line=3, cex=1))

dev.off()

# Plot aspect data --------------------------------------------------------

setwd(workdir_DEM)
DEM_Ounila<-raster("DEM_OUnila.grd")
# compute aspect from DEM data
aspect_Ounila <- terrain(DEM_Ounila,opt=c("aspect"),unit="degrees")
writeRaster(aspect_Ounila,"aspect_Ounila.grd")

setwd(workdir_Fig)

png(filename="aspect.png",width=1890,family="Calibri",height=1350,res=300)

par(mfrow = c(1,1), mar=c(5,4,4,9)+0.1)

list_labels <- c(minValue(aspect_Ounila),seq(minValue(aspect_Ounila)+22.5,maxValue(aspect_Ounila),45),maxValue(aspect_Ounila))
aspect_colors <- c(rgb(132,214,0, max=255),rgb(0,171,68, max=255),rgb(0,104,192, max=255),rgb(108,0,163, max=255),rgb(202,0,156, max=255),rgb(225,85,104, max=255),rgb(225,171,71, max=255),rgb(244,250,0, max=255),rgb(132,214,0, max=255))
aspect_names = c("North","North-East","East","South-East","South","South-West","West","North-West","North")

plot(aspect_Ounila,col=aspect_colors,breaks=list_labels, legend=FALSE, axes=TRUE, asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(side=1,seq(-7.35,-7,0.05),font=1,cex.axis=0.8)
axis(side=2,seq(31.05,31.35,0.05),font=1,cex.axis=0.8)
mtext("Longitude",side=1,line=2.5,font=1,cex=1)
mtext("Latitude",side=2,line=2.5,font=1,cex=1)
mtext("Aspect in the Ounila watershed",side=3,line=0.5,font=2,cex=1.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

plot(Ounila_catchment_reprojected,add=T)

plot(extent(-7.198685 ,-7.161118 ,31.29425,31.32102), col="black",add=T)
text(x=-7.198685+0.01 ,y=31.32102 ,labels="A",cex=0.6,pos=1,col="black",font=2,offset=0.2)

plot(extent(-7.122735,-7.071827,31.28253,31.3298), col="black",add=T)
text(x=-7.122735+0.01,y=31.3298,labels="B",cex=0.6,pos=1,col="black",font=2,offset=0.2)

plot(extent(-7.221368,-7.18137 ,31.27071,31.29299), col="black",add=T)
text(x=-7.221368+0.01,y=31.29299,labels="C",cex=0.6,pos=1,col="black",font=2,offset=0.2)

plot(extent(-7.275003 ,-7.243186,31.24208,31.26753), col="black",add=T)
text(x=-7.275003+0.01 ,y=31.26753 ,labels="D",cex=0.6,pos=1,col="black",font=2,offset=0.2)

plot(extent(-7.14728 ,-7.083645,31.23253,31.26299), col="black",add=T)
text(x=-7.14728+0.01 ,y=31.26299 ,labels="E",cex=0.6,pos=1,col="black",font=2,offset=0.2)


legend(x='right',
       legend=c("North (337.5°-360°,0°-22.5°)","North-East (22.5°-67.5°)","East (67.5°-112.5°)","South-East (112.5°-157.5°)","South (157.5°-202.5°)","South-West (202.5°-247.5°)","West (247.5°-292.5°)","North-West (292.5°-337.5°)"),
       fill=c(rgb(132,214,0, max=255),rgb(0,171,68, max=255),rgb(0,104,192, max=255),rgb(108,0,163, max=255),rgb(202,0,156, max=255),rgb(225,85,104, max=255),rgb(225,171,71, max=255),rgb(244,250,0, max=255)),
       bty='n', cex=0.8, text.font=1, xpd=NA,inset=-0.60, y.intersp=1.2)

dev.off()





