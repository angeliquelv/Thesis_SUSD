# Input: NDVI stack
# Output: rasterlayers with counts of NA  and nonNA values
# for the whole stack and various subsets and a dataframe with
# NA count per date
# Plots of mean NDVI per yaer + std
# Time-series of mean NDVI

rm(list = ls())  # clear workspace

# Install packages and load libraries -------------------------------------

# install packages if required
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}
if(!require(lubridate)){install.packages("lubridate")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(scales)){install.packages("scales")}

# load libraries
library(rgdal)
library(raster) # package for raster manipulation
library(lubridate) # package to handle dates more easily
library(bfastSpatial)
library(ggplot2)
library(scales)

# Set working directory ---------------------------------------------------

workdir_NDVI <- "G:/Thesis_data/NDVI"
workdir_GIS_Ounila = "G:/Thesis_data/Shapefiles"
workdir_Fig = "G:/Thesis_data/Figures"

setwd(workdir_NDVI)

# Open rasterstack --------------------------------------------------------

NDVI <- stack("NDVI_OLS.grd")

# Open reprojected shapefile ----------------------------------------------

setwd(workdir_GIS_Ounila)

Ounila_catchment_reprojected <- readOGR(dsn=workdir_GIS_Ounila,layer="Ounila_catchment_reprojected")

# Calculate number of  NA and nonNA data ----------------------------------

setwd(workdir_NDVI)

# calculate the number of data points for each cell in the whole rasterstack
nr_nonNA <- sum(!is.na(NDVI))
nr_nonNA <- mask(nr_nonNA,Ounila_catchment_reprojected,filename="nr_nonNA.grd",overwrite=TRUE)

# calculate the number of NAs for each cell in the whole rasterstack
nr_NA <- sum(is.na(NDVI))
nr_NA <- mask(nr_NA,Ounila_catchment_reprojected,filename="nr_NA.grd")

# calculate the number of data points for each cell in the Landsat 4 data
LT04 <- raster::subset(NDVI, grep('LT04', names(NDVI), value = T))
nr_nonNA_LT04 <- sum(!is.na(LT04))
nr_nonNA_LT04 <- mask(nr_nonNA_LT04,Ounila_catchment_reprojected,filename="nr_nonNA_LT04.grd")

# calculate the number of NAs for each cell in the Landsat 4 data
nr_NA_LT04 <- sum(is.na(LT04))
nr_NA_LT04 <- mask(nr_NA_LT04,Ounila_catchment_reprojected,filename="nr_NA_LT04.grd")

# calculate the number of data points for each cell in the Landsat 5 data
LT05 <- raster::subset(NDVI, grep('LT05', names(NDVI), value = T))
nr_nonNA_LT05 <- sum(!is.na(LT05))
nr_nonNA_LT05 <- mask(nr_nonNA_LT05,Ounila_catchment_reprojected,filename="nr_nonNA_LT05.grd")

# calculate the number of NAs for each cell in the Landsat 5 data
nr_NA_LT05 <- sum(is.na(LT05))
nr_NA_LT05 <- mask(nr_NA_LT05,Ounila_catchment_reprojected,filename="nr_NA_LT05.grd")

# calculate the number of data points for each cell in the Landsat 7 data
LE07 <- raster::subset(NDVI, grep('LE07', names(NDVI), value = T))
nr_nonNA_LE07 <- sum(!is.na(LE07))
nr_nonNA_LE07 <- mask(nr_nonNA_LE07,Ounila_catchment_reprojected,filename="nr_nonNA_LE07.grd")

# calculate the number of NAs for each cell in the Landsat 7 data
nr_NA_LE07 <- sum(is.na(LE07))
nr_NA_LE07 <- mask(nr_NA_LE07,Ounila_catchment_reprojected,filename="nr_NA_LE07.grd")

# calculate the number of data points for each cell in the Landsat 8 data
LC08 <- raster::subset(NDVI, grep('LC08', names(NDVI), value = T))
nr_nonNA_LC08 <- sum(!is.na(LC08))
nr_nonNA_LC08 <- mask(nr_nonNA_LC08,Ounila_catchment_reprojected,filename="nr_nonNA_LC08.grd")

# calculate the number of NAs for each cell in the Landsat 8 data
nr_NA_LC08 <- sum(is.na(LC08))
nr_NA_LC08 <- mask(nr_NA_LC08,Ounila_catchment_reprojected,filename="nr_NA_LC08.grd")

# add time dimension to NDVI stack
dates_names <- substr(names(NDVI),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
NDVI <- setZ(NDVI, date)

# calculate the number of data points for each cell for the period 1984/04/19 - 1999/07/10 (Landsat 4 and 5 only)
P1 <- subset(NDVI, which(getZ(NDVI) >= '1984-04-19' & (getZ(NDVI) <= '1999-07-10')))
nr_nonNA_P1 <- sum(!is.na(P1))
nr_nonNA_P1 <- mask(nr_nonNA_P1,Ounila_catchment_reprojected,filename="nr_nonNA_P1.grd")

# calculate the number of NAs for each cell for the period 1984/04/19 - 1999/07/10 (Landsat 4 and 5 only)
nr_NA_P1 <- sum(is.na(P1))
nr_NA_P1 <- mask(nr_NA_P1,Ounila_catchment_reprojected,filename="nr_NA_P1.grd")

# calculate the number of data points for each cell for the period 1999/07/10 - 2011/11/08 (Landsat 5 and 7)
P2 <- subset(NDVI, which(getZ(NDVI) >= '1999-07-10' & (getZ(NDVI) <= '2011-11-08')))
nr_nonNA_P2 <- sum(!is.na(P2))
nr_nonNA_P2 <- mask(nr_nonNA_P2,Ounila_catchment_reprojected,filename="nr_nonNA_P2.grd")

# calculate the number of NAs for each cell for the period 1999/07/10 - 2011/11/08 (Landsat 5 and 7)
nr_NA_P2 <- sum(is.na(P2))
nr_NA_P2 <- mask(nr_NA_P2,Ounila_catchment_reprojected,filename="nr_NA_P2.grd")

# calculate the number of data points for each cell for the period 2011/11/08 - 2013/03/23 (Landsat 7 only)
P3 <- subset(NDVI, which(getZ(NDVI) >= '2011-11-08' & (getZ(NDVI) <= '2013-03-23')))
nr_nonNA_P3 <- sum(!is.na(P3))
nr_nonNA_P3 <- mask(nr_nonNA_P3,Ounila_catchment_reprojected,filename="nr_nonNA_P3.grd")

# calculate the number of NAs for each cell for the period 2011/11/08 - 2013/03/23 (Landsat 7 only)
nr_NA_P3 <- sum(is.na(P3))
nr_NA_P3 <- mask(nr_NA_P3,Ounila_catchment_reprojected,filename="nr_NA_P3.grd")

# calculate the number of data points for each cell for the period 2013/03/23 - 2020/01/01 (Landsat 7 and 8)
P4 <- subset(NDVI, which(getZ(NDVI) >= '2013-03-23' & (getZ(NDVI) <= '2020-01-01')))
nr_nonNA_P4 <- sum(!is.na(P4))
nr_nonNA_P4 <- mask(nr_nonNA_P4,Ounila_catchment_reprojected,filename="nr_nonNA_P4.grd")

# calculate the number of NAs for each cell for the period 2013/03/23 - 2020/01/01 (Landsat 7 and 8)
nr_NA_P4 <- sum(is.na(P4))
nr_NA_P4 <- mask(nr_NA_P4,Ounila_catchment_reprojected,filename="nr_NA_P4.grd")


# Plot figures of NA and nonNA counts -------------------------------------

setwd(workdir_Fig)

# figure of nr of NA over whole rasterstack
png(filename="nr_NA.png",width=1890,height=1890,res=300)

plot(nr_NA)
mtext("nr of NA counts per pixel",side=3,line=0.5,font=2,cex=1.2)

dev.off()

nr_nonNA<-mask(nr_nonNA,NDVI)

# figure of nr of nonNA over whole rasterstack
png(filename="nr_nonNA.png",family="Calibri",width=2031,height=1700,res=300)
plot(nr_nonNA,asp=1,xaxt='n',yaxt='n',col=viridis_pal(direction=-1,option="viridis")(9),breaks=seq(0,1600,200),colNA="white")
title("nr of nonNA counts per pixel",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

# Figure of nr of NA for different Landsat satellites
png(filename="nr_NA_Landsats.png",width=1890,height=1890,res=300)

par(mfrow = c(2,2), mar=c(3,4,4,2)+0.2, cex.main=1.3)

plot(nr_NA_LT04)
mtext("Landsat 4",side=3,line=0.5,font=2,cex=1.2)

plot(nr_NA_LT05)
mtext("Landsat 5",side=3,line=0.5,font=2,cex=1.2)

plot(nr_NA_LE07)
mtext("Landsat 7",side=3,line=0.5,font=2,cex=1.2)

plot(nr_NA_LC08)
mtext("Landsat 8",side=3,line=0.5,font=2,cex=1.2)

dev.off()

# Figure of nr of nonNA for different Landsat satellites
png(filename="nr_nonNA_Landsats.png",width=1890,height=1890,res=300)

par(mfrow = c(2,2), mar=c(3,4,4,2)+0.1, cex.main=1.3)

plot(nr_nonNA_LT04)
mtext("Landsat 4",side=3,line=0.5,font=2,cex=1.2)

plot(nr_nonNA_LT05)
mtext("Landsat 5",side=3,line=0.5,font=2,cex=1.2)

plot(nr_nonNA_LE07)
mtext("Landsat 7",side=3,line=0.5,font=2,cex=1.2)

plot(nr_nonNA_LC08)
mtext("Landsat 8",side=3,line=0.5,font=2,cex=1.2)

dev.off()

# Figure of nr of NA for different periods
png(filename="nr_NA_periods.png",width=1890,height=1890,res=300)

par(mfrow = c(2,2), mar=c(3,4,4,2)+0.1, cex.main=1.3)

plot(nr_NA_P1)
mtext("1984/04/19 - 1999/07/10",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 5, few Landsat 4)",side=3,line=0.25,font=2,cex=0.8)

plot(nr_NA_P2)
mtext("1999/07/10 - 2011/11/08",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 5 and 7)",side=3,line=0.25,font=2,cex=0.8)

plot(nr_NA_P3)
mtext("2011/11/08 - 2013/03/23",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 7)",side=3,line=0.25,font=2,cex=0.8)

plot(nr_NA_P4)
mtext("2013/03/23 - 2020/01/01",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 7 and 8)",side=3,line=0.25,font=2,cex=0.8)

dev.off()

# Figure of nr of nonNA for different periods
png(filename="nr_nonNA_periods.png",width=1890,height=1890,res=300)

par(mfrow = c(2,2), mar=c(3,4,4,2)+0.1, cex.main=1.3)

plot(nr_nonNA_P1)
mtext("1984/04/19 - 1999/07/10",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 5, few Landsat 4)",side=3,line=0.25,font=2,cex=0.8)

plot(nr_nonNA_P2)
mtext("1999/07/10 - 2011/11/08",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 5 and 7)",side=3,line=0.25,font=2,cex=0.8)

plot(nr_nonNA_P3)
mtext("2011/11/08 - 2013/03/23",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 7)",side=3,line=0.25,font=2,cex=0.8)

plot(nr_nonNA_P4)
mtext("2013/03/23 - 2020/01/01",side=3,line=1.25,font=2,cex=1)
mtext("(Landsat 7 and 8)",side=3,line=0.25,font=2,cex=0.8)

dev.off()


# bfastSpatial::countObs() ------------------------------------------------

setwd(workdirNDVI)

NDVI_brick <- brick("NDVI_final.grd")
# add time dimension to NDVI brick
dates_names <- substr(names(NDVI_brick),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
NDVI_brick <- setZ(NDVI_brick, date)

# Use the bfastSpatial function countObs() to check the nonNA values calculated above
obs <- countObs(NDVI_brick)
writeRaster(obs,"obs.grd")
plot(obs)
plot(nr_nonNA)
hist(obs)
hist(nr_nonNA)

# Make a map of the percentage of NA
percNA <- 100 - countObs(NDVI_brick, as.perc=TRUE)
percNA <- mask(percNA,Ounila_catchment_reprojected)

setwd(workdir_Fig)

png(filename="percent_NA.png",width=1890,height=1890,res=300)

plot(percNA)
mtext("Percent NA per pixel",side=3,line=0.5,font=2,cex=1.2)

dev.off()


# Calculate and plot the nr of NA per date --------------------------------

NDVI_mask_to_100 <- mask(NDVI,Ounila_catchment_reprojected,updatevalue=100)

NA_per_date <- data.frame(freq=rep(0,nlayers(NDVI_mask_to_100)),
                          date=rep(0,nlayers(NDVI_mask_to_100)))

for (i in 1:nlayers(NDVI_mask_to_100)) {
  NA_freq <- freq(NDVI_mask_to_100[[i]],value=NA)
  NA_per_date[i, ] = c(NA_freq,date[i])
}

NA_per_date$date<-as.Date(NA_per_date$date)

setwd(workdirNDVI)
saveRDS(NA_per_date,"NA_per_date.Rda")
setwd(workdir_Fig)

ggplot(data=NA_per_date,aes(y=freq,x=date)) +
  geom_line(size=0.5,col="blue") +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y")) +
  labs(title="Count of missing values over time",x="Year",y="NA count") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
ggsave(filename="density_plot_NA_time.png", width = 8, height = 3, dpi=300)

obs_per_date_total <- data.frame(freq=rep(0,nlayers(NDVI)),
                                 date=rep(0,nlayers(NDVI)))
NA_per_date_total <- data.frame(freq=rep(0,nlayers(NDVI)),
                                date=rep(0,nlayers(NDVI)))
for (i in 1:nlayers(NDVI)) {
  NA_freq_total <- freq(NDVI[[i]],value=NA)
  obs_per_date_total[i, ] = c(ncell(NDVI)-NA_freq_total,date[i])
  NA_per_date_total[i, ] = c(NA_freq_total,date[i])
}

obs_per_date_total$date<-as.Date(obs_per_date_total$date)

setwd(workdirNDVI)
saveRDS(obs_per_date_total,"obs_per_date.Rda")
setwd(workdir_Fig)

ggplot(data=obs_per_date_total,aes(y=freq,x=date)) +
  geom_line(size=0.5,col="blue") +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y")) +
  labs(title="Count of observations over time",x="Year",y="observation count") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="density_plot_obs_time.png", width = 8, height = 3, dpi=300)

# Calculate annual mean and std of NDVI per pixel -------------------------

annual_mean_NDVI <- annualSummary(NDVI_brick,fun=mean,na.rm=TRUE)
annual_mean_NDVI <- annual_mean_NDVI[[1:36]]
names(annual_mean_NDVI) <- paste0("Mean NDVI ",1984:2019)
setwd(workdirNDVI)
writeRaster(annual_mean_NDVI,"annual_mean_NDVI.grd")
setwd(workdir_Fig)

png(filename="annual_mean_NDVI_1.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI[[1:12]],cex.main=0.9,zlim=c(0.00,0.75))

dev.off()

png(filename="annual_mean_NDVI_2.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI[[13:24]],cex.main=0.9,zlim=c(0.00,0.75))

dev.off()

png(filename="annual_mean_NDVI_3.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI[[25:36]],cex.main=0.9,zlim=c(0.00,0.75))

dev.off()

annual_std_NDVI <- annualSummary(NDVI_brick,fun=sd,na.rm=TRUE)
annual_std_NDVI <- annual_std_NDVI[[1:36]]
names(annual_std_NDVI) <- paste0("Standard Deviation ",1984:2019)
setwd(workdirNDVI)
writeRaster(annual_std_NDVI,"annual_std_NDVI.grd",overwrite=TRUE)
setwd(workdir_Fig)

png(filename="annual_std_NDVI_1.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI[[1:12]],cex.main=0.9,zlim=c(0.00,0.45))

dev.off()

png(filename="annual_std_NDVI_2.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI[[13:24]],cex.main=0.9,zlim=c(0.00,0.45))

dev.off()

png(filename="annual_std_NDVI_3.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI[[25:36]],cex.main=0.9,zlim=c(0.00,0.45))

dev.off()

# Mean NDVI per date  -----------------------------------------------------

dates_names <- substr(names(NDVI),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
NDVI <- setZ(NDVI, date)

# With the cellStats function the mean per Rasterlayer (date) is calculated
df_mean_NDVI <- data.frame(matrix(ncol = 5, nrow = nlayers(NDVI)))
colnames(df_mean_NDVI) <- c("Mean", "Std","Date","Year","Month")

for(i in seq(1:nlayers(NDVI))) {
  df_mean_NDVI$Mean[[i]] <- cellStats(NDVI[[i]],stat="mean",na.rm=TRUE)
  df_mean_NDVI$Std[[i]] <- cellStats(NDVI[[i]],stat="sd",na.rm=TRUE)
  # Add date information to separate columns
  df_mean_NDVI$Date[i] <- NDVI@z[["time"]][[i]]
  df_mean_NDVI$Year[i] <- year(NDVI@z[["time"]][[i]])
  df_mean_NDVI$Month[i] <- month(NDVI@z[["time"]][[i]])
}

df_mean_NDVI$Satellite <- substr(names(NDVI),1,4)
df_mean_NDVI$Satellite <- factor(df_mean_NDVI$Satellite)

# saveRDS(df_mean_NDVI,file="df_mean_NDVI.Rda")
df_mean_NDVI <- readRDS("df_mean_NDVI.Rda")

df_mean_NDVI$Date <- as.Date(df_mean_NDVI$Date)

# Plot the data from 1999-2019
ggplot(data=df_mean_NDVI,aes(x=Date,
                             y=Mean,
                             color=Satellite,
                             fill=Satellite),na.rm=TRUE) +
  # geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4) +
  geom_line(size=0.2) +
  geom_point(aes(x=Date,y=Mean),alpha=0.4,size=0.8) +
  coord_cartesian(expand=FALSE) +
  ylim(NA,0.2) +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
  labs(title="NDVI between 1984 and 2019",x="Year",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  scale_color_manual(values=c("#ff3399","#35978f","#00CC99","#bf812d"))
  ggsave(filename="mean_NDVI_1984_2019_per_satellite.png", width = 8, height = 3, dpi=300)

  # Plot the data from 1999-2019
  ggplot(data=df_mean_NDVI,aes(x=Date,
                               y=Mean),na.rm=TRUE) +
    geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4,fill="#35978f") +
    geom_line(size=0.2,color="#35978f") +
    coord_cartesian(expand=FALSE) +
    ylim(NA,0.2) +
    scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
    labs(title="NDVI between 1984 and 2019",x="Year",y="NDVI") +
    theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
          axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
          legend.title=element_text(size=11), legend.text=element_text(size=9)) +
    ggsave(filename="mean_NDVI_1984_2019.png", width = 8, height = 3, dpi=300)

NDVI_transformed_RMA <- brick("NDVI_final_RMA.grd")
  
dates_names <- substr(names(NDVI_transformed_RMA),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
NDVI_transformed_RMA <- setZ(NDVI_transformed_RMA, date)
  
# With the cellStats function the mean per Rasterlayer (date) is calculated
df_mean_NDVI <- data.frame(matrix(ncol = 5, nrow = nlayers(NDVI_transformed_RMA)))
colnames(df_mean_NDVI) <- c("Mean", "Std","Date","Year","Month")
  
for(i in seq(1:nlayers(NDVI_transformed_RMA))) {
  df_mean_NDVI$Mean[[i]] <- cellStats(NDVI_transformed_RMA[[i]],stat="mean",na.rm=TRUE)
  df_mean_NDVI$Std[[i]] <- cellStats(NDVI_transformed_RMA[[i]],stat="sd",na.rm=TRUE)
  # Add date information to separate columns
  df_mean_NDVI$Date[i] <- NDVI_transformed_RMA@z[["time"]][[i]]
  df_mean_NDVI$Year[i] <- year(NDVI_transformed_RMA@z[["time"]][[i]])
  df_mean_NDVI$Month[i] <- month(NDVI_transformed_RMA@z[["time"]][[i]])
}
  
df_mean_NDVI$Satellite <- substr(names(NDVI_transformed_RMA),1,4)
df_mean_NDVI$Satellite <- factor(df_mean_NDVI$Satellite)
  
saveRDS(df_mean_NDVI,file="df_mean_NDVI_RMA.Rda")
# df_mean_NDVI <- readRDS("df_mean_NDVI_RMA.Rda")
  
df_mean_NDVI$Date <- as.Date(df_mean_NDVI$Date)
  
# Plot the data from 1999-2019
ggplot(data=df_mean_NDVI,aes(x=Date,
                             y=Mean,
                             color=Satellite,
                             fill=Satellite),na.rm=TRUE) +
  # geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4) +
  geom_line(size=0.2) +
  geom_point(aes(x=Date,y=Mean),alpha=0.4,size=0.8) +
  coord_cartesian(expand=FALSE) +
  ylim(0.025,0.225) +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
  labs(title="NDVI between 1984 and 2019 ETM RMA transformed",x="Year",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  scale_color_manual(values=c("#ff3399","#35978f","#00CC99","#bf812d"))
ggsave(filename="mean_NDVI_1984_2019_per_satellite_ETM_RMA.png", width = 8, height = 3, dpi=300)


NDVI_transformed_OLS <- brick("NDVI_final_OLS.grd")
  
dates_names <- substr(names(NDVI_transformed_OLS),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
NDVI_transformed_OLS <- setZ(NDVI_transformed_OLS, date)
  
# With the cellStats function the mean per Rasterlayer (date) is calculated
df_mean_NDVI <- data.frame(matrix(ncol = 5, nrow = nlayers(NDVI_transformed_OLS)))
colnames(df_mean_NDVI) <- c("Mean", "Std","Date","Year","Month")
  
for(i in seq(1:nlayers(NDVI_transformed_OLS))) {
  df_mean_NDVI$Mean[[i]] <- cellStats(NDVI_transformed_OLS[[i]],stat="mean",na.rm=TRUE)
  df_mean_NDVI$Std[[i]] <- cellStats(NDVI_transformed_OLS[[i]],stat="sd",na.rm=TRUE)
  # Add date information to separate columns
  df_mean_NDVI$Date[i] <- NDVI_transformed_OLS@z[["time"]][[i]]
  df_mean_NDVI$Year[i] <- year(NDVI_transformed_OLS@z[["time"]][[i]])
  df_mean_NDVI$Month[i] <- month(NDVI_transformed_OLS@z[["time"]][[i]])
}
  
df_mean_NDVI$Satellite <- substr(names(NDVI_transformed_OLS),1,4)
df_mean_NDVI$Satellite <- factor(df_mean_NDVI$Satellite)
  
saveRDS(df_mean_NDVI,file="df_mean_NDVI_OLS.Rda")
# df_mean_NDVI <- readRDS("df_mean_NDVI_OLS.Rda")
  
df_mean_NDVI$Date <- as.Date(df_mean_NDVI$Date)

setwd(workdir_Fig)

# Plot the data from 1999-2019
ggplot(data=df_mean_NDVI,aes(x=Date,
                             y=Mean,
                             color=Satellite,
                             fill=Satellite),na.rm=TRUE) +
  # geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4) +
  geom_line(size=0.2) +
  geom_point(aes(x=Date,y=Mean),alpha=0.4,size=0.8) +
  coord_cartesian(expand=FALSE) +
  ylim(0.025,0.225) +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
  labs(title="NDVI between 1984 and 2019 ETM OLS transformed",x="Year",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  scale_color_manual(values=c("#ff3399","#35978f","#00CC99","#bf812d"))
ggsave(filename="mean_NDVI_1984_2019_per_satellite_ETM_OLS.png", family="Calibri",width = 8, height = 3, dpi=300)  


# Density plot of scenes over time ----------------------------------------
dates_names <- substr(names(NDVI_transformed_OLS),5,12)
dates_df<-as.data.frame(dates_names)
dates_df$dates_names <- as.Date(dates_df$dates_names ,format="%Y%m%d",origin="1970-01-01")

setwd(workdir_Fig)
# Density plot
ggplot(dates_df, aes(x=dates_names)) +  geom_density(aes(y = ..count..)) + 
  labs(x="Year",y="Number of scenes") +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(min(dates_df$dates_names),max(dates_df$dates_names))) +
ggsave(filename="scenes_count.png", family="Calibri", width = 8, height = 3, dpi=300)

# Mean NDVI per Landsat 7 and Landsat 8 -----------------------------------

setwd(workdirNDVI)
  
NDVI_brick <- brick("NDVI_final.grd")

LE07_brick <- raster::subset(NDVI_brick, grep('LE07', names(NDVI_brick), value = T))
# add time dimension to NDVI brick
dates_names <- substr(names(LE07_brick),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
LE07_brick <- setZ(LE07_brick, date)

annual_mean_NDVI_LE07 <- annualSummary(LE07_brick,fun=mean,na.rm=TRUE)
names(annual_mean_NDVI_LE07) <- paste0("Mean NDVI ",1999:2019)
setwd(workdirNDVI)
writeRaster(annual_mean_NDVI_LE07,"annual_mean_NDVI_LE07.grd",overwrite=TRUE)
setwd(workdir_Fig)
  
png(filename="annual_mean_NDVI_LE07_1.png",width=1890,height=1260,res=300)
  
plot(annual_mean_NDVI_LE07[[1:12]],cex.main=0.9,zlim=c(0.00,0.75))
  
dev.off()
  
png(filename="annual_mean_NDVI_LE07_2.png",width=1890,height=1260,res=300)
  
plot(annual_mean_NDVI_LE07[[13:21]],cex.main=0.9,zlim=c(0.00,0.75))
  
dev.off()
 
annual_std_NDVI_LE07 <- annualSummary(LE07_brick,fun=sd,na.rm=TRUE)
names(annual_std_NDVI_LE07) <- paste0("Standard Deviation ",1999:2019)
setwd(workdirNDVI)
writeRaster(annual_std_NDVI_LE07,"annual_std_NDVI_LE07.grd",overwrite=TRUE)
setwd(workdir_Fig)
  
png(filename="annual_std_NDVI_LE07_1.png",width=1890,height=1260,res=300)
  
plot(annual_std_NDVI_LE07[[1:12]],cex.main=0.9,zlim=c(0.00,0.45))
  
dev.off()
  
png(filename="annual_std_NDVI_LE07_2.png",width=1890,height=1260,res=300)
  
plot(annual_std_NDVI_LE07[[13:21]],cex.main=0.9,zlim=c(0.00,0.45))
  
dev.off()
  
LC08_brick <- raster::subset(NDVI_brick, grep('LC08', names(NDVI_brick), value = T))
# add time dimension to NDVI brick
dates_names <- substr(names(LC08_brick),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
LC08_brick <- setZ(LC08_brick, date)

annual_mean_NDVI_LC08 <- annualSummary(LC08_brick,fun=mean,na.rm=TRUE)
annual_mean_NDVI_LC08 <- annual_mean_NDVI_LC08[[1:7]]
names(annual_mean_NDVI_LC08) <- paste0("Mean NDVI ",2013:2019)
setwd(workdirNDVI)
writeRaster(annual_mean_NDVI_LC08,"annual_mean_NDVI_LC08.grd",overwrite=TRUE)
setwd(workdir_Fig)

png(filename="annual_mean_NDVI_LC08_1.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI_LC08[[1:7]],cex.main=0.9,zlim=c(0.00,0.75))

dev.off()

annual_std_NDVI_LC08 <- annualSummary(LC08_brick,fun=sd,na.rm=TRUE)
annual_std_NDVI_LC08 <- annual_std_NDVI_LC08[[1:7]]
names(annual_std_NDVI_LC08) <- paste0("Standard Deviation ",2013:2019)
setwd(workdirNDVI)
writeRaster(annual_std_NDVI_LC08,"annual_std_NDVI_LC08.grd",overwrite=TRUE)
setwd(workdir_Fig)

png(filename="annual_std_NDVI_LC08_1.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI_LC08[[1:7]],cex.main=0.9,zlim=c(0.00,0.45))

dev.off()


  

