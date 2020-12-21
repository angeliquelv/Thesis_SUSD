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
library(ggplot2)
library(scales)

# Set working directory ---------------------------------------------------

workdirNDVI <- "G:/Thesis_data/NDVI"
workdir_GIS_Ounila = "G:/Thesis_data/Shapefiles"
workdir_Fig = "G:/Thesis_data/Figures/PermaAtlas"

setwd(workdirNDVI)

NDVI <- brick("NDVI_final.grd")
NDVI <- brick("NDVI_no_neg.grd")

NDVI_PermaAtlas <- crop(NDVI,extent(-7.159,-7.1565,31.276,31.279))
# add time dimension to NDVI brick
dates_names <- substr(names(NDVI_PermaAtlas),5,12)
date <- as.Date(dates_names,format="%Y%m%d")
NDVI_PermaAtlas <- setZ(NDVI_PermaAtlas, date)

# Calculate annual mean and std of NDVI per pixel -------------------------

annual_mean_NDVI <- annualSummary(NDVI_PermaAtlas,fun=mean,na.rm=TRUE)
annual_mean_NDVI <- annual_mean_NDVI[[1:36]]
names(annual_mean_NDVI) <- paste0("Mean NDVI ",1984:2019)
setwd(workdirNDVI)
writeRaster(annual_mean_NDVI,"annual_mean_NDVI_PermaAtlas.grd")
setwd(workdir_Fig)

png(filename="annual_mean_NDVI_1.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI[[1:12]],cex.main=0.9,zlim=c(0.00,0.25))

dev.off()

png(filename="annual_mean_NDVI_2.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI[[13:24]],cex.main=0.9,zlim=c(0.00,0.25))

dev.off()

png(filename="annual_mean_NDVI_3.png",width=1890,height=1260,res=300)

plot(annual_mean_NDVI[[25:36]],cex.main=0.9,zlim=c(0.00,0.25))

dev.off()

annual_std_NDVI <- annualSummary(NDVI_PermaAtlas,fun=sd,na.rm=TRUE)
annual_std_NDVI <- annual_std_NDVI[[1:36]]
names(annual_std_NDVI) <- paste0("Standard Deviation ",1984:2019)
setwd(workdirNDVI)
writeRaster(annual_std_NDVI,"annual_std_NDVI_PermaAtlas.grd",overwrite=TRUE)
setwd(workdir_Fig)

png(filename="annual_std_NDVI_1.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI[[1:12]],cex.main=0.9,zlim=c(0.00,0.20))

dev.off()

png(filename="annual_std_NDVI_2.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI[[13:24]],cex.main=0.9,zlim=c(0.00,0.20))

dev.off()

png(filename="annual_std_NDVI_3.png",width=1890,height=1260,res=300)

plot(annual_std_NDVI[[25:36]],cex.main=0.9,zlim=c(0.00,0.20))

dev.off


# Mean NDVI per date  -----------------------------------------------------

# With the cellStats function the mean per Rasterlayer (date) is calculated
df_mean_NDVI <- data.frame(matrix(ncol = 5, nrow = nlayers(NDVI_PermaAtlas)))
colnames(df_mean_NDVI) <- c("Mean", "Std","Date","Year","Month")

for(i in seq(1:nlayers(NDVI_PermaAtlas))) {
  df_mean_NDVI$Mean[[i]] <- cellStats(NDVI_PermaAtlas[[i]],stat="mean",na.rm=TRUE)
  df_mean_NDVI$Std[[i]] <- cellStats(NDVI_PermaAtlas[[i]],stat="sd",na.rm=TRUE)
  # Add date information to separate columns
  df_mean_NDVI$Date[i] <- NDVI_PermaAtlas@z[["time"]][[i]]
  df_mean_NDVI$Year[i] <- year(NDVI_PermaAtlas@z[["time"]][[i]])
  df_mean_NDVI$Month[i] <- month(NDVI_PermaAtlas@z[["time"]][[i]])
}

# saveRDS(df_mean_NDVI,file="df_mean_NDVI.Rda")

df_mean_NDVI$Date <- as.Date(df_mean_NDVI$Date)

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

df_mean_1984 <- subset(df_mean_NDVI,Year==1984,select=c("Mean","Std","Date","Month"))
# Plot the data from 1984
ggplot(data=df_mean_1984,aes(x=as.Date(Date,format="%Y-%m-%d"),
                             y=Mean),na.rm=TRUE) +
  geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4,fill="#35978f") +
  geom_line(size=0.4,color="#35978f") +
  #coord_cartesian(expand=FALSE) +
  ylim(NA,0.2) +
  scale_x_date(breaks=date_breaks("1 month"),labels=date_format("%b")) +
  labs(title="NDVI in 1984",x="Month",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
ggsave(filename="mean_NDVI_1984.png", width = 8, height = 3, dpi=300)

df_mean_1994 <- subset(df_mean_NDVI,Year==1994,select=c("Mean","Std","Date","Month"))
# Plot the data from 1994
ggplot(data=df_mean_1994,aes(x=as.Date(Date,format="%Y-%m-%d"),
                             y=Mean),na.rm=TRUE) +
  geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4,fill="#35978f") +
  geom_line(size=0.4,color="#35978f") +
  #coord_cartesian(expand=FALSE) +
  ylim(NA,0.2) +
  scale_x_date(breaks=date_breaks("1 month"),labels=date_format("%b")) +
  labs(title="NDVI in 1994",x="Month",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="mean_NDVI_1994.png", width = 8, height = 3, dpi=300)

df_mean_2002 <- subset(df_mean_NDVI,Year==2002,select=c("Mean","Std","Date","Month"))
# Plot the data from 2002
ggplot(data=df_mean_2002,aes(x=as.Date(Date,format="%Y-%m-%d"),
                             y=Mean),na.rm=TRUE) +
  geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4,fill="#35978f") +
  geom_line(size=0.4,color="#35978f") +
  #coord_cartesian(expand=FALSE) +
  ylim(NA,0.2) +
  scale_x_date(breaks=date_breaks("1 month"),labels=date_format("%b")) +
  labs(title="NDVI in 2002",x="Month",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="mean_NDVI_2002.png", width = 8, height = 3, dpi=300)

df_mean_2014 <- subset(df_mean_NDVI,Year==2014,select=c("Mean","Std","Date","Month"))
# Plot the data from 2014
ggplot(data=df_mean_2014,aes(x=as.Date(Date,format="%Y-%m-%d"),
                             y=Mean),na.rm=TRUE) +
  geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4,fill="#35978f") +
  geom_line(size=0.4,color="#35978f") +
  #coord_cartesian(expand=FALSE) +
  ylim(NA,0.2) +
  scale_x_date(breaks=date_breaks("1 month"),labels=date_format("%b")) +
  labs(title="NDVI in 2014",x="Month",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="mean_NDVI_2014.png", width = 8, height = 3, dpi=300)

df_mean_2019 <- subset(df_mean_NDVI,Year==2019,select=c("Mean","Std","Date","Month"))
# Plot the data from 2019
ggplot(data=df_mean_2019,aes(x=as.Date(Date,format="%Y-%m-%d"),
                             y=Mean),na.rm=TRUE) +
  geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4,fill="#35978f") +
  geom_line(size=0.4,color="#35978f") +
  #coord_cartesian(expand=FALSE) +
  ylim(NA,0.2) +
  scale_x_date(breaks=date_breaks("1 month"),labels=date_format("%b")) +
  labs(title="NDVI in 2019",x="Month",y="NDVI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="mean_NDVI_2019.png", width = 8, height = 3, dpi=300)

