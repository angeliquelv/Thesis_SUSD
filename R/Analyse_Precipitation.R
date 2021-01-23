rm(list = ls())  # clear workspace

# Install packages and load libraries -------------------------------------

# install packages if required 
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}
if(!require(lubridate)){install.packages("lubridate")}
if(!require(precintcon)){install.packages("precintcon")}
if(!require(lfstat)){install.packages("lfstat")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(scales)){install.packages("scales")}
if(!require(extrafont)){install.packages("extrafont")}

# load libraries
library(rgdal)
library(raster) # package for raster manipulation
library(lubridate) # package for netcdf manipulation 
library(precintcon) 
library(lfstat)
library(ggplot2)
library(extrafont)


# Set working directory ---------------------------------------------------

workdirdata <- "G:/Thesis_data/Precipitation"
workdir_Fig = "G:/Thesis_data/Figures"

setwd(workdirdata)

# Calculate annual cumulative precipitation -------------------------------

tp <- stack("tp.grd")

# the precipitation data are monthly mean daily precipitation values
# so per month the value should be multiplied by the nr of days 
date_df<-readRDS(file="date_df.Rda")
date_list <- as.Date(date_df[,1], "%Y-%m-%d")
nr_days_month <- lubridate::days_in_month(date_list)
# Calculate the total monthly precipitation per cell in the raster layer 
monthly_precipitation <- stack()
for (i in 1:nlayers(tp)) {
  precipitation <- nr_days_month[i][[1]] * tp[[i]]
  monthly_precipitation <- stack(monthly_precipitation,precipitation)
}

rm(precipitation) # remove unused variable from workspace

# now calculate the annual precipitation values by accumulating monthly total precipitation
annual_precipitation <- stack()
for (i in 1:36) {  
  precipitation <- calc(monthly_precipitation[[((i-1)*12+1):((i-1)*12+12)]],fun=sum)
  annual_precipitation <- stack(annual_precipitation,precipitation)
}

rm(precipitation) # remove unused variable from workspace

names(annual_precipitation) <- 1984:2019

# Now calculate the mean annual precipitation per year in the whole watershed (rasterlayer)
# Results in a numeric vector with 36 values (one for each year)
annual_precipitation_Ounila <- cellStats(annual_precipitation,stat=mean)

# Now calculate the median of these annual precipitation values 
median_annual_precipitation <- median(annual_precipitation_Ounila)

# Now calculate the mean annual precipitation 
mean_annual_precipitation <- mean(annual_precipitation_Ounila)

# Plot the distribution of annual precipitation values in the Ounila watershed
plot(annual_precipitation_Ounila)
# There are 2 extremely dry years (2000 and 2001, <200 mm) and 2 extremely wet years (1989 and 2014 >500 mm)

# Calculate the mean monthly precipitation per year in the whole watershed (rasterlayer)
# Results in a numeric vector with 432 values (one for each month)
monthly_precipitation_Ounila <- cellStats(monthly_precipitation,stat=mean)
monthly_precipitation_Ounila <- cbind(date_df$Date,date_df$Year,date_df$Month,monthly_precipitation_Ounila)
colnames(monthly_precipitation_Ounila) <- c("date","year","month","precipitation")
# now add the hydrological year to the dataframe
monthly_precipitation_Ounila <- as.data.frame(monthly_precipitation_Ounila)
saveRDS(monthly_precipitation_Ounila,"monthly_precipitation_Ounila.Rda")
monthly_precipitation_Ounila$water_year<-water_year(as.Date(substr(rownames(monthly_precipitation_Ounila),2,11),format="%Y.%m.%d"),origin="usgs",ax.POSIX=TRUE)

setwd(workdir_Fig)


# Calculate SPI -----------------------------------------------------------
# Calculate Standardized Precipitation Index SPI per month
monthly_precipitation_Ounila <- as.data.frame(monthly_precipitation_Ounila)
monthly_precipitation_year <- cbind(monthly_precipitation_Ounila$year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation)
monthly_precipitation_year <- as.data.frame(monthly_precipitation_year)
monthly_precipitation_year <- as.precintcon.monthly(monthly_precipitation_year)
SPI_month <- spi(monthly_precipitation_year)
SPI_year <- spi.per.year(monthly_precipitation_year)

monthly_precipitation_Ounila <- as.data.frame(monthly_precipitation_Ounila)
monthly_precipitation_water_year <- cbind(monthly_precipitation_Ounila$water_year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation)
monthly_precipitation_water_year <- as.data.frame(monthly_precipitation_water_year)
monthly_precipitation_water_year <- as.precintcon.monthly(monthly_precipitation_water_year)
SPI_water_year <- spi.per.year(monthly_precipitation_water_year)
SPI_water_year$year <- 1984:2020

monthly_precipitation_year_log <- monthly_precipitation_year
monthly_precipitation_year_log$precipitation <- log(monthly_precipitation_year_log$precipitation)

monthly_precipitation_year_sqrt <- monthly_precipitation_year
monthly_precipitation_year_sqrt$precipitation <- sqrt(monthly_precipitation_year_sqrt$precipitation)

monthly_precipitation_year_cbrt <- monthly_precipitation_year
monthly_precipitation_year_cbrt$precipitation <- (monthly_precipitation_year_cbrt$precipitation)^(1/3)

# Calculate Standardized Precipitation Index SPI per month
monthly_precipitation_Ounila <- as.data.frame(monthly_precipitation_Ounila)
monthly_precipitation_year_cbrt <- cbind(monthly_precipitation_Ounila$year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation)
monthly_precipitation_year_cbrt <- as.data.frame(monthly_precipitation_year_cbrt)
monthly_precipitation_year_cbrt <- as.precintcon.monthly(monthly_precipitation_year_cbrt)
SPI_month_cbrt <- spi(monthly_precipitation_year_cbrt,period=12)
SPI_year_cbrt <- spi.per.year(monthly_precipitation_year_cbrt)

monthly_precipitation_Ounila <- as.data.frame(monthly_precipitation_Ounila)
monthly_precipitation_water_year <- cbind(monthly_precipitation_Ounila$water_year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation)
monthly_precipitation_water_year <- as.data.frame(monthly_precipitation_water_year)
monthly_precipitation_water_year <- as.precintcon.monthly(monthly_precipitation_water_year)

monthly_precipitation_water_year_cbrt <- monthly_precipitation_water_year
monthly_precipitation_water_year_cbrt$precipitation <- (monthly_precipitation_water_year_cbrt$precipitation)^(1/3)

SPI_water_year_cbrt <- spi.per.year(monthly_precipitation_water_year_cbrt,period=12)
SPI_water_year_cbrt$year <- 1985:2020


# Calculate RAI -----------------------------------------------------------
# Calculate Rainfall Anomaly Index (RAI) per year
monthly_precipitation_year <- cbind(as.data.frame(monthly_precipitation_Ounila$year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation))
monthly_precipitation_year <- as.data.frame(monthly_precipitation_year)
monthly_precipitation_year <- as.precintcon.monthly(monthly_precipitation_year)
RAI_year <- rai(monthly_precipitation_year,granularity='a')

RAI_year$sign <- RAI_year$rai>0
setwd(workdir_Fig)
ggplot(data=RAI_year,aes(y=rai,x=year,fill=sign)) +
  geom_col(show.legend=FALSE) +
  labs(title="Rainfall Anomaly Index per year",x="Year",y="RAI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  scale_x_continuous(breaks = round(seq(min(RAI_year$year),max(RAI_year$year),by =2),1)) +
  # scale_fill_manual(values =c())
  ggsave(filename="RAI_year.png", family="Calibri", width = 8, height = 3, dpi=300)

RAI_monthly<- rai(monthly_precipitation_year,granularity='m')
RAI_monthly$sign <- RAI_monthly$rai>0

RAI_monthly$date <- paste0("1-",RAI_monthly$month,"-",RAI_monthly$year)
RAI_monthly$date <- as.Date(RAI_monthly$date,format="%d-%m-%Y")

ggplot(data=RAI_monthly,aes(y=rai,x=date,fill=sign)) +
  geom_col(show.legend=FALSE) +
  coord_cartesian(expand=FALSE) +
  ylim(-3.5,7.5) +
  scale_x_date(breaks=date_breaks("2 years"),labels=date_format("%Y"),limits=c(min(RAI_monthly$date),max(RAI_monthly$date))) +
  labs(title="Rainfall Anomaly Index per month",x="Date",y="RAI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="RAI_monthly.png", family="Calibri", width = 8, height = 3, dpi=300)   

# Now calculate RAI per hydrological year 
monthly_precipitation_water_year <- cbind(monthly_precipitation_Ounila$water_year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation)
monthly_precipitation_water_year <- as.data.frame(monthly_precipitation_water_year)
monthly_precipitation_water_year <- as.precintcon.monthly(monthly_precipitation_water_year)
RAI_water_year <- rai(monthly_precipitation_water_year,granularity='a')

RAI_water_year$sign <- RAI_water_year$rai>0

# remove the year 1984 and 2020 because they are not complete years (Due to conversion to hydrological years)
RAI_water_year <- RAI_water_year[2:36,]
RAI_water_year$year <- 1985:2019
ggplot(data=RAI_water_year,aes(y=rai,x=year,fill=sign)) +
  geom_col(show.legend=FALSE) +
  labs(title="Rainfall Anomaly Index per hydrological year",x="Year",y="RAI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        axis.text.x=element_text(size=9,color="black"),axis.text.y=element_text(size=9,color="black"),
        axis.title.x=element_text(size=11,color="black"),axis.title.y=element_text(size=11,color="black")) +
  scale_x_continuous(breaks = round(seq(min(RAI_water_year$year),max(RAI_water_year$year),by =2),1)) +
  scale_fill_manual(values =c("#BF812D","#35978F"))
  ggsave(filename="RAI_water_year.png", family="Calibri", width = 6.76, height = 2.5, dpi=300)

RAI_monthly_water_year <- rai(monthly_precipitation_water_year,granularity='m')
RAI_monthly_water_year$sign <- RAI_monthly_water_year$rai>0

RAI_monthly_water_year$date <- paste0("1-",RAI_monthly_water_year$month,"-",RAI_monthly_water_year$year)
RAI_monthly_water_year$date <- as.Date(RAI_monthly_water_year$date,format="%d-%m-%Y")

RAI_monthly_water_year <- RAI_monthly_water_year[10:429,]

ggplot(data=RAI_monthly_water_year,aes(y=rai,x=date,fill=sign)) +
  geom_col(show.legend=FALSE) +
  coord_cartesian(expand=FALSE) +
  ylim(-3.5,7.5) +
  scale_x_date(breaks=date_breaks("2 years"),labels=date_format("%Y"),limits=c(min(RAI_monthly_water_year$date),max(RAI_monthly_water_year$date))) +
  labs(title="Rainfall Anomaly Index per month in hydrological years",x="Date",y="RAI") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
  ggsave(filename="RAI_monthly_water_year.png", family="Calibri", width = 8, height = 3, dpi=300)   
  
  ggplot(data=monthly_precipitation_Ounila,aes(y=precipitation,x=date)) +
  geom_col(show.legend=FALSE) +
  labs(title="Precipitation time series ",x="Year",y="precipitation") +
  theme(plot.title=element_text(hjust=0.5,size=12,face="bold"),
        axis.text.x=element_text(size=9),axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),axis.title.y=element_text(size=11),
        legend.title=element_text(size=11), legend.text=element_text(size=9)) +
    scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
    # scale_fill_manual(values =c())
  ggsave(filename="Precipitation_timeseries.png", width = 8, height = 3, dpi=300)