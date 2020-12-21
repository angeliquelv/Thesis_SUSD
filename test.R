##Check examples ResInd

library(resInd)
library(bfastSpatial)

workdirNDVI <- "G:/Thesis_data/NDVI"

setwd(workdirNDVI)

NDVI<-brick("NDVI_final_OLS.grd")

# When excluding negative values the warning In sqrt(fr) : NaNs produced is avoided in SOME cases 
### Turns out it actually doesnt matter that much, I just get this sometimes

dates<-as.Date(substr(names(NDVI),5,12),format="%Y%m%d")

# pixels on which I will test
plot(NDVI[[1]],addfun=points(x=c(-7.144,-7.303,-7.148,-7.279,-7.103),y=c(31.092,31.228,31.274,31.293,31.337),col="red",pch=22,cex=1,bg="red"))
text(-7.144,31.092,1,pos=4,offset=0.5,font=2,cex=0.8)
text(-7.303,31.228,2,pos=4,offset=0.5,font=2,cex=0.8)
text(-7.148,31.274,3,pos=4,offset=0.5,font=2,cex=0.8)
text(-7.279,31.293,4,pos=4,offset=0.5,font=2,cex=0.8)
text(-7.103,31.337,5,pos=4,offset=0.5,font=2,cex=0.8)

raster::cellFromXY(NDVI,cbind(c(-7.144,-7.303,-7.148,-7.279,-7.103),c(31.092,31.228,31.274,31.293,31.337)))##Rund function resIndSpatial on test raster set
# Output, these are the targetCells
# [1] 1438029  754302  523424  427738  207815  986018  157074


##Select target pixel from raster stack
targcell1 <- 1438029
targcell2 <- 754302
targcell3 <- 523424
targcell4 <- 427738
targcell5 <- 207815

x1 <- as.vector(NDVI[targcell1])
x2 <- as.vector(NDVI[targcell2])
x3 <- as.vector(NDVI[targcell3])
x4 <- as.vector(NDVI[targcell4])
x5 <- as.vector(NDVI[targcell5])


##Extract vector of NDVI values for targcell and plot time series.
plot(dates,x, ylab='NDVI', xlab='year', main='targcell', ylim=c(0,1))

decimal_date(as.Date("1-10-1998",format="%d-%m-%Y"))
decimal_date(as.Date("30-09-2002",format="%d-%m-%Y"))

windows()

 y1 <- resInd(x1, dates, type='irregular', sc=1, order=2, 
              formula = response ~ (trend + harmon),
              h=0.05, plevel=0.05,
              dr=c(1998.748,2002.745), drd=1998.748, 
              s=3,plot=TRUE)
 y2 <- resInd(x2, dates, type='irregular', sc=1, order=2, 
              formula = response ~ (trend + harmon),
              h=0.05, plevel=0.05,
              dr=c(1998.748,2002.745), drd=1998.748, 
              s=3,plot=TRUE)
 y3 <- resInd(x3, dates, type='irregular', sc=1, order=2, 
              formula = response ~ (trend + harmon),
              h=0.05, plevel=0.05,
              dr=c(1998.748,2002.745), drd=1998.748, 
              s=3,plot=TRUE)
 y4 <- resInd(x4, dates, type='irregular', sc=1, order=2, 
              formula = response ~ (trend + harmon),
              h=0.05, plevel=0.05,
              dr=c(1998.748,2002.745), drd=1998.748, 
              s=3,plot=TRUE)
 y5 <- resInd(x5, dates, type='irregular', sc=1, order=2, 
              formula = response ~ (trend + harmon),
              h=0.05, plevel=0.05,
              dr=c(1998.748,2002.745), drd=1998.748, 
              s=3,plot=TRUE)

 
 