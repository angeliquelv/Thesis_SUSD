# Input: output and NDVI rasterbrick

rm(list = ls())  # clear workspace

# Install required packages and load libraries  ---------------------------

# Install requried packages
if(!require(raster)){install.packages("raster")}
if(!require(parallel)){install.packages("parallel")}
if(!require(strucchange)){install.packages("strucchange")}
if(!require(MASS)){install.packages("MASS")}
if(!require(zoo)){install.packages("zoo")}
if(!require(rgdal)){install.packages("rgdal")}
if(!require(rgeos)){install.packages("rgeos")}
if(!require(gdalUtils)){install.packages("gdalUtils")}
if(!require(bfast)){install.packages("bfast")}
if(!require(bfastSpatial)){install.packages("bfastSpatial")}
if(!require(extrafont)){install.packages("extrafont")}
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(lubridate)){install.packages("lubridate")}
if(!require(scales)){install.packages("scales")}
if(!require(viridis)){install.packages("viridis")}
if(!require(colortools)){install.packages("colortools")}
if(!require(gridExtra)){install.packages("gridExtra")}
if(!require(colorspace)){install.packages("colorspace")}
if(!require(precintcon)){install.packages("precintcon")}
if(!require(SpatialPack)){install.packages("spatialEco")}
if(!require(lfstat)){install.packages("lfstat")}


# Load libraries 
library(raster) # package for raster manipulation
library(parallel)
library(strucchange)
library(MASS)
library(zoo)
library(rgdal)
library(rgeos)
library(gdalUtils)
library(bfast)
library(bfastSpatial)
library(extrafont)
library(RColorBrewer)
library(ggplot2)
library(lubridate)
library(scales)
library(viridis)
library(colortools)
library(gridExtra)
library(colorspace)
library(precintcon)
library(spatialEco)
library(lfstat)

# Set working directories and load data  ----------------------------------

workdir_Fig = "G:/Thesis_data/Figures"
workdir_output = "G:/Thesis_data/output_resIndSpatial"
workdir_NDVI = "G:/Thesis_data/NDVI"
workdir_Prec <- "G:/Thesis_data/Precipitation"
workdir_LandUse = "G:/Thesis_data/Land_Use"

setwd(workdir_output)

output <- brick("output.grd")

setwd(workdir_NDVI)

NDVI <- brick("NDVI_final_OLS.grd")

df_mean_NDVI <- readRDS("df_mean_NDVI_OLS.Rda")

setwd(workdir_LandUse)

tif_file <- list.files(path = workdir_LandUse, pattern = '.tif$')

Land_Use <- raster(tif_file[1])
Land_Use<-projectRaster(Land_Use,output,method="ngb")
Land_Use_Ounila<-mask(Land_Use,output$BPNumb)

setwd(workdir_Fig)


# Colours -----------------------------------------------------------------

pos_pal<-c("grey85",rev(paste0(viridis(n=6,option="D"))))
neg_pal <-c("grey85",rev(paste0(viridis(n=6,option="plasma"))))

bp_pal_vis <- c("grey85","#575756",
                rev(paste0(viridis(n=5,option="D"))[1:4]),
                rev(paste0(viridis(n=6,option="plasma"))[3:6]))

bp_pal <- c("#575756",
            rev(paste0(viridis(n=5,option="D"))[3:4]),viridis(n=6,option="plasma")[4],
            rev(paste0(viridis(n=6,option="plasma"))[5:6]),viridis(n=5,option="D")[2],
            viridis(n=5,option="D")[1],viridis(n=6,option="plasma")[3],
            "grey85")

pizza(c(paste0(viridis(n=6,option="D")[2:5]),"#003C30",
         "#F0CE0A","#F05B0A","#E8190C","#FF1398","#543005","#575756","#9B9B9B","grey85"))

# Example plots breakpoint typologies  ------------------------------------

df_output <- as.data.frame(output, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no breakpoints at all 
df_output_no_bp <- subset(df_output,BPNumb==0) # 11241 pixels
nrow(df_output_no_bp)

# observations with no drought breakpoint
df_output_no_dbp <- subset(df_output,DBP==0) # 178045 pixels
nrow(df_output_no_dbp)

# observations with drought breakpoint
df_output_dbp <- subset(df_output,DBP==1) # 629846 pixels
nrow(df_output_dbp)

# Example plots without breakpoints:
# Run the following function for each categorie on one of the pixels without breakpoitns
# Within the function change the name of the png file and the location of the significance asterisk (*) as necessary, 
# to create a pretty picture
resInd <- function(x, dates, type='irregular', sc=1, order=3,
                   formula = response ~ (trend + harmon), h=0.15,
                   plevel=0.05, dr, drd, s=3, NV=NA, plot=FALSE) {
  #Set output vector to NA
  resind <- NA
  #Set own NoData value to differentiate between "no data pixels" and "no breakpoint pixels"
  NV=NV
  #Apply bfastts() and multiply data by scaling factor if required
  bfts <- bfastts(x*sc,dates,type=type)
  ##If not all oberservation NA (e.g. pixels outside area of interest) run procedure
  if(!all(is.na(bfts))){
    #Create data frame from bfastts
    bpp <- bfastpp(bfts, order=order, stl=("none"), na.action=na.omit)
    #Calculate initial NDVI: First s years after start of observations
    Ini <- subset(bpp, time >= bpp$time[1] & time <= bpp$time[1]+s)
    MIni <- mean(Ini$response)
    ##MOSUM test for structural stability
    Mos <- efp(formula, data=bpp, type="OLS-MOSUM", h=h, rescale=FALSE)
    #Calculate test statistics
    sct <- sctest(Mos)
    ##Fit breakpoints and segmented model if empirical fluctuation test was significant
    if (sct$p.value <= plevel){
      #Fit the breakpoints
      bpoints <- breakpoints(formula=formula, data=bpp, h=h)
      p <- bpoints$breakpoints
      if((length(p)==1 && is.na(p))) {
        #If no breakpoint p = NA. Set breakpoint number to zero. Give warning.
        warning('No breakpoint found')
        bpnumb <- 0
        #Set breakpoint parameters to NA (needed for further calculations)
        bd <- NA
        cid <- NA
        #Fit unsegmented model
        m <- rlm (formula=formula, data=bpp, maxit=100)
        bpp$segment <- 1 #even if you have only one segment you need to specify
        #this, for trend calculation later on
      } else {
        #If there is a breakpoint extract breakpoint parameters
        #breakpoint number
        bpnumb <- length(p)
        #breakdates
        bd <- bpp$time[p]
        #confidence invervals
        ci <- confint(bpoints, breakpoints=TRUE)
        bda <- ci[[1]]
        cid <- bpp$time[c(bda[,1],bda[,3])]
        ##Fit segmented model
        #Create column "segment" that indices segment in the dataframe "bpp"
        bpp$segment <- breakfactor(bpoints)
        #RLM fit; I increased the default maxiterations (20 -> 100)
        m <- rlm (response ~ segment/(trend+harmon), data=bpp, maxit=100)
      }
    } else {
      ##If MOSUM not significant set breakpoint number and DBP to zero and other
      #output variables to NV
      bpnumb <- 0
      #Fit unsegmented model
      m <- rlm (formula=formula, data=bpp, maxit=100)
      bpp$segment <- 1 ##even if you have only one segment you need to specify this,
      #for trend calculation later on
      warning('No significant deviation from structural stability (MOSUM test).No breakpoints fitted.')
    }
    #Predict values and add column "prediction" in in the dataframe "bpp"
    bpp$prediction <- predict(m,newdata=bpp)
    #Add trend prediction for each segment. Corrected hight of trend line for
    #irregular data: sets harmonic term based on mean DOY of observations/segment
    # (instead of trend$harmon[] <- 0, which assumes regular data);
    # idea based on email exchange with Achim Zeileis
    for(i in 1:length(unique(bpp$segment))) {
      seg <- unique(bpp$segment)[i]
      #trend <- subset(bpp, segment == seg)
      trend <- bpp[bpp$segment == seg,]
      trend$days <- as.numeric(substr(formatC(trend$time,format='f',digits=3), 6, 8))
      dmean <- mean(trend$days)
      har <- dmean/365
      trend$harmon[] <- rep(
        c(cos(2 * pi * har * 1:order), sin(2 * pi * har * 1:order)),
        each = nrow(trend))
      #Add trendprediction column in bpp matrix
      bpp$trendprediction[bpp$segment == seg] <- predict(m, newdata = trend)
    }
    ##Add column for residuals to dataframe bpp
    # bpp$residuals <- as.numeric(bpp$response)-as.numeric(bpp$prediction)
    ##Extract coefficients
    coef <- coefficients(m)
    ##Save intercept
    Int <- coef[1]
    ##Extract trends
    trends <- coef[which(grepl('trend', names(coef)))]
    ##Calculate confidence intervals around trends
    ci_t <- confint.default(m)
    ci_trend <- ci_t[which(grepl('trend', rownames(ci_t))),]
    ##Extract Amplitudes for each segment. Depending on parameter "order" one gets
    #multiple sine and cosine terms
    cosines <- coef[which(grepl('harmoncos', names(coef)))]
    sines <- coef[which(grepl('harmonsin', names(coef)))]
    amps <- sqrt(cosines^2 + sines^2)
    names(amps) <- rep(1:length(unique(bpp$segment)), order)
    # Calculate (drought) resilience indicators -------------------------------
    ##If no (significant) breakpoint was found set breakpoint position,
    # breakpoint timing, and recovery trend to "NV"; set DBP to 0
    if(bpnumb==0){
      DBP <- 0
      DBP_Type <- NV
      bpt <- NV
      tlag <- NV
      trendrecov <- NV
      pretrend <- NV
      preNDVI <- NV
      MagA <- NV
      MagR <- NV
      MagTA <- NV
      MagTR <- NV
      AmpDiff <- NV
      #Set breakpoint position to NA (needed for further calculations)
      bpd <- NA
    } else {
      ##If breakpoint, check if one breakpoint ocurred around drought period
      #("drought breakpoint")
      bpd <- match(bd[bd >= dr[1] & bd <= dr[2]], bd) #bpd gives you the position of bp.
      ##If drought breakpoint occured set DBP to 1 and calculate breakpoint timing,
      #trend of recovery, trend in segment before BP, preNDVI,
      #  Magnitude of change (MagA & MagR) based on mean NDVI of s years before and after BP.
      #If there is more than one BP in the specified time interval, the first one is chosen
      if (length(bpd) > 1) {
        bpd <- bpd[1]
      }
      if (length(bpd) > 0) {
        DBP <- 1
        #Calculate BP timing
        bpt <- bd[bpd]
        #Calculate tlag (days between drought reference year and bpt in days,
        # assuming 365 days/year)
        tlag <- (bpt-drd)*365
        #Calculate trend slopes
        seg2 <- (bpd+1)
        trendrecov <- trends[seg2]
        seg1 <- (bpd)
        pretrend <- trends[seg1]
        #Calculate Confidence Intervals trends slopes
        ci_trendrecov <- ci_trend[seg2,]
        ci_pretrend <- ci_trend[seg1,]
        #Calculate preNDVI, Magnitude of change (MagA & MagR) based on mean NDVI
        #of s years before and after BP. Use subset of bpp$response that falls
        #within date vector
        ti <- c(bpt-s, bpt+s) # +/- s years around breakpoint
        S1 <- subset(bpp, time >= ti[1] & time <= bpt)
        S2 <- subset(bpp, time > bpt & time <= ti[2])
        preNDVI <- mean(S1$response)
        M2 <- mean(S2$response)
        MagA <- (M2-preNDVI) #absolute change magnitude
        MagR <- MagA/preNDVI #relative change magnitude
        #Calculate Magnitude of change based on corrected trend prediction (MagTA & MagTR)
        w <- which(bpp$time==bpt) #gives you position of breakpoint in trend dataframe
        y1 <- bpp$trendprediction[w] #trend prediction value at time of breakpoint
        y2 <- bpp$trendprediction[w+1] #trend prediction value of observation after breakpoint
        MagTA <- y2-y1 #absolute breakpoint magnitude
        MagTR <- (y2-y1)/y1 #relative breakpoint magnitude
        #Calculate Difference in mean amplitudes between segments
        mean_ampsseg1 <- mean(amps[which(grepl(seg1,names(amps)))]) #mean amplitude in segment before bp
        mean_ampsseg2 <- mean(amps[which(grepl(seg2,names(amps)))]) #mean amplitude in segment after bp
        AmpDiff <- (mean_ampsseg2-mean_ampsseg1)/mean_ampsseg1 #relative difference between mean amplitudes
        # Determine the breakpoint typology
        # no significant trend before and after breakpoint
        if (prod(ci_pretrend)<0 && prod(ci_trendrecov)<0) {
          DBP_Type <- 0
        }
        # interrupted increase, significant increase before and after breakpoint
        if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend>0 && trendrecov>0) {
          DBP_Type <- 1
        }
        # increase after period of no significant trend
        if (prod(ci_pretrend)<0 && prod(ci_trendrecov)>0 && trendrecov>0) {
          DBP_Type <- 2
        }
        # period of no significant trend after increase
        if (prod(ci_pretrend)>0 && prod(ci_trendrecov)<0  && pretrend>0) {
          DBP_Type <- 3
        }
        # interrupted decrease, significant decrease before and after breakpoint
        if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend<0 && trendrecov<0) {
          DBP_Type <- 4
        }
        # decrease after period of no significant trend
        if (prod(ci_pretrend)<0 && prod(ci_trendrecov)>0 && trendrecov<0) {
          DBP_Type <- 5
        }
        # period of no significant trend after decrease
        if (prod(ci_pretrend)>0 && prod(ci_trendrecov)<0 && pretrend<0) {
          DBP_Type <- 6
        }
        # positive reversal
        if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend<0 && trendrecov>0) {
          DBP_Type <- 7
        }
        # negative reversal
        if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend>0 && trendrecov<0) {
          DBP_Type <- 8
        }
      } else {
        ##If no breakpoint around drought occured set patrameters to NA
        warning('No drought breakpoint found')
        DBP <- 0
        DBP_Type <- NV
        bpt <- NV
        tlag <- NV
        trendrecov <- NV
        pretrend <- NV
        preNDVI <- NV
        MagA <- NV
        MagR <- NV
        MagTA <- NV
        MagTR <- NV
        AmpDiff <- NV
      }
    }
  }else {
    ##If all values NA (pixels without NDVI values!) set output variables to NA and give warning
    warning('No observations in time series')
    bpnumb <- NA
    MIni <- NA
    Int <- NA
    DBP <- NA
    DBP_Type <- NA
    bpt <- NA
    tlag <- NA
    trendrecov <- NA
    pretrend <- NA
    preNDVI <- NA
    MagA <- NA
    MagR <- NA
    MagTA <- NA
    MagTR <- NA
    AmpDiff <- NA
  }
  # Plotting ----------------------------------------------------------------
  if(plot) {
    plot.new()
    png(filename="name.png", family="Calibri",width=1350, height=900, res=300)
    par(mfrow = c(1,1), mar=c(0,3,0,0)+0.1)
    plot_bpp<-subset(bpp,time >= 1997.000 & time <= 2007.000)
    plot(plot_bpp$time,plot_bpp$response,xaxt='n',ann=FALSE,cex.axis=1,type='l',col="grey",lwd=3);
    mtext("NDVI",side=2,line=2,font=1,cex=1.2)
    points(plot_bpp$time, plot_bpp$prediction,col='blue',type='l',lwd=3);
    #Add vertical line at breakdates
    text(x=2001.000,y=0.25,"*",cex=2,font=2)
    # plot trend
    lines(zoo::zoo(plot_bpp$trendprediction[plot_bpp$segment == unique(bpp$segment)],
                   plot_bpp$time[plot_bpp$segment == unique(bpp$segment)]), col = 'red',lwd=3)

    dev.off()
  }
  # Save and return output --------------------------------------------------
  #Save all indicators:
  resind <- cbind(bpnumb, MIni, as.numeric(Int), DBP, DBP_Type, bpt, tlag,
                  as.numeric(trendrecov), as.numeric(pretrend), preNDVI, MagA,
                  MagR, MagTA, MagTR, AmpDiff)
  colnames(resind) <- c('BPNumb', 'Initial NDVI', 'Intercept', 'DBP','DBP_Type','BpTime',
                        'Timelag', 'RecTrend', 'PreTrend', 'PreNDVI', 'MagObsA',
                        'MagObsR', 'MagTrendA', 'MagTrendR', 'AmpDiffR')
  return(resind)
}

# change the number between brackets to the cell nr you want to make the bfast plot for 
x <- as.vector(NDVI[627])

resInd(x,dates,type='irregular', sc=1, order=2,
       formula = response ~ (trend + harmon),
       h=0.10, plevel=0.05,
       dr=c(1998.748,2002.745), drd=1998.748,
       s=3,plot=TRUE)

# Example plot no trend 
df_output_no_bp_no_trend <- subset(df_output_no_bp,Trend_noBP_lowerCI*Trend_noBP_upperCI<0) # 2709 pixels 
nrow(df_output_no_bp_no_trend)
rownames(df_output_nobp_notrend)[1:10] # show first 10 cellnumbers that have no trend and no breakpoint

# Example plot stable increasing trend 
df_output_no_bp_pos_trend <- subset(df_output_no_bp,Trend_noBP_lowerCI*Trend_noBP_upperCI>0 & Trend_noBP>0) # 3651 pixels 
nrow(df_output_no_bp_pos_trend)
rownames(df_output_no_bp_pos_trend)[1:10] # show first 10 cellnumbers that have a stable increasing trend and no breakpoint

# Example plot stable decreasing trend 
df_output_no_bp_neg_trend <- subset(df_output_no_bp,Trend_noBP_lowerCI*Trend_noBP_upperCI>0 & Trend_noBP<0) # 4881 pixels 
nrow(df_output_no_bp_neg_trend)
rownames(df_output_no_bp_neg_trend)[1:10] # show first 10 cellnumbers that a stable decreasing trend and no breakpoint

# Example plots with breakpoints:
# Run the following function for each categorie on one of the pixels that has a drought breakpoint following that category,
# Within the function change the name of the png file and the location of the significane asterisk (*) as necessary, 
# to create a pretty picture
resInd <- function(x, dates, type='irregular', sc=1, order=3,
                   formula = response ~ (trend + harmon), h=0.15,
                   plevel=0.05, dr, drd, s=3, NV=NA, plot=FALSE) {
  #Set output vector to NA
  resind <- NA
  #Set own NoData value to differentiate between "no data pixels" and "no breakpoint pixels"
  NV=NV
  #Apply bfastts() and multiply data by scaling factor if required
  bfts <- bfastts(x*sc,dates,type=type)
  ##If not all oberservation NA (e.g. pixels outside area of interest) run procedure
  if(!all(is.na(bfts))){
    #Create data frame from bfastts
    bpp <- bfastpp(bfts, order=order, stl=("none"), na.action=na.omit)
    #Calculate initial NDVI: First s years after start of observations
    Ini <- subset(bpp, time >= bpp$time[1] & time <= bpp$time[1]+s)
    MIni <- mean(Ini$response)
    ##MOSUM test for structural stability
    Mos <- efp(formula, data=bpp, type="OLS-MOSUM", h=h, rescale=FALSE)
    #Calculate test statistics
    sct <- sctest(Mos)
    ##Fit breakpoints and segmented model if empirical fluctuation test was significant
    if (sct$p.value <= plevel){
      #Fit the breakpoints
      bpoints <- breakpoints(formula=formula, data=bpp, h=h)
      p <- bpoints$breakpoints
      if((length(p)==1 && is.na(p))) {
        #If no breakpoint p = NA. Set breakpoint number to zero. Give warning.
        warning('No breakpoint found')
        bpnumb <- 0
        #Set breakpoint parameters to NA (needed for further calculations)
        bd <- NA
        cid <- NA
        #Fit unsegmented model
        m <- rlm (formula=formula, data=bpp, maxit=100)
        bpp$segment <- 1 #even if you have only one segment you need to specify
        #this, for trend calculation later on
        } else {
          #If there is a breakpoint extract breakpoint parameters
          #breakpoint number
          bpnumb <- length(p)
          #breakdates
          bd <- bpp$time[p]
          #confidence invervals
          ci <- confint(bpoints, breakpoints=TRUE)
          bda <- ci[[1]]
          cid <- bpp$time[c(bda[,1],bda[,3])]
          ##Fit segmented model
          #Create column "segment" that indices segment in the dataframe "bpp"
          bpp$segment <- breakfactor(bpoints)
          #RLM fit; I increased the default maxiterations (20 -> 100)
          m <- rlm (response ~ segment/(trend+harmon), data=bpp, maxit=100)
          }
      } else {
        ##If MOSUM not significant set breakpoint number and DBP to zero and other
        #output variables to NV
        bpnumb <- 0
        #Fit unsegmented model
        m <- rlm (formula=formula, data=bpp, maxit=100)
        bpp$segment <- 1 ##even if you have only one segment you need to specify this,
        #for trend calculation later on
        warning('No significant deviation from structural stability (MOSUM test).No breakpoints fitted.')
        }
    #Predict values and add column "prediction" in in the dataframe "bpp"
    bpp$prediction <- predict(m,newdata=bpp)
    #Add trend prediction for each segment. Corrected hight of trend line for
    #irregular data: sets harmonic term based on mean DOY of observations/segment
    # (instead of trend$harmon[] <- 0, which assumes regular data);
    # idea based on email exchange with Achim Zeileis
    for(i in 1:length(unique(bpp$segment))) {
      seg <- unique(bpp$segment)[i]
      #trend <- subset(bpp, segment == seg)
      trend <- bpp[bpp$segment == seg,]
      trend$days <- as.numeric(substr(formatC(trend$time,format='f',digits=3), 6, 8))
      dmean <- mean(trend$days)
      har <- dmean/365
      trend$harmon[] <- rep(
        c(cos(2 * pi * har * 1:order), sin(2 * pi * har * 1:order)),
        each = nrow(trend))
      #Add trendprediction column in bpp matrix
      bpp$trendprediction[bpp$segment == seg] <- predict(m, newdata = trend)
      }
    ##Add column for residuals to dataframe bpp
    # bpp$residuals <- as.numeric(bpp$response)-as.numeric(bpp$prediction)
    ##Extract coefficients
    coef <- coefficients(m)
    ##Save intercept
    Int <- coef[1]
    ##Extract trends
    trends <- coef[which(grepl('trend', names(coef)))]
    ##Calculate confidence intervals around trends
    ci_t <- confint.default(m)
    ci_trend <- ci_t[which(grepl('trend', rownames(ci_t))),]
    ##Extract Amplitudes for each segment. Depending on parameter "order" one gets
    #multiple sine and cosine terms
    cosines <- coef[which(grepl('harmoncos', names(coef)))]
    sines <- coef[which(grepl('harmonsin', names(coef)))]
    amps <- sqrt(cosines^2 + sines^2)
    names(amps) <- rep(1:length(unique(bpp$segment)), order)
    # Calculate (drought) resilience indicators -------------------------------
    ##If no (significant) breakpoint was found set breakpoint position,
    # breakpoint timing, and recovery trend to "NV"; set DBP to 0
    if(bpnumb==0){
      DBP <- 0
      DBP_Type <- NV
      bpt <- NV
      tlag <- NV
      trendrecov <- NV
      pretrend <- NV
      preNDVI <- NV
      MagA <- NV
      MagR <- NV
      MagTA <- NV
      MagTR <- NV
      AmpDiff <- NV
      #Set breakpoint position to NA (needed for further calculations)
      bpd <- NA
      } else {
        ##If breakpoint, check if one breakpoint ocurred around drought period
        #("drought breakpoint")
        bpd <- match(bd[bd >= dr[1] & bd <= dr[2]], bd) #bpd gives you the position of bp.
        ##If drought breakpoint occured set DBP to 1 and calculate breakpoint timing,
        #trend of recovery, trend in segment before BP, preNDVI,
        #  Magnitude of change (MagA & MagR) based on mean NDVI of s years before and after BP.
        #If there is more than one BP in the specified time interval, the first one is chosen
        if (length(bpd) > 1) {
          bpd <- bpd[1]
          }
        if (length(bpd) > 0) {
          DBP <- 1
          #Calculate BP timing
          bpt <- bd[bpd]
          #Calculate tlag (days between drought reference year and bpt in days,
          # assuming 365 days/year)
          tlag <- (bpt-drd)*365
          #Calculate trend slopes
          seg2 <- (bpd+1)
          trendrecov <- trends[seg2]
          seg1 <- (bpd)
          pretrend <- trends[seg1]
          #Calculate Confidence Intervals trends slopes
          ci_trendrecov <- ci_trend[seg2,]
          ci_pretrend <- ci_trend[seg1,]
          #Calculate preNDVI, Magnitude of change (MagA & MagR) based on mean NDVI
          #of s years before and after BP. Use subset of bpp$response that falls
          #within date vector
          ti <- c(bpt-s, bpt+s) # +/- s years around breakpoint
          S1 <- subset(bpp, time >= ti[1] & time <= bpt)
          S2 <- subset(bpp, time > bpt & time <= ti[2])
          preNDVI <- mean(S1$response)
          M2 <- mean(S2$response)
          MagA <- (M2-preNDVI) #absolute change magnitude
          MagR <- MagA/preNDVI #relative change magnitude
          #Calculate Magnitude of change based on corrected trend prediction (MagTA & MagTR)
          w <- which(bpp$time==bpt) #gives you position of breakpoint in trend dataframe
          y1 <- bpp$trendprediction[w] #trend prediction value at time of breakpoint
          y2 <- bpp$trendprediction[w+1] #trend prediction value of observation after breakpoint
          MagTA <- y2-y1 #absolute breakpoint magnitude
          MagTR <- (y2-y1)/y1 #relative breakpoint magnitude
          #Calculate Difference in mean amplitudes between segments
          mean_ampsseg1 <- mean(amps[which(grepl(seg1,names(amps)))]) #mean amplitude in segment before bp
          mean_ampsseg2 <- mean(amps[which(grepl(seg2,names(amps)))]) #mean amplitude in segment after bp
          AmpDiff <- (mean_ampsseg2-mean_ampsseg1)/mean_ampsseg1 #relative difference between mean amplitudes
          # Determine the breakpoint typology
          # no significant trend before and after breakpoint
          if (prod(ci_pretrend)<0 && prod(ci_trendrecov)<0) {
            DBP_Type <- 0
            }
          # interrupted increase, significant increase before and after breakpoint
          if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend>0 && trendrecov>0) {
            DBP_Type <- 1
            }
          # increase after period of no significant trend
          if (prod(ci_pretrend)<0 && prod(ci_trendrecov)>0 && trendrecov>0) {
            DBP_Type <- 2
            }
          # period of no significant trend after increase
          if (prod(ci_pretrend)>0 && prod(ci_trendrecov)<0  && pretrend>0) {
            DBP_Type <- 3
            }
          # interrupted decrease, significant decrease before and after breakpoint
          if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend<0 && trendrecov<0) {
            DBP_Type <- 4
            }
          # decrease after period of no significant trend
          if (prod(ci_pretrend)<0 && prod(ci_trendrecov)>0 && trendrecov<0) {
            DBP_Type <- 5
            }
          # period of no significant trend after decrease
          if (prod(ci_pretrend)>0 && prod(ci_trendrecov)<0 && pretrend<0) {
            DBP_Type <- 6
            }
          # positive reversal
          if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend<0 && trendrecov>0) {
            DBP_Type <- 7
            }
          # negative reversal
          if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend>0 && trendrecov<0) {
            DBP_Type <- 8
            }
          } else {
            ##If no breakpoint around drought occured set patrameters to NA
              warning('No drought breakpoint found')
            DBP <- 0
            DBP_Type <- NV
            bpt <- NV
            tlag <- NV
            trendrecov <- NV
            pretrend <- NV
            preNDVI <- NV
            MagA <- NV
            MagR <- NV
            MagTA <- NV
            MagTR <- NV
            AmpDiff <- NV
          }
        }
    }else {
      ##If all values NA (pixels without NDVI values!) set output variables to NA and give warning
      warning('No observations in time series')
      bpnumb <- NA
      MIni <- NA
      Int <- NA
      DBP <- NA
      DBP_Type <- NA
      bpt <- NA
      tlag <- NA
      trendrecov <- NA
      pretrend <- NA
      preNDVI <- NA
      MagA <- NA
      MagR <- NA
      MagTA <- NA
      MagTR <- NA
      AmpDiff <- NA
      }
  # Plotting ----------------------------------------------------------------
  if(plot) {
    plot.new()
    png(filename="name.png", family="Calibri",width=1350, height=900, res=300)
    par(mfrow = c(1,1), mar=c(0,3,0,0)+0.1)
    plot_bpp<-subset(bpp,time >= (bpt-5) & time <= (bpt+5))
    plot(plot_bpp$time,plot_bpp$response,xaxt='n',ann=FALSE,cex.axis=1,type='l',col="grey",lwd=3);
    mtext("NDVI",side=2,line=2,font=1,cex=1.2)
    points(plot_bpp$time, plot_bpp$prediction,col='blue',type='l',lwd=3);
    #Add vertical line at breakdates
    abline(v = bd, lty = 3, col="black", lwd=3);
    text(x=1997.000,y=0.25,"*",cex=2,font=2)
    text(x=1997.000,y=0.22,"*",cex=2,font=2)
    # plot pretrend
    lines(zoo::zoo(plot_bpp$trendprediction[plot_bpp$segment == unique(bpp$segment)[seg1]],
                   plot_bpp$time[plot_bpp$segment == unique(bpp$segment)[seg1]]), col = 'red',lwd=3)
    # plot recovtrend
    lines(zoo::zoo(plot_bpp$trendprediction[plot_bpp$segment == unique(bpp$segment)[seg2]],
                   plot_bpp$time[plot_bpp$segment == unique(bpp$segment)[seg2]]), col = 'red',lwd=3)
    dev.off()
    }
  # Save and return output --------------------------------------------------
  #Save all indicators:
  resind <- cbind(bpnumb, MIni, as.numeric(Int), DBP, DBP_Type, bpt, tlag,
                  as.numeric(trendrecov), as.numeric(pretrend), preNDVI, MagA,
                  MagR, MagTA, MagTR, AmpDiff)
  colnames(resind) <- c('BPNumb', 'Initial NDVI', 'Intercept', 'DBP','DBP_Type','BpTime',
                        'Timelag', 'RecTrend', 'PreTrend', 'PreNDVI', 'MagObsA',
                        'MagObsR', 'MagTrendA', 'MagTrendR', 'AmpDiffR')
  return(resind)
}

# change the number between brackets to the cell nr you want to make the bfast plot for 
x <- as.vector(NDVI[11842])

resInd(x,dates,type='irregular', sc=1, order=2,
       formula = response ~ (trend + harmon),
       h=0.10, plevel=0.05,
       dr=c(1998.748,2002.745), drd=1998.748,
       s=3,plot=TRUE)

# Example plot Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_type0 <- subset(df_output,DBPType==0) # 11484 pixels
nrow(df_output_type0)
df_output_type0 <- subset(df_output_type0,BPNumb==1)
rownames(df_output_type0)[1:10] # show first 10 cellnumbers that have a type 0 drought breakpoint

# Example plot Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_type1 <- subset(df_output,DBPType==1) # 38705 pixels
nrow(df_output_type1)
df_output_type1 <- subset(df_output_type1,BPNumb==1)
rownames(df_output_type1)[1:10] # show first 10 cellnumbers that have a type 1 drought breakpoint

# Example plot Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_type2 <- subset(df_output,DBPType==2) # 54137 pixels
nrow(df_output_type2)
df_output_type2 <- subset(df_output_type2,BPNumb==1)
rownames(df_output_type2)[1:10] # show first 10 cellnumbers that have a type 2 drought breakpoint

# Example plot Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_type3 <- subset(df_output,DBPType==3) # 12382 pixels
nrow(df_output_type3)
df_output_type3 <- subset(df_output_type3,BPNumb==1)
rownames(df_output_type3)[1:10] # show first 10 cellnumbers that have a type 3 drought breakpoint

# Example plot Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_type4 <- subset(df_output,DBPType==4) # 136780 pixels
nrow(df_output_type4)
df_output_type4 <- subset(df_output_type4,BPNumb==1)
rownames(df_output_type4)[1:10] # show first 10 cellnumbers that have a type 4 drought breakpoint

# Example plot Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_type5 <- subset(df_output,DBPType==5) # 27226 pixels
nrow(df_output_type5)
df_output_type5 <- subset(df_output_type5,BPNumb==1)
rownames(df_output_type5)[1:10] # show first 10 cellnumbers that have a type 5 drought breakpoint

# Example plot Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_type6 <- subset(df_output,DBPType==6) # 42435 pixels
nrow(df_output_type6)
df_output_type6 <- subset(df_output_type6,BPNumb==1)
rownames(df_output_type6)[1:10] # show first 10 cellnumbers that have a type 6 drought breakpoint

# Example plot Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_type7 <- subset(df_output,DBPType==7) # 287668 pixels
nrow(df_output_type7)
df_output_type7 <- subset(df_output_type7,BPNumb==1)
rownames(df_output_type7)[1:10] # show first 10 cellnumbers that have a type 7 drought breakpoint

# Example plot Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_type8 <- subset(df_output,DBPType==8) # 19029 pixels
nrow(df_output_type8)
df_output_type8 <- subset(df_output_type8,BPNumb==1)
rownames(df_output_type8)[1:10] # show first 10 cellnumbers that have a type 8 drought breakpoint



# Breakpoints and precipitation over time  --------------------------------

valid <- subset(output, 1) # use the first layer (BPNumb) to check how many valid observations 
# BPNumb is only NA when there is no valid data for that pixel

#Extract number of valid Pixels - needed later
vP <- sum(!is.na(values(valid))) # 807891 pixels

bdates_output <- subset(output, 17:25) # subset of breakdates
lowerci_output <- subset(output, 26:34) # subset of lower ci limit around breakdate
upperci_output <- subset(output, 35:43) # subset of upper ci limit around breakdate

# create dataframes
bdates_df <- as.data.frame(bdates_output, xy=FALSE, na.rm=FALSE)
lowerci_df <- as.data.frame(lowerci_output, xy=FALSE, na.rm=FALSE)
upperci_df <- as.data.frame(upperci_output, xy=FALSE, na.rm=FALSE)

#Only include rows that don't have only NA values to get rid of NO data pixels
bdates_df <- bdates_df[rowSums(is.na(bdates_df)) != (ncol(bdates_df)), ]
lowerci_df <- lowerci_df[rowSums(is.na(lowerci_df)) != (ncol(lowerci_df)), ]
lowerci_df <- lowerci_df[rowSums(is.na(lowerci_df)) != (ncol(lowerci_df)), ]

#Extract dates
x <- bdates_df[!is.na(bdates_df)]
y <- bdates_df[!is.na(bdates_df)]
z <- bdates_df[!is.na(bdates_df)]

#Create new data frame for results
d <- as.data.frame(x)

df_mean_NDVI$Date <- as.Date(df_mean_NDVI$Date)
brewer.pal("BrBG",n=11)
# "#543005" "#8C510A" "#BF812D" "#DFC27D" "#F6E8C3" "#F5F5F5" "#C7EAE5" "#80CDC1" "#35978F" "#01665E" "#003C30"
pal<-c("#01665E","#8C510A","#DFC27D","#80CDC1")
# pal<-c("black","grey85","#575756","#9B9B9B")
# Plot the data from 1999-2019
p1<-ggplot(data=df_mean_NDVI,aes(x=Date,
                             y=Mean,
                             color=Satellite,
                             fill=Satellite),na.rm=TRUE) +
  # geom_ribbon(aes(ymin= Mean - Std, ymax= Mean + Std),alpha=0.4) +
  geom_line(size=0.24) +
  geom_point(aes(x=Date,y=Mean),alpha=0.4,size=0.5) +
  # coord_cartesian(expand=FALSE) +
  ylim(0.025,0.225) +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
  labs(y="NDVI") +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        plot.title=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_text(size=9,color="black"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=11,color="black"),
        legend.position=c(0.99,0.05),legend.direction="horizontal",legend.justification=c("right","bottom"),
        legend.margin=margin(3,3,3,3),legend.background=element_rect(fill="grey95"),
        legend.title=element_text(size=9,color="black",face="bold"), legend.text=element_text(size=9,color="black")) +
  scale_color_manual(values=pal)

# Density plot
d$x<-date_decimal(d$x)
d$x<-as.Date(d$x,origin="1970-01-01")
p2<-ggplot(d, aes(x=x)) +  geom_density(aes(y = ..density..),color="#575756") + 
  labs(y="Breakpoint density") +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
  scale_y_continuous(labels = scales::comma) +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        axis.text.x=element_blank(),axis.text.y=element_text(size=9,color="black"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=11,color="black")) 


setwd(workdir_Prec)

monthly_precipitation_Ounila <- readRDS("monthly_precipitation_Ounila.Rda")
monthly_precipitation_Ounila$date<-as.Date(monthly_precipitation_Ounila$date)

setwd(workdir_Fig)

p3<-ggplot(data=monthly_precipitation_Ounila,aes(y=precipitation,x=date)) +
  geom_col(show.legend=FALSE,fill="#575756",size=1) +
  labs(x="Year",y="Precipitation (mm)") +
  scale_x_date(breaks=date_breaks("2 year"),labels=date_format("%Y"), limits=c(as.Date("1984-04-19","%Y-%m-%d"), as.Date("2020-01-01","%Y-%m-%d"))) +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        axis.text.x=element_text(size=9,color="black"),axis.text.y=element_text(size=9,color="black"),
        axis.title.x=element_text(size=11,color="black"),axis.title.y=element_text(size=11,color="black")) 

monthly_precipitation_year <- cbind(monthly_precipitation_Ounila$year,monthly_precipitation_Ounila$month,monthly_precipitation_Ounila$precipitation)
monthly_precipitation_year <- as.data.frame(monthly_precipitation_year)
monthly_precipitation_year <- as.precintcon.monthly(monthly_precipitation_year)

RAI_monthly<- rai(monthly_precipitation_year,granularity='m')
RAI_monthly$sign <- RAI_monthly$rai>0

RAI_monthly$date <- monthly_precipitation_Ounila$date
RAI_monthly$date <-as.Date(RAI_monthly$date)

p3<-ggplot(data=RAI_monthly,aes(y=rai,x=date,fill=sign)) +
  geom_col(show.legend=FALSE) +
  # ylim(-3.5,7.5) +
  scale_x_date(breaks=date_breaks("2 years"),labels=date_format("%Y"),limits=c(min(RAI_monthly$date),max(RAI_monthly$date))) +
  labs(x="Year",y="RAI") +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        axis.text.x=element_text(size=9,color="black"),axis.text.y=element_text(size=9,color="black"),
        axis.title.x=element_text(size=11,color="black"),axis.title.y=element_text(size=11,color="black")) +
  scale_fill_manual(values =c("#8C510A","#01665E"))

png(filename="NDVI_bpcount_P.png", family="Calibri",width=2031, height=2100, res=300)

grid.arrange(p1,p2,p3)

dev.off()

# Total number of breakpoints ---------------------------------------------
png(filename="BPNumb.png", family="Calibri",width=2200, height=1600, res=300)
par(mfrow = c(1,1), mar=c(4,4,4,10)+0.1)
plot(output$BPNumb,col=purple_pal,zlim=c(0,8),asp=1,legend=FALSE,main="Number of breakpoints per pixel",xlab="Longitude",ylab="Latitude")
# Map with legend next to it
legend(x='right', legend =  c(0:8),
       fill = purple_pal,bty='n', cex=1, xpd=NA,inset=-0.25)
dev.off()

output<-addLayer(output,output$BPNumb)

output$BPNumb.2[output$BPNumb.2==0] <- NA

output<-addLayer(output,output$Trend_noBP)

output$Trend_noBP.2[output$Trend_noBP.2<0] <- -9
output$Trend_noBP.2[output$Trend_noBP.2>0] <- 10
output$Trend_noBP.2[output$Trend_noBP.2==-9] <- 9

output$significance <- output$Trend_noBP_lowerCI*output$Trend_noBP_upperCI
output$significance[output$significance<0] <- NA
output$Trend_noBP.2 <- mask(output$Trend_noBP.2,output$significance,maskvalue=NA,updatevalue=11)

output$BPNumb.2 <- merge(output$BPNumb.2,output$Trend_noBP.2)

# output$BPNumb.2 can have values between 1 and 11 (1-8 is nr of breakpoints, 9 is decreasing trend, 10 is increasing trend, 11 is no trend)
png(filename="BPNumb_and_stable_trends.png", family="Calibri",width=2031, height=1330, res=300)
par(mfrow = c(1,1), mar=c(3,3,2,7)+0.2)
pal<-c(paste0(brewer.pal("BuPu",n=8)),"#DFC27D","#80CDC1","#9B9B9B")
# mask out cropland 
output$BPNumb.2 <- mask(output$BPNumb.2,Land_Use_Ounila,maskvalue=4)
plot(output$BPNumb.2,col=pal,zlim=c(1,11),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Number of breakpoints per pixel",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
# Map with legend next to it
legend(x='right', legend = c(1:8,"no bp, negative trend","no bp, postive trend","no bp, no trend"),
       fill = pal,bty='n', cex=1, xpd=NA,inset=-0.45)
dev.off()

# Number of breakpoints per land use type
bpnumb_grassland <- mask(output$BPNumb.2,Land_Use_Ounila,maskvalue=3,inverse=TRUE)
bpnumb_sparse_vegetation <- mask(output$BPNumb.2,Land_Use_Ounila,maskvalue=6,inverse=TRUE)
bpnumb_bare <- mask(output$BPNumb.2,Land_Use_Ounila,maskvalue=7,inverse=TRUE)

png(filename="BPNumb_per_LU_type.png", family="Calibri",width=2031, height=1600, res=300)
par(mfrow = c(2,2), mar=c(3.2,3,1.2,0))
pal<-c(paste0(brewer.pal("BuPu",n=8)),"#DFC27D","#80CDC1","#9B9B9B")
plot(bpnumb_grassland,col=pal,zlim=c(1,11),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Grassland",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb_sparse_vegetation,col=pal,zlim=c(1,11),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Sparse vegetation",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb_bare,col=pal,zlim=c(1,11),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Bare soil",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend = c(1:8,"no bp, negative trend","no bp, postive trend","no bp, no trend"),
       fill = pal,bty='n',cex=1, xpd=NA,inset=c(-2,0))
dev.off()

# Number of positive breakpoints ------------------------------------------
png(filename="BPNumb_positive_total.png", family="Calibri",width=2000, height=1500, res=300)
par(mar=c(3,3,1,0)+0.20)
output$BPNumb_pos<-output$BPNumbType1+output$BPNumbType2+output$BPNumbType6+output$BPNumbType7
# mask out cropland
output$BPNumb_pos <- mask(output$BPNumb_pos,Land_Use_Ounila,maskvalue=4)
plot(output$BPNumb_pos,col=pos_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Total number of positive breakpoints",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:6),fill = pos_pal,bty='n', cex=1, xpd=NA,inset=-0.12)
dev.off()

# number of breakpoints per positive breakpoint category
png(filename="BPNumb_positive.png", family="Calibri",width=2031, height=1550, res=300)
par(mfrow = c(2,2), mar=c(3,3,1,1)+0.20)
# 1 interrupted increase
# values     : 0, 6  (min, max)
# mask out cropland
output$BPNumbType1 <- mask(output$BPNumbType1,Land_Use_Ounila,maskvalue=4)
# 2	A significant increase after the breakpoint and no significant trend before the breakpoint
# values     : 0, 3  (min, max)
output$BPNumbType2 <- mask(output$BPNumbType2,Land_Use_Ounila,maskvalue=4)
# 6	A significant decrease before the breakpoint and no significant trend after the breakpoint
# values     : 0, 3  (min, max)
output$BPNumbType6 <- mask(output$BPNumbType6,Land_Use_Ounila,maskvalue=4)
# 7	A positive reversal: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
# values     : 0, 3  (min, max)
output$BPNumbType7 <- mask(output$BPNumbType7,Land_Use_Ounila,maskvalue=4)

plot(output$BPNumbType1,col=pos_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Interrupted increase",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:6),fill = pos_pal,bty='n', cex=1, xpd=NA,inset=-0.20)
plot(output$BPNumbType2,col=pos_pal[1:4],zlim=c(0,3),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Increase after no trend",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:3),fill = pos_pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.20)
plot(output$BPNumbType6,col=pos_pal[1:4],zlim=c(0,3),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("No trend after decrease",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:3),fill = pos_pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.20)
plot(output$BPNumbType7,col=pos_pal[1:4],zlim=c(0,3),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Positive reversal",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:3),fill = pos_pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.20)
dev.off()

# Number of breakpoints per land use type
bpnumb_grassland <- mask(output$BPNumb_pos,Land_Use_Ounila,maskvalue=3,inverse=TRUE)
bpnumb_sparse_vegetation <- mask(output$BPNumb_pos,Land_Use_Ounila,maskvalue=6,inverse=TRUE)
bpnumb_bare <- mask(output$BPNumb_pos,Land_Use_Ounila,maskvalue=7,inverse=TRUE)

png(filename="BPNumb_pos_per_LU_type.png", family="Calibri",width=2031, height=1600, res=300)
par(mfrow = c(2,2), mar=c(3.2,3,1.2,0))
pal<-pos_pal
plot(bpnumb_grassland,col=pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Grassland",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb_sparse_vegetation,col=pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Sparse vegetation",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb_bare,col=pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Bare soil",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend = c(0:6),fill = pos_pal,
       bty='n',cex=1, xpd=NA,inset=c(-2,0))
dev.off()


# Number of negative breakpoints ------------------------------------------
png(filename="BPNumb_negative_total.png", family="Calibri",width=2031, height=1500, res=300)
par(mar=c(3,3,1,0)+0.20)
output$BPNumb_neg<-output$BPNumbType4+output$BPNumbType5+output$BPNumbType3+output$BPNumbType8
# mask out cropland
output$BPNumb_neg <- mask(output$BPNumb_neg,Land_Use_Ounila,maskvalue=4)
plot(output$BPNumb_neg,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Total number of negative breakpoints",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:6),fill = neg_pal,bty='n', cex=1, xpd=NA,inset=-0.12)
dev.off()

# nr of breakpoints per negative breakpoint category
png(filename="BPNumb_negative.png", family="Calibri",width=2031, height=1550, res=300)
par(mfrow = c(2,2), mar=c(3,3,1,1)+0.20)
# 4	An interrupted decrease: both the trend before and after the breakpoint are negative and significant
# values     : 0, 6  (min, max)
output$BPNumbType4 <- mask(output$BPNumbType4,Land_Use_Ounila,maskvalue=4)
# 5	A significant decrease after the breakpoint and no significant trend before the breakpoint
# values     : 0, 3  (min, max)
output$BPNumbType5 <- mask(output$BPNumbType5,Land_Use_Ounila,maskvalue=4)
# 3	A significant increase before the breakpoint and no significant trend after the breakpoint
# values     : 0, 3  (min, max)
output$BPNumbType3 <- mask(output$BPNumbType3,Land_Use_Ounila,maskvalue=4)
# 8 negative reversal: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
# values     : 0, 4  (min, max)
output$BPNumbType8 <- mask(output$BPNumbType8,Land_Use_Ounila,maskvalue=4)

plot(output$BPNumbType4,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Interrupted decrease",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:6),fill = neg_pal,bty='n', cex=1, xpd=NA,inset=-0.20)
plot(output$BPNumbType5,col=neg_pal[1:4],zlim=c(0,3),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Decrease after no trend",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:3),fill = neg_pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.20)
plot(output$BPNumbType3,col=neg_pal[1:4],zlim=c(0,3),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("No trend after increase",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:3),fill = neg_pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.20)
plot(output$BPNumbType8,col=neg_pal[1:5],zlim=c(0,3),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Negative reversal",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c(0:4),fill = neg_pal[1:5],bty='n', cex=1, xpd=NA,inset=-0.20)
dev.off()

# Number of breakpoints per land use type
bpnumb_grassland <- mask(output$BPNumb_neg,Land_Use_Ounila,maskvalue=3,inverse=TRUE)
bpnumb_sparse_vegetation <- mask(output$BPNumb_neg,Land_Use_Ounila,maskvalue=6,inverse=TRUE)
bpnumb_bare <- mask(output$BPNumb_neg,Land_Use_Ounila,maskvalue=7,inverse=TRUE)

png(filename="BPNumb_neg_per_LU_type.png", family="Calibri",width=2031, height=1600, res=300)
par(mfrow = c(2,2), mar=c(3.2,3,1.2,0))
pal<-neg_pal
plot(bpnumb_grassland,col=pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Grassland",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb_sparse_vegetation,col=pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Sparse vegetation",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb_bare,col=pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("Bare soil",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend = c(0:6),fill = neg_pal,
       bty='n',cex=1, xpd=NA,inset=c(-2,0))
dev.off()

# Correlation between positive and negative breakpoints -------------------
# check whether the total breakpoint count follows a normal distribution 
qqnorm(df_bpoints$BPNumb)

# skewness:
e1071::skewness(df_bpoints$BPNumb)
# -0.1428425
# kurtosis
e1071::kurtosis(df_bpoints$BPNumb)
# -1.088179

# jarcque-bera test: 
normtest::ajb.norm.test(x=df_bpoints$BPNumb)
# Adjusted Jarque-Bera test for normality
# 
# data:  df_bpoints$BPNumb
# AJB = 42608, p-value < 2.2e-16

# conclusion: the data is not normally distributed 

# quickly check whether a log(x+1) or sqrt distribution improve normality:
df_bpoints$BPNumb_log_plus_1<-log(df_bpoints$BPNumb+1)
e1071::kurtosis(df_bpoints$BPNumb_log_plus_1)
# 0.4170488
e1071::skewness(df_bpoints$BPNumb_log_plus_1)
# -0.9244926
normtest::ajb.norm.test(x=df_bpoints$BPNumb_log_plus_1)
# Adjusted Jarque-Bera test for normality
# 
# data:  df_bpoints$BPNumb_log_plus_1
# AJB = 120939, p-value < 2.2e-16
df_bpoints$BPNumb_sqrt<-sqrt(df_bpoints$BPNumb)
e1071::skewness(df_bpoints$BPNumb_sqrt)
# -0.8408911
e1071::kurtosis(df_bpoints$BPNumb_sqrt)
# 0.5793668
normtest::ajb.norm.test(x=df_bpoints$BPNumb_sqrt)
# Adjusted Jarque-Bera test for normality
# 
# data:  df_bpoints$BPNumb_sqrt
# AJB = 106511, p-value < 2.2e-16
# this is also not the case, so: use a parametric tests

# test assumptions for a parametric test:

# calculate Tjostheims correlation coefficient
bpoints.cor_total_pos<-cor.spatial(df_bpoints$BPNumb,df_bpoints$BPNumb_pos,coords=df_bpoints[,1:2])
# [1] 0.3684128
# attr(,"variance")
# [1] 1.237774e-06
bpoints.cor_total_neg<-cor.spatial(df_bpoints$BPNumb,df_bpoints$BPNumb_neg,coords=df_bpoints[,1:2])
# [1] 0.4557195
# attr(,"variance")
# [1] 1.237774e-06
bpoints.cor_pos_neg<-cor.spatial(df_bpoints$BPNumb_pos,df_bpoints$BPNumb_neg,coords=df_bpoints[,1:2])
# [1] 0.1846345
# attr(,"variance")
# [1] 1.237774e-06

df_bpoints_cor<-subset(df_bpoints,BPNumb>0)
bpoints.cor_pos_neg<-cor.spatial(df_bpoints_cor$BPNumb_pos,df_bpoints_cor$BPNumb_neg,coords=df_bpoints_cor[,1:2])

bpnumb.cor<-rasterCorrelation(output$BPNumb_neg,output$BPNumb_pos,type="spearman",s=9)

png(filename="correlation.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)

plot(bpnumb.cor,breaks=seq(-1,1,0.25),col=(brewer.pal("BrBG",n=8)),zlim=c(-1,1),asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("Correlation between number of positive and negative breakpoints",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
plot(bpnumb.cor,legend.only=TRUE,breaks=seq(-1,1,0.25),col=(brewer.pal("BrBG",n=8)),legend.shrink=1,
     axis.args=list(at=c(-1.0,-0.5,0.0,0.5,1.0),
                    labels=c(-1.0,-0.5,0.0,0.5,1.0), 
                    cex.axis=0.8),
     legend.args=list(text='Spearman correlation coefficient', side=4, font=2, line=2.5, cex=1.1))

dev.off()

# try this correlatiobn again but then only with BPNumb>0
# acutally the only really interesting correlation is that between the negative and positive breakpoints
# its evident that total breakpoints is correlated with the two
# so only check negative and positive when a breakpoint occurs (so BPNumb>0, because when BPNumb==0 the nr of negative and positive breakpoints is ofcourse the same)

### now write text -- also include the calculated correlation!! 
## then: just put in the mf small plots (total/pos/neg bpoints; reflect on usefulness and relation to lu type maybe plots with pos, neg and lu.. or just look at the lu or see if overlay works this time with the imrpoved lu)
# correlation matrix postive/negative/total nr of breakpoints per lu class

# zoomed in plots of positive/negative/total(? met lu overlay?)

# Drought breakpoint typology ---------------------------------------------

output$nodbp <- mask(output$DBP,output$DBP,maskvalue=1,updatevalue=NA)
output$nodbp[output$nodbp==0]<-9
output$DBPType <- merge(output$DBPType,output$nodbp)
output$DBPType <- mask(output$DBPType,Land_Use_Ounila,maskvalue=4)

png(filename="DBPTypology.png", family="Calibri",width=2031, height=1270, res=300)
par(mfrow = c(1,1), mar=c(3,3,2,9)+0.2)

plot(output$DBPType,legend=FALSE,col=bp_pal,zlim=c(0,9),asp=1,xaxt='n',yaxt='n')
title("Typology of drought breakpoints",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c("no breakpoint","breakpoint, no trend",
                             "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
                             "interrupted decrease","decrease after no trend","no trend after increase","negative reversal"),
       fill = bp_pal_vis, bty='n', cex=1, xpd=NA,inset=-0.55)

dev.off()

names = c("no breakpoint","breakpoint, no trend",
          "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
          "interrupted decrease","decrease after no trend","no trend after increase","negative reversal")

masked_output<- mask(output,Land_Use_Ounila,maskvalue=4)
df_output <- as.data.frame(masked_output, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no flood breakpoint
df_output_no_dbp <- subset(df_output,DBP==0) # 478991 pixels
nrow(df_output_no_dbp)

# observations with flood breakpoint
df_output_dbp <- subset(df_output,DBP==1) # 328900 pixels
nrow(df_output_dbp)

# total number of observations
nobs <- nrow(df_output_no_dbp)+nrow(df_output_dbp) # 807891 pixels
perc_no_dbp <- nrow(df_output_no_dbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_dbp_type0 <- subset(df_output_dbp,DBPType==0) # 588 pixels
nrow(df_output_dbp_type0)
perc_dbp_type0 <- nrow(df_output_dbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_dbp_type1 <- subset(df_output_dbp,DBPType==1) # 588 pixels
nrow(df_output_dbp_type1)
perc_dbp_type1 <- nrow(df_output_dbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_dbp_type2 <- subset(df_output_dbp,DBPType==2) # 588 pixels
nrow(df_output_dbp_type2)
perc_dbp_type2 <- nrow(df_output_dbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_dbp_type3 <- subset(df_output_dbp,DBPType==3) # 588 pixels
nrow(df_output_dbp_type3)
perc_dbp_type3 <- nrow(df_output_dbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_dbp_type4 <- subset(df_output_dbp,DBPType==4) # 588 pixels
nrow(df_output_dbp_type4)
perc_dbp_type4 <- nrow(df_output_dbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_dbp_type5 <- subset(df_output_dbp,DBPType==5) # 588 pixels
nrow(df_output_dbp_type5)
perc_dbp_type5 <- nrow(df_output_dbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_dbp_type6 <- subset(df_output_dbp,DBPType==6) # 588 pixels
nrow(df_output_dbp_type6)
perc_dbp_type6 <- nrow(df_output_dbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_dbp_type7 <- subset(df_output_dbp,DBPType==7) # 588 pixels
nrow(df_output_dbp_type7)
perc_dbp_type7 <- nrow(df_output_dbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_dbp_type8 <- subset(df_output_dbp,DBPType==8) # 588 pixels
nrow(df_output_dbp_type8)
perc_dbp_type8 <- nrow(df_output_dbp_type8)/nobs

perc_no_dbp+perc_dbp_type0+perc_dbp_type1+perc_dbp_type2+perc_dbp_type3+perc_dbp_type4+perc_dbp_type5+perc_dbp_type6+perc_dbp_type7+perc_dbp_type8

percentages = c(perc_no_dbp,perc_dbp_type0,
                perc_dbp_type1,perc_dbp_type2,perc_dbp_type6,perc_dbp_type7,
                perc_dbp_type4,perc_dbp_type5,perc_dbp_type3,perc_dbp_type8)*100

df_dbp <- data.frame(names,percentages)
df_dbp$names <- factor(df_dbp$names,levels=names)

# For now just plot the barplot sepaerately for convencience
ggplot(data=df_dbp,aes(x=0,y=percentages)) +
  geom_col(aes(fill=names),show.legend=FALSE) +
  annotate(geom="text",x=c(0.01,0.01,0.01,0.01,0.01,0.01),y = c(90,74,68.5,62.5,42,15),label = c("22.3","4.5","6.8","5.3","35.4","17.2"), color = c("black","white","white","white","white","black")) +
  scale_fill_manual(values = bp_pal_vis) + 
  scale_y_continuous(breaks=c(0,25,50,75,100),labels=c("100","75","50","25","0"),trans="reverse") +
  ylab("Percentage of pixels") +
  scale_x_continuous(expand=c(0,0)) +
  coord_flip() +
  theme(axis.title.x=element_text(color="black"),
        axis.text.x=element_text(color="black"),
        axis.ticks.x=element_line(color="black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(filename="Perc_bar_DBP.png", family="Calibri", width = 6.8, height = 1, dpi=300)


# crop output to the 3 major land use classes
output_grassland <- mask(masked_output,Land_Use_Ounila,inverse=TRUE,maskvalue=3,updatevalue=NA)
output_sparse_vegetation <- mask(masked_output,Land_Use_Ounila,inverse=TRUE,maskvalue=6,updatevalue=NA)
output_bare <- mask(masked_output,Land_Use_Ounila,inverse=TRUE,maskvalue=7,updatevalue=NA)

# plot a percentage bar per land use class:
# grassland:
df_output_grassland <- as.data.frame(output_grassland, xy=TRUE, na.rm=FALSE) # 1725748 pixels

### make sure that totally empty values arent included -- is ok i think becaues DBP==1/0
# observations with no flood breakpoint
df_output_grassland_no_dbp <- subset(df_output_grassland,DBP==0) # 61194 pixels
nrow(df_output_grassland_no_dbp)

# observations with flood breakpoint
df_output_grassland_dbp <- subset(df_output_grassland,DBP==1) # 219687 pixels
nrow(df_output_grassland_dbp)

# total number of observations
nobs <- nrow(df_output_grassland_no_dbp)+nrow(df_output_grassland_dbp)
perc_no_dbp <- nrow(df_output_grassland_no_dbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_grassland_dbp_type0 <- subset(df_output_grassland_dbp,DBPType==0) # 1073 pixels
nrow(df_output_grassland_dbp_type0)
perc_dbp_type0 <- nrow(df_output_grassland_dbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_grassland_dbp_type1 <- subset(df_output_grassland_dbp,DBPType==1) # 12131 pixels
nrow(df_output_grassland_dbp_type1)
perc_dbp_type1 <- nrow(df_output_grassland_dbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_grassland_dbp_type2 <- subset(df_output_grassland_dbp,DBPType==2) # 14283 pixels
nrow(df_output_grassland_dbp_type2)
perc_dbp_type2 <- nrow(df_output_grassland_dbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_grassland_dbp_type3 <- subset(df_output_grassland_dbp,DBPType==3) # 4679 pixels
nrow(df_output_grassland_dbp_type3)
perc_dbp_type3 <- nrow(df_output_grassland_dbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_grassland_dbp_type4 <- subset(df_output_grassland_dbp,DBPType==4) # 37061 pixels
nrow(df_output_grassland_dbp_type4)
perc_dbp_type4 <- nrow(df_output_grassland_dbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_grassland_dbp_type5 <- subset(df_output_grassland_dbp,DBPType==5) # 4077 pixels
nrow(df_output_grassland_dbp_type5)
perc_dbp_type5 <- nrow(df_output_grassland_dbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_grassland_dbp_type6 <- subset(df_output_grassland_dbp,DBPType==6) # 11341 pixels
nrow(df_output_grassland_dbp_type6)
perc_dbp_type6 <- nrow(df_output_grassland_dbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_grassland_dbp_type7 <- subset(df_output_grassland_dbp,DBPType==7) # 128613 pixels
nrow(df_output_grassland_dbp_type7)
perc_dbp_type7 <- nrow(df_output_grassland_dbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_grassland_dbp_type8 <- subset(df_output_grassland_dbp,DBPType==8) # 6429 pixels
nrow(df_output_grassland_dbp_type8)
perc_dbp_type8 <- nrow(df_output_grassland_dbp_type8)/nobs

perc_no_dbp+perc_dbp_type0+perc_dbp_type1+perc_dbp_type2+perc_dbp_type3+perc_dbp_type4+perc_dbp_type5+perc_dbp_type6+perc_dbp_type7+perc_dbp_type8

percentages = c(perc_no_dbp,perc_dbp_type0,
                perc_dbp_type1,perc_dbp_type2,perc_dbp_type6,perc_dbp_type7,
                perc_dbp_type4,perc_dbp_type5,perc_dbp_type3,perc_dbp_type8)*100

df_dbp_grassland <- data.frame(names,percentages)
df_dbp_grassland$names <- factor(df_dbp_grassland$names,levels=names)
df_dbp_grassland$landuse <- "grassland"

# sparse_vegetation:
df_output_sparse_vegetation <- as.data.frame(output_sparse_vegetation, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no flood breakpoint
df_output_sparse_vegetation_no_dbp <- subset(df_output_sparse_vegetation,DBP==0) # 47218 pixels
nrow(df_output_sparse_vegetation_no_dbp)

# observations with flood breakpoint
df_output_sparse_vegetation_dbp <- subset(df_output_sparse_vegetation,DBP==1) # 111906 pixels
nrow(df_output_sparse_vegetation_dbp)

# total number of observations
nobs <- nrow(df_output_sparse_vegetation_no_dbp)+nrow(df_output_sparse_vegetation_dbp) 
perc_no_dbp <- nrow(df_output_sparse_vegetation_no_dbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_sparse_vegetation_dbp_type0 <- subset(df_output_sparse_vegetation_dbp,DBPType==0) # 980 pixels
nrow(df_output_sparse_vegetation_dbp_type0)
perc_dbp_type0 <- nrow(df_output_sparse_vegetation_dbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_sparse_vegetation_dbp_type1 <- subset(df_output_sparse_vegetation_dbp,DBPType==1) # 5463 pixels
nrow(df_output_sparse_vegetation_dbp_type1)
perc_dbp_type1 <- nrow(df_output_sparse_vegetation_dbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_sparse_vegetation_dbp_type2 <- subset(df_output_sparse_vegetation_dbp,DBPType==2) # 8372 pixels
nrow(df_output_sparse_vegetation_dbp_type2)
perc_dbp_type2 <- nrow(df_output_sparse_vegetation_dbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_sparse_vegetation_dbp_type3 <- subset(df_output_sparse_vegetation_dbp,DBPType==3) # 2086 pixels
nrow(df_output_sparse_vegetation_dbp_type3)
perc_dbp_type3 <- nrow(df_output_sparse_vegetation_dbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_sparse_vegetation_dbp_type4 <- subset(df_output_sparse_vegetation_dbp,DBPType==4) # 26305 pixels
nrow(df_output_sparse_vegetation_dbp_type4)
perc_dbp_type4 <- nrow(df_output_sparse_vegetation_dbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_sparse_vegetation_dbp_type5 <- subset(df_output_sparse_vegetation_dbp,DBPType==5) # 2072 pixels
nrow(df_output_sparse_vegetation_dbp_type5)
perc_dbp_type5 <- nrow(df_output_sparse_vegetation_dbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_sparse_vegetation_dbp_type6 <- subset(df_output_sparse_vegetation_dbp,DBPType==6) # 5840 pixels
nrow(df_output_sparse_vegetation_dbp_type6)
perc_dbp_type6 <- nrow(df_output_sparse_vegetation_dbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_sparse_vegetation_dbp_type7 <- subset(df_output_sparse_vegetation_dbp,DBPType==7) # 58981 pixels
nrow(df_output_sparse_vegetation_dbp_type7)
perc_dbp_type7 <- nrow(df_output_sparse_vegetation_dbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_sparse_vegetation_dbp_type8 <- subset(df_output_sparse_vegetation_dbp,DBPType==8) # 1807 pixels
nrow(df_output_sparse_vegetation_dbp_type8)
perc_dbp_type8 <- nrow(df_output_sparse_vegetation_dbp_type8)/nobs

perc_no_dbp+perc_dbp_type0+perc_dbp_type1+perc_dbp_type2+perc_dbp_type3+perc_dbp_type4+perc_dbp_type5+perc_dbp_type6+perc_dbp_type7+perc_dbp_type8

percentages = c(perc_no_dbp,perc_dbp_type0,
                perc_dbp_type1,perc_dbp_type2,perc_dbp_type6,perc_dbp_type7,
                perc_dbp_type4,perc_dbp_type5,perc_dbp_type3,perc_dbp_type8)*100

df_dbp_sparse_vegetation <- data.frame(names,percentages)
df_dbp_sparse_vegetation$names <- factor(df_dbp_sparse_vegetation$names,levels=names)
df_dbp_sparse_vegetation$landuse <- "sparse vegetation"

# bare:
df_output_bare <- as.data.frame(output_bare, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no flood breakpoint
df_output_bare_no_dbp <- subset(df_output_bare,DBP==0) # 65711 pixels
nrow(df_output_bare_no_dbp)

# observations with flood breakpoint
df_output_bare_dbp <- subset(df_output_bare,DBP==1) # 277619 pixels
nrow(df_output_bare_dbp)

# total number of observations
nobs <- nrow(df_output_bare_no_dbp)+nrow(df_output_bare_dbp) 
perc_no_dbp <- nrow(df_output_bare_no_dbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_bare_dbp_type0 <- subset(df_output_bare_dbp,DBPType==0) # 9224 pixels
nrow(df_output_bare_dbp_type0)
perc_dbp_type0 <- nrow(df_output_bare_dbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_bare_dbp_type1 <- subset(df_output_bare_dbp,DBPType==1) # 17638 pixels
nrow(df_output_bare_dbp_type1)
perc_dbp_type1 <- nrow(df_output_bare_dbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_bare_dbp_type2 <- subset(df_output_bare_dbp,DBPType==2) # 30283 pixels
nrow(df_output_bare_dbp_type2)
perc_dbp_type2 <- nrow(df_output_bare_dbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_bare_dbp_type3 <- subset(df_output_bare_dbp,DBPType==3) # 4423 pixels
nrow(df_output_bare_dbp_type3)
perc_dbp_type3 <- nrow(df_output_bare_dbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_bare_dbp_type4 <- subset(df_output_bare_dbp,DBPType==4) # 71538 pixels
nrow(df_output_bare_dbp_type4)
perc_dbp_type4 <- nrow(df_output_bare_dbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_bare_dbp_type5 <- subset(df_output_bare_dbp,DBPType==5) # 20468 pixels
nrow(df_output_bare_dbp_type5)
perc_dbp_type5 <- nrow(df_output_bare_dbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_bare_dbp_type6 <- subset(df_output_bare_dbp,DBPType==6) # 24160 pixels
nrow(df_output_bare_dbp_type6)
perc_dbp_type6 <- nrow(df_output_bare_dbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_bare_dbp_type7 <- subset(df_output_bare_dbp,DBPType==7) # 90376 pixels
nrow(df_output_bare_dbp_type7)
perc_dbp_type7 <- nrow(df_output_bare_dbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_bare_dbp_type8 <- subset(df_output_bare_dbp,DBPType==8) # 9509 pixels
nrow(df_output_bare_dbp_type8)
perc_dbp_type8 <- nrow(df_output_bare_dbp_type8)/nobs

perc_no_dbp+perc_dbp_type0+perc_dbp_type1+perc_dbp_type2+perc_dbp_type3+perc_dbp_type4+perc_dbp_type5+perc_dbp_type6+perc_dbp_type7+perc_dbp_type8

percentages = c(perc_no_dbp,perc_dbp_type0,
                perc_dbp_type1,perc_dbp_type2,perc_dbp_type6,perc_dbp_type7,
                perc_dbp_type4,perc_dbp_type5,perc_dbp_type3,perc_dbp_type8)*100

df_dbp_bare <- data.frame(names,percentages)
df_dbp_bare$names <- factor(df_dbp_bare$names,levels=names)
df_dbp_bare$landuse <- "bare"

df_dbp <- rbind(df_dbp_bare,df_dbp_grassland,df_dbp_sparse_vegetation)
df_dbp$landuse <- factor(df_dbp$landuse)

setwd(workdir_Fig)

ggplot(data=df_dbp,aes(x=landuse,y=percentages)) +
  geom_col(aes(fill=names),show.legend=FALSE) +
  annotate(geom="text",x=c(3,2,1,1,1,3,2,1,3,2,1),y = c(85,90,90,69,61,45,45,45,12,12,21),label = c("29.7","21.8","19.1","8.8","7.0","37.1","45.8","26.3","16.5","13.2","20.8"), color =c("black","black","black","white","white","white","white","white","black","black","black")) +
  scale_fill_manual(values = bp_pal_vis) + 
  scale_y_continuous(breaks=c(0,25,50,75,100),labels=c("100","75","50","25","0"),trans="reverse") +
  ylab("Percentage of pixels") +
  scale_x_discrete() +
  coord_flip() +
  theme(axis.title.x=element_text(color="black"),
        axis.text.x=element_text(color="black"),
        axis.ticks.x=element_line(color="black"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color="black",size=9),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(filename="Perc_bar_DBP_landuse.png", family="Calibri", width = 6.8, height = 2.3, dpi=300)

# Flood breakpoint typology ---------------------------------------------

output$nofbp <- mask(output$FBP,output$FBP,maskvalue=1,updatevalue=NA)
output$nofbp[output$nofbp==0]<-9
output$FBPType <- merge(output$FBPType,output$nofbp)
output$FBPType <- mask(output$FBPType,Land_Use_Ounila,maskvalue=4)

png(filename="FBPtypology.png", family="Calibri",width=2031, height=1270, res=300)
par(mfrow = c(1,1), mar=c(3,3,2,9)+0.2)

plot(output$FBPType,legend=FALSE,col=bp_pal,zlim=c(0,9),asp=1,xaxt='n',yaxt='n')
title("Typology of flood breakpoints",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")
legend(x='right', legend = c("no breakpoint","breakpoint, no trend",
                             "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
                             "interrupted decrease","decrease after no trend","no trend after increase","negative reversal"),
       fill = bp_pal_vis, bty='n', cex=1, xpd=NA,inset=-0.55)

dev.off()

names = c("no breakpoint","breakpoint, no trend",
          "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
          "interrupted decrease","decrease after no trend","no trend after increase","negative reversal")

df_output <- as.data.frame(masked_output, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no flood breakpoint
df_output_no_fbp <- subset(df_output,FBP==0) # 478991 pixels
nrow(df_output_no_fbp)

# observations with flood breakpoint
df_output_fbp <- subset(df_output,FBP==1) # 328900 pixels
nrow(df_output_fbp)

# total number of observations
nobs <- nrow(df_output_no_fbp)+nrow(df_output_fbp) # 807891 pixels
perc_no_fbp <- nrow(df_output_no_fbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_fbp_type0 <- subset(df_output_fbp,FBPType==0) # 588 pixels
nrow(df_output_fbp_type0)
perc_fbp_type0 <- nrow(df_output_fbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_fbp_type1 <- subset(df_output_fbp,FBPType==1) # 588 pixels
nrow(df_output_fbp_type1)
perc_fbp_type1 <- nrow(df_output_fbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_fbp_type2 <- subset(df_output_fbp,FBPType==2) # 588 pixels
nrow(df_output_fbp_type2)
perc_fbp_type2 <- nrow(df_output_fbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_fbp_type3 <- subset(df_output_fbp,FBPType==3) # 588 pixels
nrow(df_output_fbp_type3)
perc_fbp_type3 <- nrow(df_output_fbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_fbp_type4 <- subset(df_output_fbp,FBPType==4) # 588 pixels
nrow(df_output_fbp_type4)
perc_fbp_type4 <- nrow(df_output_fbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_fbp_type5 <- subset(df_output_fbp,FBPType==5) # 588 pixels
nrow(df_output_fbp_type5)
perc_fbp_type5 <- nrow(df_output_fbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_fbp_type6 <- subset(df_output_fbp,FBPType==6) # 588 pixels
nrow(df_output_fbp_type6)
perc_fbp_type6 <- nrow(df_output_fbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_fbp_type7 <- subset(df_output_fbp,FBPType==7) # 588 pixels
nrow(df_output_fbp_type7)
perc_fbp_type7 <- nrow(df_output_fbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_fbp_type8 <- subset(df_output_fbp,FBPType==8) # 588 pixels
nrow(df_output_fbp_type8)
perc_fbp_type8 <- nrow(df_output_fbp_type8)/nobs

perc_no_fbp+perc_fbp_type0+perc_fbp_type1+perc_fbp_type2+perc_fbp_type3+perc_fbp_type4+perc_fbp_type5+perc_fbp_type6+perc_fbp_type7+perc_fbp_type8

percentages = c(perc_no_fbp,perc_fbp_type0,
                perc_fbp_type1,perc_fbp_type2,perc_fbp_type6,perc_fbp_type7,
                perc_fbp_type4,perc_fbp_type5,perc_fbp_type3,perc_fbp_type8)*100

df_fbp <- data.frame(names,percentages)
df_fbp$names <- factor(df_fbp$names,levels=names)

# For now just plot the barplot sepaerately for convencience
ggplot(data=df_fbp,aes(x=0,y=percentages)) +
  geom_col(aes(fill=names),show.legend=FALSE) +
  annotate(geom="text",x=c(0.01,0.01),y = c(70,20),label = c("59.0","35.5"), color = c("black","black")) +
  scale_fill_manual(values = bp_pal_vis) + 
  scale_y_continuous(breaks=c(0,25,50,75,100),labels=c("100","75","50","25","0"),trans="reverse") +
  ylab("Percentage of pixels") +
  scale_x_continuous(expand=c(0,0)) +
  coord_flip() +
  theme(axis.title.x=element_text(color="black"),
        axis.text.x=element_text(color="black"),
        axis.ticks.x=element_line(color="black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(filename="Perc_bar_FBP.png", family="Calibri", width = 6.8, height = 1, dpi=300)

# plot a percentage bar per land use class:
# grassland:
df_output_grassland <- as.data.frame(output_grassland, xy=TRUE, na.rm=FALSE) # 1725748 pixels

### make sure that totally empty values arent included -- is ok i think becaues fbp==1/0
# observations with no flood breakpoint
df_output_grassland_no_fbp <- subset(df_output_grassland,FBP==0) # 188263 pixels
nrow(df_output_grassland_no_fbp)

# observations with flood breakpoint
df_output_grassland_fbp <- subset(df_output_grassland,FBP==1) # 92618 pixels
nrow(df_output_grassland_fbp)

# total number of observations
nobs <- nrow(df_output_grassland_no_fbp)+nrow(df_output_grassland_fbp) 
perc_no_fbp <- nrow(df_output_grassland_no_fbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_grassland_fbp_type0 <- subset(df_output_grassland_fbp,FBPType==0) # 339 pixels
nrow(df_output_grassland_fbp_type0)
perc_fbp_type0 <- nrow(df_output_grassland_fbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_grassland_fbp_type1 <- subset(df_output_grassland_fbp,FBPType==1) # 2773 pixels
nrow(df_output_grassland_fbp_type1)
perc_fbp_type1 <- nrow(df_output_grassland_fbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_grassland_fbp_type2 <- subset(df_output_grassland_fbp,FBPType==2) # 399 pixels
nrow(df_output_grassland_fbp_type2)
perc_fbp_type2 <- nrow(df_output_grassland_fbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_grassland_fbp_type3 <- subset(df_output_grassland_fbp,FBPType==3) # 1221 pixels
nrow(df_output_grassland_fbp_type3)
perc_fbp_type3 <- nrow(df_output_grassland_fbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_grassland_fbp_type4 <- subset(df_output_grassland_fbp,FBPType==4) # 76327 pixels
nrow(df_output_grassland_fbp_type4)
perc_fbp_type4 <- nrow(df_output_grassland_fbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_grassland_fbp_type5 <- subset(df_output_grassland_fbp,FBPType==5) # 2787 pixels
nrow(df_output_grassland_fbp_type5)
perc_fbp_type5 <- nrow(df_output_grassland_fbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_grassland_fbp_type6 <- subset(df_output_grassland_fbp,FBPType==6) # 4004 pixels
nrow(df_output_grassland_fbp_type6)
perc_fbp_type6 <- nrow(df_output_grassland_fbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_grassland_fbp_type7 <- subset(df_output_grassland_fbp,FBPType==7) # 1549 pixels
nrow(df_output_grassland_fbp_type7)
perc_fbp_type7 <- nrow(df_output_grassland_fbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_grassland_fbp_type8 <- subset(df_output_grassland_fbp,FBPType==8) # 3219 pixels
nrow(df_output_grassland_fbp_type8)
perc_fbp_type8 <- nrow(df_output_grassland_fbp_type8)/nobs

perc_no_fbp+perc_fbp_type0+perc_fbp_type1+perc_fbp_type2+perc_fbp_type3+perc_fbp_type4+perc_fbp_type5+perc_fbp_type6+perc_fbp_type7+perc_fbp_type8

percentages = c(perc_no_fbp,perc_fbp_type0,
                perc_fbp_type1,perc_fbp_type2,perc_fbp_type6,perc_fbp_type7,
                perc_fbp_type4,perc_fbp_type5,perc_fbp_type3,perc_fbp_type8)*100

df_fbp_grassland <- data.frame(names,percentages)
df_fbp_grassland$names <- factor(df_fbp_grassland$names,levels=names)
df_fbp_grassland$landuse <- "grassland"

# sparse_vegetation:
df_output_sparse_vegetation <- as.data.frame(output_sparse_vegetation, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no flood breakpoint
df_output_sparse_vegetation_no_fbp <- subset(df_output_sparse_vegetation,FBP==0) # 107052 pixels
nrow(df_output_sparse_vegetation_no_fbp)

# observations with flood breakpoint
df_output_sparse_vegetation_fbp <- subset(df_output_sparse_vegetation,FBP==1) # 52072 pixels
nrow(df_output_sparse_vegetation_fbp)

# total number of observations
nobs <- nrow(df_output_sparse_vegetation_no_fbp)+nrow(df_output_sparse_vegetation_fbp) 
perc_no_fbp <- nrow(df_output_sparse_vegetation_no_fbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_sparse_vegetation_fbp_type0 <- subset(df_output_sparse_vegetation_fbp,FBPType==0) # 88 pixels
nrow(df_output_sparse_vegetation_fbp_type0)
perc_fbp_type0 <- nrow(df_output_sparse_vegetation_fbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_sparse_vegetation_fbp_type1 <- subset(df_output_sparse_vegetation_fbp,FBPType==1) # 2221 pixels
nrow(df_output_sparse_vegetation_fbp_type1)
perc_fbp_type1 <- nrow(df_output_sparse_vegetation_fbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_sparse_vegetation_fbp_type2 <- subset(df_output_sparse_vegetation_fbp,FBPType==2) # 136 pixels
nrow(df_output_sparse_vegetation_fbp_type2)
perc_fbp_type2 <- nrow(df_output_sparse_vegetation_fbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_sparse_vegetation_fbp_type3 <- subset(df_output_sparse_vegetation_fbp,FBPType==3) # 776 pixels
nrow(df_output_sparse_vegetation_fbp_type3)
perc_fbp_type3 <- nrow(df_output_sparse_vegetation_fbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_sparse_vegetation_fbp_type4 <- subset(df_output_sparse_vegetation_fbp,FBPType==4) # 45218 pixels
nrow(df_output_sparse_vegetation_fbp_type4)
perc_fbp_type4 <- nrow(df_output_sparse_vegetation_fbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_sparse_vegetation_fbp_type5 <- subset(df_output_sparse_vegetation_fbp,FBPType==5) # 1049 pixels
nrow(df_output_sparse_vegetation_fbp_type5)
perc_fbp_type5 <- nrow(df_output_sparse_vegetation_fbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_sparse_vegetation_fbp_type6 <- subset(df_output_sparse_vegetation_fbp,FBPType==6) # 1288 pixels
nrow(df_output_sparse_vegetation_fbp_type6)
perc_fbp_type6 <- nrow(df_output_sparse_vegetation_fbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_sparse_vegetation_fbp_type7 <- subset(df_output_sparse_vegetation_fbp,FBPType==7) # 532 pixels
nrow(df_output_sparse_vegetation_fbp_type7)
perc_fbp_type7 <- nrow(df_output_sparse_vegetation_fbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_sparse_vegetation_fbp_type8 <- subset(df_output_sparse_vegetation_fbp,FBPType==8) # 764 pixels
nrow(df_output_sparse_vegetation_fbp_type8)
perc_fbp_type8 <- nrow(df_output_sparse_vegetation_fbp_type8)/nobs

perc_no_fbp+perc_fbp_type0+perc_fbp_type1+perc_fbp_type2+perc_fbp_type3+perc_fbp_type4+perc_fbp_type5+perc_fbp_type6+perc_fbp_type7+perc_fbp_type8

percentages = c(perc_no_fbp,perc_fbp_type0,
                perc_fbp_type1,perc_fbp_type2,perc_fbp_type6,perc_fbp_type7,
                perc_fbp_type4,perc_fbp_type5,perc_fbp_type3,perc_fbp_type8)*100

df_fbp_sparse_vegetation <- data.frame(names,percentages)
df_fbp_sparse_vegetation$names <- factor(df_fbp_sparse_vegetation$names,levels=names)
df_fbp_sparse_vegetation$landuse <- "sparse vegetation"

# bare:
df_output_bare <- as.data.frame(output_bare, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# observations with no flood breakpoint
df_output_bare_no_fbp <- subset(df_output_bare,FBP==0) # 166335 pixels
nrow(df_output_bare_no_fbp)

# observations with flood breakpoint
df_output_bare_fbp <- subset(df_output_bare,FBP==1) # 176995 pixels
nrow(df_output_bare_fbp)

# total number of observations
nobs <- nrow(df_output_bare_no_fbp)+nrow(df_output_bare_fbp) 
perc_no_fbp <- nrow(df_output_bare_no_fbp)/nobs

# percentage breakpoints following Type 0: 'no significant trends' typology: both the trend before and after the breakpoint are non-significant
df_output_bare_fbp_type0 <- subset(df_output_bare_fbp,FBPType==0) # 38 pixels
nrow(df_output_bare_fbp_type0)
perc_fbp_type0 <- nrow(df_output_bare_fbp_type0)/nobs

# percentage breakpoints following Type 1: 'interrupted increase' typology: both the trend before and after the breakpoint are positive and significant
df_output_bare_fbp_type1 <- subset(df_output_bare_fbp,FBPType==1) # 1602 pixels
nrow(df_output_bare_fbp_type1)
perc_fbp_type1 <- nrow(df_output_bare_fbp_type1)/nobs

# percentage breakpoints following Type 2: 'increase after no trend' typology: significant increase after the breakpoint and no significant trend before the breakpoint
df_output_bare_fbp_type2 <- subset(df_output_bare_fbp,FBPType==2) # 174 pixels
nrow(df_output_bare_fbp_type2)
perc_fbp_type2 <- nrow(df_output_bare_fbp_type2)/nobs

# percentage breakpoints following Type 3: 'no trend after increase' typology: significant increase before the breakpoint and no significant trend after the breakpoint
df_output_bare_fbp_type3 <- subset(df_output_bare_fbp,FBPType==3) # 1199 pixels
nrow(df_output_bare_fbp_type3)
perc_fbp_type3 <- nrow(df_output_bare_fbp_type3)/nobs

# percentage breakpoints following Type 4: 'interrupted decrease ' typology: both the trend before and after the breakpoint are negative and significant
df_output_bare_fbp_type4 <- subset(df_output_bare_fbp,FBPType==4) # 157545 pixels
nrow(df_output_bare_fbp_type4)
perc_fbp_type4 <- nrow(df_output_bare_fbp_type4)/nobs

# percentage breakpoints following Type 5: 'decrease after no trend' typology: significant decrease after the breakpoint and no significant trend before the breakpoint 
df_output_bare_fbp_type5 <- subset(df_output_bare_fbp,FBPType==5) # 3319 pixels
nrow(df_output_bare_fbp_type5)
perc_fbp_type5 <- nrow(df_output_bare_fbp_type5)/nobs

# percentage breakpoints following Type 6: 'no trend after decrease' typology: significant decrease before the breakpoint and no significant trend after the breakpoint
df_output_bare_fbp_type6 <- subset(df_output_bare_fbp,FBPType==6) # 2206 pixels
nrow(df_output_bare_fbp_type6)
perc_fbp_type6 <- nrow(df_output_bare_fbp_type6)/nobs

# percentage breakpoints following Type 7: 'positive reversal' typology: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
df_output_bare_fbp_type7 <- subset(df_output_bare_fbp,FBPType==7) # 779 pixels
nrow(df_output_bare_fbp_type7)
perc_fbp_type7 <- nrow(df_output_bare_fbp_type7)/nobs

# percentage breakpoints following Type 8: 'negative reversal' typology: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.
df_output_bare_fbp_type8 <- subset(df_output_bare_fbp,FBPType==8) # 10133 pixels
nrow(df_output_bare_fbp_type8)
perc_fbp_type8 <- nrow(df_output_bare_fbp_type8)/nobs

perc_no_fbp+perc_fbp_type0+perc_fbp_type1+perc_fbp_type2+perc_fbp_type3+perc_fbp_type4+perc_fbp_type5+perc_fbp_type6+perc_fbp_type7+perc_fbp_type8

percentages = c(perc_no_fbp,perc_fbp_type0,
                perc_fbp_type1,perc_fbp_type2,perc_fbp_type6,perc_fbp_type7,
                perc_fbp_type4,perc_fbp_type5,perc_fbp_type3,perc_fbp_type8)*100

df_fbp_bare <- data.frame(names,percentages)
df_fbp_bare$names <- factor(df_fbp_bare$names,levels=names)
df_fbp_bare$landuse <- "bare"

df_fbp <- rbind(df_fbp_bare,df_fbp_grassland,df_fbp_sparse_vegetation)
df_fbp$landuse <- factor(df_fbp$landuse)

setwd(workdir_Fig)

ggplot(data=df_fbp,aes(x=landuse,y=percentages)) +
  geom_col(aes(fill=names),show.legend=FALSE) +
  annotate(geom="text",x=c(3,2,1,3,2,1),y = c(65,65,75,18,18,28),label = c("67.3","67.0","48.4","28.4","27.2","45.9"), color ="black") +
  scale_fill_manual(values = bp_pal_vis) + 
  scale_y_continuous(breaks=c(0,25,50,75,100),labels=c("100","75","50","25","0"),trans="reverse") +
  ylab("Percentage of pixels") +
  scale_x_discrete() +
  coord_flip() +
  theme(axis.title.x=element_text(color="black"),
        axis.text.x=element_text(color="black"),
        axis.ticks.x=element_line(color="black"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color="black",size=9),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(filename="Perc_bar_FBP_landuse.png", family="Calibri", width = 6.8, height = 2.3, dpi=300)

# Change in NDVI ----------------------------------------------------------

png(filename="NDVIchange.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45)
pal<-brewer.pal("BrBG",n=11)[3:10]

# output$NDVIchange <- output$End.NDVI-output$Initial.NDVI
output$NDVIchange <- (output$End.NDVI-output$Initial.NDVI)
plot(output$NDVIchange,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n')
title("Change in NDVI",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="NDVIchange1.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(-0.4,-0.3,-0.2,-0.1,0.00,0.1,0.2,0.3,0.4,0.5)
pal<-c(paste0(brewer.pal("BrBG",n=11)[2:5]),paste0(brewer.pal("BrBG",n=11)[7:11]))

# output$NDVIchange <- output$End.NDVI-output$Initial.NDVI
output$NDVIchange <- (output$End.NDVI-output$Initial.NDVI)
plot(output$NDVIchange,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n')
title("Change in NDVI",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

# potential areas of interest 
# A:
# class      : Extent 
# xmin       : -7.198685 
# xmax       : -7.161118 
# ymin       : 31.29425 
# ymax       : 31.32102 

map_aoi_A <- get_googlemap(center=c(-(7.198685+7.161118)/2,(31.29425+31.32102)/2),
                           zoom = 14,
                           size=c(640,640),
                           scale=2,
                           format="png8",
                           maptype="satellite",
                           language = "en-EN",
                           filename="satellite_aoi_A.png",
                           color="color")
ggmap(map_aoi_A) +
  geom_rect(aes(xmin=-7.198685,xmax=-7.161118,ymin=31.29425,ymax=31.32102),color="red",fill="transparent",size=2) +
  coord_fixed(ratio=1) 
ggsave("map_aoi_A.png", family="Calibri", width = 3.4, height = 3, dpi=300)

# B:
# class      : Extent 
# xmin       : -7.122735 
# xmax       : -7.071827 
# ymin       : 31.28253 
# ymax       : 31.3298
map_aoi_B <- get_googlemap(center=c(-(7.122735+7.071827)/2,(31.28253+31.3298)/2),
                           zoom = 13,
                           size=c(640,640),
                           scale=2,
                           format="png8",
                           maptype="satellite",
                           language = "en-EN",
                           filename="satellite_aoi_B.png",
                           color="color")
ggmap(map_aoi_B) +
  geom_rect(aes(xmin=-7.122735,xmax=-7.071827,ymin=31.28253,ymax=31.3298),color="red",fill="transparent",size=2)  +
  coord_fixed(ratio=1) 
ggsave("map_aoi_B.png", family="Calibri", width = 3.4, height = 3, dpi=300)
# C:
# class      : Extent 
# xmin       : -7.221368 
# xmax       : -7.18137 
# ymin       : 31.27071 
# ymax       : 31.29299 
map_aoi_C <- get_googlemap(center=c(-(7.221368+7.18137)/2,(31.27071+31.29299)/2),
                           zoom = 14,
                           size=c(640,640),
                           scale=2,
                           format="png8",
                           maptype="satellite",
                           language = "en-EN",
                           filename="satellite_aoi_C.png",
                           color="color")
ggmap(map_aoi_C) +
  geom_rect(aes(xmin=-7.221368,xmax=-7.18137,ymin=31.27071,ymax=31.29299),color="red",fill="transparent",size=2) +
  coord_fixed(ratio=1) 
ggsave("map_aoi_C.png", family="Calibri", width = 3.4, height = 3, dpi=300)
# D:
# class      : Extent 
# xmin       : -7.275003 
# xmax       : -7.243186 
# ymin       : 31.24208 
# ymax       : 31.26753 
map_aoi_D <- get_googlemap(center=c(-(7.275003+7.243186)/2,(31.24208+31.26753)/2),
                           zoom = 14,
                           size=c(640,640),
                           scale=2,
                           format="png8",
                           maptype="satellite",
                           language = "en-EN",
                           filename="satellite_aoi_D.png",
                           color="color")
ggmap(map_aoi_D) +
  geom_rect(aes(xmin=-7.275003,xmax=-7.243186,ymin=31.24208,ymax=31.26753),color="red",fill="transparent",size=2) +
  coord_fixed(ratio=1) 
ggsave("map_aoi_D.png", family="Calibri", width = 3.4, height = 3, dpi=300)
# E:
# class      : Extent 
# xmin       : -7.14728 
# xmax       : -7.083645 
# ymin       : 31.23253 
# ymax       : 31.26299
map_aoi_E <- get_googlemap(center=c(-(7.14728+7.083645)/2,(31.23253+31.26299)/2),
                           zoom = 13,
                           size=c(640,640),
                           scale=2,
                           format="png8",
                           maptype="satellite",
                           language = "en-EN",
                           filename="satellite_aoi_E.png",
                           color="color")
ggmap(map_aoi_E) +
  geom_rect(aes(xmin=-7.14728,xmax=-7.083645,ymin=31.23253,ymax=31.26299),color="red",fill="transparent",size=2) +
  coord_fixed(ratio=1) 
ggsave("map_aoi_E.png", family="Calibri", width = 3.4, height = 3, dpi=300)


png(filename="NDVIchangepercent.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(-150,-50,-30,-10,10,30,50,150,250,350)
pal<-brewer.pal("BrBG",n=11)[3:11]

output$NDVIchangepercent <- ((output$End.NDVI-output$Initial.NDVI)/output$Initial.NDVI)*100
plot(output$NDVIchangepercent,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
# title("Change in NDVI",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
plot(output$NDVIchangepercent,legend.only=TRUE,breaks=brks,col=pal,legend.shrink=1,
     axis.args=list(at=c(-150,-47.5,-27,-7,13,32.5,52.5,151,251,350),
                    labels=brks, 
                    cex.axis=0.8),
     legend.args=list(text='Percentage change in NDVI (%)', side=4, font=2, line=2.5, cex=1.1))


plot(extent(-7.198685 ,-7.161118 ,31.29425,31.32102), col="red",add=T)
text(x=-7.198685+0.005 ,y=31.32102 ,labels="A",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.122735,-7.071827,31.28253,31.3298), col="red",add=T)
text(x=-7.122735+0.005,y=31.3298,labels="B",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.221368,-7.18137 ,31.27071,31.29299), col="red",add=T)
text(x=-7.221368+0.005,y=31.29299,labels="C",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.275003 ,-7.243186,31.24208,31.26753), col="red",add=T)
text(x=-7.275003+0.005 ,y=31.26753 ,labels="D",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.14728 ,-7.083645,31.23253,31.26299), col="red",add=T)
text(x=-7.14728+0.005 ,y=31.26299 ,labels="E",cex=0.8,pos=1,col="red",font=2,offset=0.2)

dev.off()

# NDVI change with croplands masked out 

png(filename="NDVIchangepercent_withoutLU.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(-150,-50,-30,-10,10,30,50,150,250,350)
pal<-brewer.pal("BrBG",n=11)[3:11]

output$NDVIchangepercent <- ((output$End.NDVI-output$Initial.NDVI)/output$Initial.NDVI)*100
output$NDVIchangepercent_withoutLU <- mask(output$NDVIchangepercent,Land_Use_Ounila,maskvalue=4)

plot(output$NDVIchangepercent_withoutLU,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE,colNA="black")
# title("Change in NDVI",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
plot(output$NDVIchangepercent_withoutLU,legend.only=TRUE,breaks=brks,col=pal,legend.shrink=1,
     axis.args=list(at=c(-150,-47.5,-27,-7,13,32.5,52.5,151,251,350),
                    labels=brks, 
                    cex.axis=0.8),
     legend.args=list(text='Percentage change in NDVI (%)', side=4, font=2, line=2.5, cex=1.1))

plot(extent(-7.198685 ,-7.161118 ,31.29425,31.32102), col="red",add=T)
text(x=-7.198685+0.005 ,y=31.32102 ,labels="A",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.122735,-7.071827,31.28253,31.3298), col="red",add=T)
text(x=-7.122735+0.005,y=31.3298,labels="B",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.221368,-7.18137 ,31.27071,31.29299), col="red",add=T)
text(x=-7.221368+0.005,y=31.29299,labels="C",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.275003 ,-7.243186,31.24208,31.26753), col="red",add=T)
text(x=-7.275003+0.005 ,y=31.26753 ,labels="D",cex=0.8,pos=1,col="red",font=2,offset=0.2)

plot(extent(-7.14728 ,-7.083645,31.23253,31.26299), col="red",add=T)
text(x=-7.14728+0.005 ,y=31.26299 ,labels="E",cex=0.8,pos=1,col="red",font=2,offset=0.2)

dev.off()

# zoomed_ndvi -------------------------------------------------------------

setwd(workdir_LandUse)
Land_Use_Ounila<-raster("Land_Use_Ounila.tif")
LU_Ounila_A <- crop(Land_Use_Ounila,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
LU_Ounila_B <- crop(Land_Use_Ounila,extent(-7.122735,-7.071827,31.28253,31.3298))
LU_Ounila_C <- crop(Land_Use_Ounila,extent(-7.221368,-7.18137 ,31.27071,31.29299))
LU_Ounila_D <- crop(Land_Use_Ounila,extent(-7.275003 ,-7.243186,31.24208,31.26753))
LU_Ounila_E <- crop(Land_Use_Ounila,extent(-7.14728 ,-7.083645,31.23253,31.26299))

setwd(workdir_Fig)

NDVIchangepercent_A <- crop(output$NDVIchangepercent,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
NDVIchangepercent_B <- crop(output$NDVIchangepercent,extent(-7.122735,-7.071827,31.28253,31.3298))
NDVIchangepercent_C <- crop(output$NDVIchangepercent,extent(-7.221368,-7.18137 ,31.27071,31.29299))
NDVIchangepercent_D <- crop(output$NDVIchangepercent,extent(-7.275003 ,-7.243186,31.24208,31.26753))
NDVIchangepercent_E <- crop(output$NDVIchangepercent,extent(-7.14728 ,-7.083645,31.23253,31.26299))

png(filename="NDVIchange_LU_A.png", family="Calibri",width=1320, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(-50,-30,-10,10,30,50,150)
pal<-brewer.pal("BrBG",n=11)[4:9]

plot(NDVIchangepercent_A,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("A.",line=0.5)
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255))

plot(LU_Ounila_A, legend = FALSE, col = Land_Use_Palette, zlim=c(3,8), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.198685+0.001,31.29425+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="NDVIchange_LU_B.png", family="Calibri",width=1100, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(-150,-50,-30,-10,10,30,50,150)
pal<-brewer.pal("BrBG",n=11)[3:9]

plot(NDVIchangepercent_B,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("B.",line=0.5)
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                      rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255))

plot(LU_Ounila_B, legend = FALSE, col = Land_Use_Palette, zlim=c(1,8), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.122735+0.002,31.28253+0.003),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="NDVIchange_LU_C.png", family="Calibri",width=1540, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(-30,-10,10,30,50,150,250)
pal<-brewer.pal("BrBG",n=11)[5:10]

plot(NDVIchangepercent_C,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("C.",line=0.5)
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255))

plot(LU_Ounila_C, legend = FALSE, col = Land_Use_Palette, zlim=c(3,8), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.221368+0.002,31.27071+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="NDVIchange_LU_D.png", family="Calibri",width=1200, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(-150,-50,-30,-10,10,30,50,150,250)
pal<-brewer.pal("BrBG",n=11)[3:10]

plot(NDVIchangepercent_D,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("D.",line=0.5)
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                      rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255))

plot(LU_Ounila_D, legend = FALSE, col = Land_Use_Palette, zlim=c(1,8), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.275003+0.002,31.24208+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="NDVIchange_LU_E.png", family="Calibri",width=1630, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(-30,-10,10,30,50,150)
pal<-brewer.pal("BrBG",n=11)[5:9]

plot(NDVIchangepercent_E,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("E.",line=0.5)
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255),
                      rgb(255,255,255, max=255),rgb(0,70,200, max=255))

plot(LU_Ounila_E, legend = FALSE, col = Land_Use_Palette, zlim=c(3,10), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.14728+0.002,31.23253+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="zoomed_in_NDVIchange_legend.png", family="Calibri",width=2031, height=1100, res=300)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
brks=c(-150,-50,-30,-10,10,30,50,150,250)
pal<-brewer.pal("BrBG",n=11)[3:10]
plot(NDVIchangepercent_A,legend.only=TRUE,breaks=brks,col=pal,legend.shrink=1,
     smallplot=c(0.15,.18, 0.2,0.9),
     axis.args=list(at=c(-150,-47.5,-27,-7,13,32.5,52.5,151,250),
                    labels=brks, 
                    cex.axis=1,
                    tck=-0.5),
     legend.args=list(text='Percentage change in NDVI (%)', side=4, font=2, line=2.5, cex=1.1))
dev.off()

png(filename="zoomed_in_NDVIchange_LUlegend.png", family="Calibri",width=2031, height=1100, res=300)
Land_Use_Palette_Vis  <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                           rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                           rgb(255,235,175, max=255),
                           rgb(255,245,215, max=255),rgb(195,20,0, max=255),
                           rgb(0,70,200, max=255))

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend =  c("Tree cover areas","Shrub cover areas",
                                      "Grassland","Cropland",
                                      "Lichen Mosses / Sparse Vegetation",
                                      "Bare areas", "Built up areas",
                                      "Open water"),
       fill = Land_Use_Palette_Vis,bty='n',cex=1, xpd=NA,inset=c(-2,0))
dev.off()

png(filename="NDVIchangepercent_LU_C.png", family="Calibri",width=1016, height=1100, res=300)
par(mfrow = c(2,1), mar=c(3,2.8,0.8,0))
brks=c(-150,-50,-30,-10,10,30,50,150)
pal<-brewer.pal("BrBG",n=11)[3:9]

plot(NDVIchangepercent_C,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
title("C.Tighza and degraded land",line=0.35,cex.main=0.85)
axis(2,seq(31.29,31.34,0.01),cex.axis=0.55,tck=-0.02)
axis(1,seq(-7.13,-7.06,0.01),cex.axis=0.55,tck=-0.02)
mtext("Latitude",side=2,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                      rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255))

plot(LU_Ounila_C, legend = FALSE, col = Land_Use_Palette, zlim=c(1,8), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.29,31.34,0.01),cex.axis=0.55,tck=-0.02)
axis(1,seq(-7.13,-7.06,0.01),cex.axis=0.55,tck=-0.02)
mtext("Latitude",side=2,line=2,cex=0.8)
mtext("Longitude",side=1,line=2,cex=0.8)

dev.off()


png(filename="NDVIchangepercent_LU_D.png", family="Calibri",width=1215, height=850, res=300)
par(mfrow = c(1,2), mar=c(3,2.6,0.8,0))
brks=c(-150,-50,-30,-10,10,30,50,150,250)
pal<-brewer.pal("BrBG",n=11)[3:10]

plot(NDVIchangepercent_D,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
mtext("D. Confluence of river channels", side = 3, font=2, line =-0.7, outer = TRUE,cex=0.85)

# title("",line=0.35)
axis(2,seq(31.03,31.11,0.01),cex.axis=0.55,tck=-0.02)
axis(1,seq(-7.17,-7.13,0.01),cex.axis=0.55,tck=-0.02)
mtext("Latitude",side=2,line=2,cex=0.8)
mtext("Longitude",side=1,line=2,cex=0.8)

Land_Use_Palette <- c(rgb(0,160,0, max=255),rgb(150,100,0, max=255),
                      rgb(255,180,0, max=255),rgb(255,255,100, max=255),
                      rgb(0,220,130, max=255),rgb(255,235,175, max=255),
                      rgb(255,245,215, max=255),rgb(195,20,0, max=255))

plot(LU_Ounila_D, legend = FALSE, col = Land_Use_Palette, zlim=c(1,8), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
axis(2,seq(31.03,31.11,0.01),cex.axis=0.55,tck=-0.02)
axis(1,seq(-7.17,-7.13,0.01),cex.axis=0.55,tck=-0.02)
mtext("Longitude",side=1,line=2,cex=0.8)

dev.off()


###

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
bty='n', cex=1, xpd=NA,inset=-0.85)

# DBPPreNDVI --------------------------------------------------------------

pal<-brewer.pal("BrBG",n=11)[7:11]
output$DBPPreNDVI <- mask(output$DBPPreNDVI,Land_Use_Ounila,maskvalue=4)

png(filename="DBPPreNDVI.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(0.00,0.15,0.30,0.45,0.60,0.75)

plot(output$DBPPreNDVI,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n')
title("Mean NDVI during 3 years before the drought breakpoint",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

# FBPPreNDVI --------------------------------------------------------------

pal<-brewer.pal("BrBG",n=11)[7:11]
output$FBPPreNDVI <- mask(output$FBPPreNDVI,Land_Use_Ounila,maskvalue=4)

png(filename="FBPPreNDVI.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(0.00,0.15,0.30,0.45,0.60,0.75)

plot(output$FBPPreNDVI,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n')
title("Mean NDVI during 3 years before the flood breakpoint",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

# Initial NDVI --------------------------------------------------------------

pal<-brewer.pal("BrBG",n=11)[7:10]

png(filename="Initial_NDVI.png", family="Calibri",width=2031, height=1500, res=300)
par(mfrow = c(1,1), mar=c(3,3,1,1)+0.2)
brks=c(0.00,0.15,0.30,0.45,0.60)

plot(output$Initial.NDVI,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n')
title("NDVI during first three years",line=0.5)
axis(2,seq(31.0,31.6,0.1),cex.axis=0.8)
axis(1,seq(-7.4,-7.0,0.1),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
raster::scalebar(d=10,xy=c(-7.39,31.05),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

# Plot of drought breakpoint typologies in areas of interest --------------
output$nodbp <- mask(output$DBP,output$DBP,maskvalue=1,updatevalue=NA)
output$nodbp[output$nodbp==0]<-9
output$DBPType <- merge(output$DBPType,output$nodbp)
output$DBPType <- mask(output$DBPType,Land_Use_Ounila,maskvalue=4)
output$DBPPreNDVI <- mask(output$DBPPreNDVI,Land_Use_Ounila,maskvalue=4)

DBPType_A <- crop(output$DBPType,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
DBPType_B <- crop(output$DBPType,extent(-7.122735,-7.071827,31.28253,31.3298))
DBPType_C <- crop(output$DBPType,extent(-7.221368,-7.18137 ,31.27071,31.29299))
DBPType_D <- crop(output$DBPType,extent(-7.275003 ,-7.243186,31.24208,31.26753))
DBPType_E <- crop(output$DBPType,extent(-7.14728 ,-7.083645,31.23253,31.26299))

DBPPreNDVI_A <- crop(output$DBPPreNDVI,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
DBPPreNDVI_B <- crop(output$DBPPreNDVI,extent(-7.122735,-7.071827,31.28253,31.3298))
DBPPreNDVI_C <- crop(output$DBPPreNDVI,extent(-7.221368,-7.18137 ,31.27071,31.29299))
DBPPreNDVI_D <- crop(output$DBPPreNDVI,extent(-7.275003 ,-7.243186,31.24208,31.26753))
DBPPreNDVI_E <- crop(output$DBPPreNDVI,extent(-7.14728 ,-7.083645,31.23253,31.26299))

png(filename="DBPType_PreDBPNDVI_A.png", family="Calibri",width=1320, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30)
pal<-brewer.pal("BrBG",n=11)[7:8]

plot(DBPType_A, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("A.",line=0.5)
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(DBPPreNDVI_A,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.198685+0.001,31.29425+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="DBPType_PreDBPNDVI_B.png", family="Calibri",width=1100, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(DBPType_B, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("B.",line=0.5)
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(DBPPreNDVI_B,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.122735+0.002,31.28253+0.003),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="DBPType_PreDBPNDVI_C.png", family="Calibri",width=1540, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(DBPType_C, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("C.",line=0.5)
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(DBPPreNDVI_C,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.221368+0.002,31.27071+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="DBPType_PreDBPNDVI_D.png", family="Calibri",width=1200, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(DBPType_D, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("D.",line=0.5)
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(DBPPreNDVI_D,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.275003+0.002,31.24208+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="DBPType_PreDBPNDVI_E.png", family="Calibri",width=1630, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(DBPType_E, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("E.",line=0.5)
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(DBPPreNDVI_E,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.14728+0.002,31.23253+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="zoomed_in_preDBPNDVI_legend.png", family="Calibri",width=2031, height=1100, res=300)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:10]
plot(DBPPreNDVI_A,legend.only=TRUE,breaks=brks,col=pal,legend.shrink=1,
     smallplot=c(0.25,0.28, 0.22,0.92),
     axis.args=list(at=brks,
                    labels=brks, 
                    cex.axis=1,
                    tck=-0.5),
     legend.args=list(text='Mean NDVI 3 years before breakpoint', side=4, font=2, line=2.5, cex=1.1))
dev.off()

png(filename="zoomed_in_bp_legend.png", family="Calibri",width=2031, height=1100, res=300)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend =  c("no breakpoint","breakpoint, no trend",
                                      "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
                                      "interrupted decrease","decrease after no trend","no trend after increase","negative reversal"),
       fill = bp_pal_vis,bty='n',cex=1, xpd=NA,inset=c(-2,0))

dev.off()


# Plot of Flood breakpoint typology in areas of interest ------------------
output$nofbp <- mask(output$FBP,output$FBP,maskvalue=1,updatevalue=NA)
output$nofbp[output$nofbp==0]<-9
output$FBPType <- merge(output$FBPType,output$nofbp)
output$FBPType <- mask(output$FBPType,Land_Use_Ounila,maskvalue=4)
output$FBPPreNDVI <- mask(output$FBPPreNDVI,Land_Use_Ounila,maskvalue=4)

FBPType_A <- crop(output$FBPType,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
FBPType_B <- crop(output$FBPType,extent(-7.122735,-7.071827,31.28253,31.3298))
FBPType_C <- crop(output$FBPType,extent(-7.221368,-7.18137 ,31.27071,31.29299))
FBPType_D <- crop(output$FBPType,extent(-7.275003 ,-7.243186,31.24208,31.26753))
FBPType_E <- crop(output$FBPType,extent(-7.14728 ,-7.083645,31.23253,31.26299))

FBPPreNDVI_A <- crop(output$FBPPreNDVI,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
FBPPreNDVI_B <- crop(output$FBPPreNDVI,extent(-7.122735,-7.071827,31.28253,31.3298))
FBPPreNDVI_C <- crop(output$FBPPreNDVI,extent(-7.221368,-7.18137 ,31.27071,31.29299))
FBPPreNDVI_D <- crop(output$FBPPreNDVI,extent(-7.275003 ,-7.243186,31.24208,31.26753))
FBPPreNDVI_E <- crop(output$FBPPreNDVI,extent(-7.14728 ,-7.083645,31.23253,31.26299))

png(filename="FBPType_PreFBPNDVI_A.png", family="Calibri",width=1320, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30)
pal<-brewer.pal("BrBG",n=11)[7:8]

plot(FBPType_A, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("A.",line=0.5)
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(FBPPreNDVI_A,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.198685+0.001,31.29425+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="FBPType_PreFBPNDVI_B.png", family="Calibri",width=1100, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(FBPType_B, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("B.",line=0.5)
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(FBPPreNDVI_B,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.122735+0.002,31.28253+0.003),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="FBPType_PreFBPNDVI_C.png", family="Calibri",width=1540, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(FBPType_C, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("C.",line=0.5)
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(FBPPreNDVI_C,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.221368+0.002,31.27071+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="FBPType_PreFBPNDVI_D.png", family="Calibri",width=1200, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30,0.45)
pal<-brewer.pal("BrBG",n=11)[7:9]

plot(FBPType_D, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("D.",line=0.5)
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(FBPPreNDVI_D,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.275003+0.002,31.24208+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="FBPType_PreFBPNDVI_E.png", family="Calibri",width=1630, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

brks=c(0.00,0.15,0.30)
pal<-brewer.pal("BrBG",n=11)[7:8]

plot(FBPType_E, legend = FALSE, col = bp_pal, zlim=c(0,9), asp=1,
     xaxt="n",yaxt="n",
     xlab="", ylab="",main="")
title("E.",line=0.5)
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(FBPPreNDVI_E,breaks=brks,col=pal,asp=1,xaxt='n',yaxt='n',legend=FALSE)
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.14728+0.002,31.23253+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="zoomed_in_preFBPNDVI_legend.png", family="Calibri",width=2031, height=1100, res=300)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
brks=c(0.00,0.15,0.30,0.45,0.60)
pal<-brewer.pal("BrBG",n=11)[7:10]
plot(FBPPreNDVI_A,legend.only=TRUE,breaks=brks,col=pal,legend.shrink=1,
     smallplot=c(0.25,0.28, 0.22,0.92),
     axis.args=list(at=brks,
                    labels=brks, 
                    cex.axis=1,
                    tck=-0.5),
     legend.args=list(text='Mean NDVI 3 years before breakpoint', side=4, font=2, line=2.5, cex=1.1))
dev.off()

png(filename="zoomed_in_bp_legend.png", family="Calibri",width=2031, height=1100, res=300)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend =  c("no breakpoint","breakpoint, no trend",
                                      "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
                                      "interrupted decrease","decrease after no trend","no trend after increase","negative reversal"),
       fill = bp_pal_vis,bty='n',cex=1, xpd=NA,inset=c(-2,0))

dev.off()


# Plot of positve and negative breakpoints in areas of interest -----------

output$BPNumb_pos<-output$BPNumbType1+output$BPNumbType2+output$BPNumbType6+output$BPNumbType7
output$BPNumb_neg<-output$BPNumbType4+output$BPNumbType5+output$BPNumbType3+output$BPNumbType8

output$BPNumb_pos <- mask(output$BPNumb_pos,Land_Use_Ounila,maskvalue=4)
output$BPNumb_neg <- mask(output$BPNumb_neg,Land_Use_Ounila,maskvalue=4)

BPNumb_pos_A <- crop(output$BPNumb_pos,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
BPNumb_pos_B <- crop(output$BPNumb_pos,extent(-7.122735,-7.071827,31.28253,31.3298))
BPNumb_pos_C <- crop(output$BPNumb_pos,extent(-7.221368,-7.18137 ,31.27071,31.29299))
BPNumb_pos_D <- crop(output$BPNumb_pos,extent(-7.275003 ,-7.243186,31.24208,31.26753))
BPNumb_pos_E <- crop(output$BPNumb_pos,extent(-7.14728 ,-7.083645,31.23253,31.26299))

BPNumb_neg_A <- crop(output$BPNumb_neg,extent(-7.198685 ,-7.161118 ,31.29425,31.32102))
BPNumb_neg_B <- crop(output$BPNumb_neg,extent(-7.122735,-7.071827,31.28253,31.3298))
BPNumb_neg_C <- crop(output$BPNumb_neg,extent(-7.221368,-7.18137 ,31.27071,31.29299))
BPNumb_neg_D <- crop(output$BPNumb_neg,extent(-7.275003 ,-7.243186,31.24208,31.26753))
BPNumb_neg_E <- crop(output$BPNumb_neg,extent(-7.14728 ,-7.083645,31.23253,31.26299))

png(filename="bpnumb_pos_neg_A.png", family="Calibri",width=1320, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

plot(BPNumb_pos_A,col=pos_pal[1:6],zlim=c(0,5),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("A.",line=0.5)
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(BPNumb_neg_A,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
axis(2,seq(31.29,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.20,-7.16,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.198685+0.001,31.29425+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="bpnumb_pos_neg_B.png", family="Calibri",width=1100, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

plot(BPNumb_pos_B,col=pos_pal[1:6],zlim=c(0,5),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("B.",line=0.5)
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(BPNumb_neg_B,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
axis(2,seq(31.28,31.33,0.01),cex.axis=0.8)
axis(1,seq(-7.13,-7.07,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.122735+0.002,31.28253+0.003),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="bpnumb_pos_neg_C.png", family="Calibri",width=1540, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

plot(BPNumb_pos_C,col=pos_pal[1:6],zlim=c(0,5),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("C.",line=0.5)
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(BPNumb_neg_C,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
axis(2,seq(31.27,31.30,0.01),cex.axis=0.8)
axis(1,seq(-7.24,-7.18,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.221368+0.002,31.27071+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="bpnumb_pos_neg_D.png", family="Calibri",width=1200, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

plot(BPNumb_pos_D,col=pos_pal[1:6],zlim=c(0,5),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("D.",line=0.5)
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(BPNumb_neg_D,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
axis(2,seq(31.24,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.28,-7.24,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.275003+0.002,31.24208+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="bpnumb_pos_neg_E.png", family="Calibri",width=1630, height=1700, res=300)
par(mfrow = c(2,1), mar=c(3,3,1,0)+0.2)

plot(BPNumb_pos_E,col=pos_pal[1:6],zlim=c(0,5),asp=1,legend=FALSE,xaxt='n',yaxt='n')
title("E.",line=0.5)
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2,cex=0.8)

plot(BPNumb_neg_E,col=neg_pal,zlim=c(0,6),asp=1,legend=FALSE,xaxt='n',yaxt='n')
axis(2,seq(31.23,31.27,0.01),cex.axis=0.8)
axis(1,seq(-7.15,-7.08,0.01),cex.axis=0.8)
mtext("Latitude",side=2,line=2.2,cex=0.8)
mtext("Longitude",side=1,line=2.2,cex=0.8)
raster::scalebar(d=1,xy=c(-7.14728+0.002,31.23253+0.002),type="bar",lonlat=TRUE,divs=2,cex=0.7,below="kilometres")

dev.off()

png(filename="zoomed_in_pos_legend.png", family="Calibri",width=2031, height=1100, res=300)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend = c(0:5),
       fill = pos_pal[1:6],bty='n',cex=1, xpd=NA,inset=c(-2,0))

dev.off()


png(filename="zoomed_in_neg_legend.png", family="Calibri",width=2031, height=1100, res=300)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x="center",ncol=1, legend = c(0:6),
       fill = neg_pal,bty='n',cex=1, xpd=NA,inset=c(-2,0))

dev.off()

