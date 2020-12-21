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
if(!require(gdalUtils)){install.packages("gdalUtils")}
if(!require(bfast)){install.packages("bfast")}
if(!require(bfastSpatial)){install.packages("bfastSpatial")}
if(!require(extrafont)){install.packages("extrafont")}
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
if(!require(ggplot2)){install.packages("ggplot2")}

# Load libraries 
library(raster) # package for raster manipulation
library(parallel)
library(strucchange)
library(MASS)
library(zoo)
library(rgdal)
library(gdalUtils)
library(bfast)
library(bfastSpatial)
library(extrafont)
library(RColorBrewer)
library(ggplot2)

# import fonts 

# Set working directories and load data  ----------------------------------

workdir_Fig = "G:/Thesis_data/Figures"
workdir_output = "G:/Thesis_data/output_resIndSpatial"
workdir_NDVI = "G:/Thesis_data/NDVI"

setwd(workdir_output)

output <- brick("output.grd")

setwd(workdir_NDVI)

NDVI <- brick("NDVI_final_OLS.grd")

setwd(workdir_Fig)

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



# Total number of breakpoints ---------------------------------------------

png(filename="BPNumb.png", family="Calibri",width=2200, height=1600, res=300)
par(mfrow = c(1,1), mar=c(4,4,4,10)+0.1)
plot(output$BPNumb,col=colorRampPalette(c("#FFFFFF","#6600CC"))(10)[2:10],zlim=c(0,8),legend=FALSE,main="Number of breakpoints per pixel",xlab="Longitude",ylab="Latitude")
# Map with legend next to it
legend(x='right', legend =  c(0:8),
       fill = colorRampPalette(c("#FFFFFF","#6600CC"))(10)[2:10],bty='n', cex=1, xpd=NA,inset=-0.25)
dev.off()
# Plot bar with percentages of total cells beneath it
# add another (white - for no bpoint (this is also registered in dbp (o or 1)))

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
png(filename="BPNumb_and_stable_trends.png", family="Calibri",width=2200, height=1600, res=300)
par(mfrow = c(1,1), mar=c(4,4,4,10)+0.1)
colorRampPalette(c("#FFFFFF","#6600CC"))(9)
pal<-c( "#EBDFF8", "#D8BFF2", "#C59FEB", "#B27FE5", "#9F5FDF", "#8C3FD8", "#791FD2", "#6600CC","#CC9445","#71BB76","#9B9B9B")
plot(output$BPNumb.2,col=pal,zlim=c(1,11),legend=FALSE,main="Number of breakpoints per pixel",xlab="Longitude",ylab="Latitude")
# Map with legend next to it
legend(x='right', legend = c(1:8,"no bp, negative trend","no bp, postive trend","no bp, no trend"),
       fill = pal,bty='n', cex=1, xpd=NA,inset=-0.55)
dev.off()

# Number of positive breakpoints ------------------------------------------

png(filename="BPNumb_positive.png", family="Calibri",width=2400, height=2000, res=300)
par(mfrow = c(2,2), mar=c(3,3,1,2)+0.20)
# 1 interrupted increase
# values     : 0, 6  (min, max)
# 2	A significant increase after the breakpoint and no significant trend before the breakpoint
# values     : 0, 3  (min, max)
# 6	A significant decrease before the breakpoint and no significant trend after the breakpoint
# values     : 0, 3  (min, max)
# 7	A positive reversal: trend before breakpoint was negative, trend after breakpoints is positive. Both trends are significant.
# values     : 0, 3  (min, max)
pal <- colorRampPalette(c("#FFFFFF","#0A2CF0"))(8)[2:8] # web colors
plot(output$BPNumbType1,col=pal,zlim=c(0,6),legend=FALSE)
title("Interrupted increase",line=0.5)
mtext("Latitude",side=2,line=2.2)
legend(x='right', legend = c(0:6),fill = pal,bty='n', cex=1, xpd=NA,inset=-0.15)
plot(output$BPNumbType2,col=pal[1:4],zlim=c(0,3),legend=FALSE)
title("Increase after no trend",line=0.5)
legend(x='right', legend = c(0:3),fill = pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.15)
plot(output$BPNumbType6,col=pal[1:4],zlim=c(0,3),legend=FALSE)
title("No trend after decrease",line=0.5)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
legend(x='right', legend = c(0:3),fill = pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.15)
plot(output$BPNumbType7,col=pal[1:4],zlim=c(0,3),legend=FALSE)
title("Positive reversal",line=0.5)
mtext("Longitude",side=1,line=2.2)
legend(x='right', legend = c(0:3),fill = pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.15)
dev.off()


# Number of negative breakpoints ------------------------------------------

png(filename="BPNumb_negative.png", family="Calibri",width=2400, height=2000, res=300)
par(mfrow = c(2,2), mar=c(3,3,1,2)+0.20)
# 4	An interrupted decrease: both the trend before and after the breakpoint are negative and significant
# values     : 0, 6  (min, max)
# 5	A significant decrease after the breakpoint and no significant trend before the breakpoint
# values     : 0, 3  (min, max)
# 3	A significant increase before the breakpoint and no significant trend after the breakpoint
# values     : 0, 3  (min, max)
# 8 negative reversal: trend before breakpoint was positive, trend after breakpoint is negative. Both trends are significant.\
# values     : 0, 4  (min, max)
pal <- colorRampPalette(c("#FFFFFF","#FF1398"))(8)[2:8]
plot(output$BPNumbType1,col=pal,zlim=c(0,6),legend=FALSE)
title("Interrupted decrease",line=0.5)
mtext("Latitude",side=2,line=2.2)
legend(x='right', legend = c(0:6),fill = pal,bty='n', cex=1, xpd=NA,inset=-0.15)
plot(output$BPNumbType2,col=pal[1:4],zlim=c(0,3),legend=FALSE)
title("Decrease after no trend",line=0.5)
legend(x='right', legend = c(0:3),fill = pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.15)
plot(output$BPNumbType6,col=pal[1:4],zlim=c(0,3),legend=FALSE)
title("No trend after increase",line=0.5)
mtext("Latitude",side=2,line=2.2)
mtext("Longitude",side=1,line=2.2)
legend(x='right', legend = c(0:3),fill = pal[1:4],bty='n', cex=1, xpd=NA,inset=-0.15)
plot(output$BPNumbType7,col=pal[1:5],zlim=c(0,3),legend=FALSE)
title("Negative reversal",line=0.5)
mtext("Longitude",side=1,line=2.2)
legend(x='right', legend = c(0:4),fill = pal[1:5],bty='n', cex=1, xpd=NA,inset=-0.15)
dev.off()

# Drought breakpoint typology ---------------------------------------------

png(filename="DBPTypology.png", family="Calibri",width=2400, height=1600, res=300)
par(mfrow = c(1,1), mar=c(5,4,4,16)+0.1)

# webcolors -- use the other pallete (see illustrator file) for print colors
bp_type_pal <- c("#575756",
                 "#81FA45","#36F58A","#E8190C",
                 "#FFC414","#FF5001","#0A98F2",
                 "#0A2CF0","#FF1398")

plot(output$DBPType,legend=FALSE,col=bp_type_pal,zlim=c(0,8),
     main="Typology of drought breakpoints",xlab="Longitude",ylab="Latitude")

bp_type_pal_vis <- c("#575756",
                     "#81FA45","#36F58A","#0A98F2","#0A2CF0",
                     "#FFC414","#FF5001","#E8190C","#FF1398")

legend(x='right', legend = c("no trend",
                             "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
                             "interrupted decrease","decrease after no trend","no trend after increase","negative reversal"),
       fill = bp_type_pal_vis, bty='n', cex=1, xpd=NA,inset=-0.65)

dev.off()

# Flood breakpoint typology ---------------------------------------------

png(filename="FBPtypology.png", family="Calibri",width=2400, height=1600, res=300)
par(mfrow = c(1,1), mar=c(5,4,4,16)+0.1)

# webcolors -- use the other pallete (see illustrator file) for print colors
bp_type_pal <- c("#575756",
                 "#81FA45","#36F58A","#E8190C",
                 "#FFC414","#FF5001","#0A98F2",
                 "#0A2CF0","#FF1398")

plot(output$FBPType,legend=FALSE,col=bp_type_pal,zlim=c(0,8),
     main="Typology of flood breakpoints",xlab="Longitude",ylab="Latitude")

bp_type_pal_vis <- c("#575756",
                     "#81FA45","#36F58A","#0A98F2","#0A2CF0",
                     "#FFC414","#FF5001","#E8190C","#FF1398")

legend(x='right', legend = c("no trend",
                             "interrupted increase","increase after no trend","no trend after decrease","positive reversal",
                             "interrupted decrease","decrease after no trend","no trend after increase","negative reversal"),
       fill = bp_type_pal_vis, bty='n', cex=1, xpd=NA,inset=-0.65)

dev.off()
