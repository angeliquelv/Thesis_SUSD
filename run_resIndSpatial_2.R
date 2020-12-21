# Install packages and load libraries -------------------------------------

# Most packages are already pre-installed on Lisa
# bfastSpatial should be installed manually to my personal directory /home/avermeer/R/library,
# this is I did interactively on the home node by running devtools::install_github('loicdtx/bfastSpatial',lib="/home/avermeer/R/library" in R
# option 3 (do not update packages) was chosen
# the dependencies gdalUtils, bfast and rgdal are also installed to the personal directory
# these 3 have to be loaded manually with the directory defined (these will not automatically be found)


# load libraries
library(raster) # package for raster manipulation
library(parallel)
library(strucchange)
library(MASS)
library(zoo)
library(rgdal,lib.loc="/home/avermeer/R/library")
library(gdalUtils,lib.loc="/home/avermeer/R/library")
library(bfast,lib.loc="/home/avermeer/R/library")
library(bfastSpatial,lib.loc="/home/avermeer/R/library")

# Set input and output directory ------------------------------------------

workdirIN <- paste0(getwd(),"/input_dir_2")
workdirOUT <- getwd()

# Updated resIndSpatial function ------------------------------------------

resIndSpatial <- function(x, dates, type='irregular', sc=1, order=3,
                          formula = response ~ (trend + harmon), h=0.15,
                          plevel=0.05, dr, s=3, fr, NV=NA, mc.cores=1) {

  # maximum number of breaks following h
  max_breaks <- floor((1/h)-1)
  namesbd <- c(paste0("bd",1:max_breaks))
  namesloCIbd <- c(paste0("loCIbd",1:max_breaks))
  namesupCIbd <- c(paste0("upCIbd",1:max_breaks))
  
  fun <- function(x){

    ##Set output vector to NA
    resind <- NA

    ##Set own NoData value to differentiate between "no data pixels" and
    #"no breakpoint pixels"
    NV=NV

    #Apply bfastts() and multiply by data by scaling factor if required
    bfts <- bfastts(x*sc,dates,type=type)

    ##If not all oberservation NA (e.g. pixels outside area of interest) run procedure
    if(!all(is.na(bfts))){
      #Create data frame from bfastts
      bpp <- bfastpp(bfts, order=order, stl=("none"), na.action=na.omit)

      ##Calculate initial NDVI: First s years after start of observations
      Ini <- subset(bpp, time >= bpp$time[1] & time <= bpp$time[1]+s)
      MIni <- mean(Ini$response)

      ##Calculate End NDVI: Last s years after before end of observations
      End <- subset(bpp, time >= bpp$time[length(bpp$time)]-s & time <= bpp$time[length(bpp$time)])
      MEnd <- mean(End$response)

      ##MOSUM test for stuctural stability
      Mos <- efp(formula, data=bpp, type="OLS-MOSUM", h=h, rescale=FALSE)
      ##Calculate test statistics
      sct <- sctest(Mos)

      ##Fit breakpoints and segmented model if empirical fluctuation test was significant
      if (sct$p.value <= plevel){
        ##Fit the breakpoints
        #This is written in a tryCatch because sometimes an error in the RSS table occurs
        #It happens more in pixels with a smaller h and with many NA values
        #Probably due to singularities in the regression
        #However, it's unpredictable, so this tryCatch filters out those pixels
        #All procedures are terminated for these pixels and all indicators set to NA
        bpoints <- tryCatch({
          breakpoints(formula=formula, data=bpp, h=h)
        },
        error=function(cond) {
          message("Original error message:")
          message("Error in my.RSS.table[as.character(i), 3:4] <- c(pot.index[opt], break.RSS[opt]) :")
          message(cond[1])
          return(NA)
        }
        )
        # if the error did not occur and breakpoinits have been calculated, continue the procedure
        if(length(bpoints)==1){
          # if bpoints length=1 (NA) set bfts to NA too
          bfts <- NA
        } else {
          p <- bpoints$breakpoints
          if((length(p)==1 && is.na(p))) {
            #If no breakpoint p is NA. set breakpoint number to zero. Give warning.
            warning('No breakpoint found')
            bpnumb <- 0
            #Set breakpoint parameters to NA (needed for further calculations)
            bd <- NA
            cid <- NA
            # set breakdates to NV
            bd_vector <- vector()
            for (i in 1:max_breaks){
              bd_vector[i] <- assign(namesbd[i], NV)
            }
            # set lower CI around breakdates to NV
            loCIbd_vector <- vector()
            for (i in 1:max_breaks){
              loCIbd_vector[i] <- assign(namesloCIbd[i], NV)
            }
            # set upper CI around breakdates to NV
            upCIbd_vector <- vector()
            for (i in 1:max_breaks){
              upCIbd_vector[i] <- assign(namesupCIbd[i], NV)
            }
            
            #Fit unsegmented model
            m <- rlm (formula=formula, data=bpp, maxit=100)
            bpp$segment <- 1 ##even if you have only one segment you need to specify this,
            #for trend calculation later on
          } else {
            #If there is a breakpoint extract breakpoint parameters
            #breakpoint number
            bpnumb <- length(p)
            #breakdates
            bd <- bpp$time[p]
            #confidence invervals
            ci <- confint(bpoints, breakpoints=TRUE)
            bda <- ci[[1]]

            #Calculate the confidence interval around the breakpoint timing
            tryCatch({
              CI_lower <- bpp$time[bda[,1]]
              CI_upper <- bpp$time[bda[,3]]
            },
            error=function(cond) {
              message("CI breakdate outside time range, Original error message:")
              message(cond[1])
              return(NA)
            }
            )
            
            ## If breakpoint extract breakdates
            bd_vector <- vector()
            for (i in 1:max_breaks){
              bd_vector[i] <- assign(namesbd[i], bd[i])
            }
            ## If breakpoint extract lower CI limit
            loCIbd_vector <- vector()
            for (i in 1:max_breaks){
              loCIbd_vector[i] <- assign(namesloCIbd[i], CI_lower[i])
            }   
            ## If breakpoint extract upper CI limit
            upCIbd_vector <- vector()
            for (i in 1:max_breaks){
              upCIbd_vector[i] <- assign(namesupCIbd[i], CI_upper[i])
            }
            
            ##Fit segmented model
            #Create column "segment" that indices segment in the dataframe "bpp"
            bpp$segment <- breakfactor(bpoints)
            #RLM fit; I increased the default maxiterations (20 -> 100)
            m <- rlm (response ~ segment/(trend+harmon), data=bpp, maxit=100)
          }
        }
      } else {
        ##If MOSUM not significant set breakpoint number and DBP to zero and other
        #output variables to NV
        bpnumb <- 0
        # set breakdates to NV
        bd_vector <- vector()
        for (i in 1:max_breaks){
          bd_vector[i] <- assign(namesbd[i], NV)
        }
        # set lower CI around breakdates to NV
        loCIbd_vector <- vector()
        for (i in 1:max_breaks){
          loCIbd_vector[i] <- assign(namesloCIbd[i], NV)
        }
        # set upper CI around breakdates to NV
        upCIbd_vector <- vector()
        for (i in 1:max_breaks){
          upCIbd_vector[i] <- assign(namesupCIbd[i], NV)
        }
        #Fit unsegmented model
        m <- rlm (formula=formula, data=bpp, maxit=100)
        bpp$segment <- 1 ##even if you have only one segment you need to specify this,
        #for trend calculation later on
        warning('No significant deviation from structural stability (MOSUM
              test). No breakpoints fitted.')
      }
      if(!all(is.na(bfts))){

        #Predict values and add column "prediction" in in the dataframe "bpp"
        bpp$prediction <- predict(m,newdata=bpp)

        #Add trend prediction for each segment. Corrected hight of trend line for
        #irregular data: sets harmonic term based on mean DOY of observations/segment
        #(instead of trend$harmon[] <- 0, which assumes regular data);
        #idea based on email exchange with Achim Zeileis
        for(i in 1:length(unique(bpp$segment))) {
          seg <- unique(bpp$segment)[i]
          # trend <- subset(bpp, segment == seg)
          trend <- bpp[bpp$segment==seg, ]
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
        
        if (length(trends)==1) {
          Trend_nobp <- trends[1]
          Trend_nobp_lowerCI <- ci_trend[1]
          Trend_nobp_upperCI <- ci_trend[2]
        } else {
          Trend_nobp <- NV
          Trend_nobp_lowerCI <- NV
          Trend_nobp_upperCI <- NV
        }

        ##Extract Amplitudes for each segment. Depending on parameter "order" one
        # gets multiple sine and cosine terms
        cosines <- coef[which(grepl('harmoncos', names(coef)))]
        sines <- coef[which(grepl('harmonsin', names(coef)))]
        amps <- sqrt(cosines^2 + sines^2)
        names(amps) <- rep(1:length(unique(bpp$segment)), order)

        if(bpnumb!=0) {
          # Determine the type of breakpoint that occured,
          # count of breakpoints following type 1,2,3,4,5 or 6:
          # Initialise all as zero
          bpnumb_type0 <- 0
          bpnumb_type1 <- 0
          bpnumb_type2 <- 0
          bpnumb_type3 <- 0
          bpnumb_type4 <- 0
          bpnumb_type5 <- 0
          bpnumb_type6 <- 0
          bpnumb_type7 <- 0
          bpnumb_type8 <- 0

          for (i in 1:bpnumb) {
            trendrecov <- trends[i+1]
            pretrend <- trends[i]
            ci_trendrecov <- ci_trend[i+1,]
            ci_pretrend <- ci_trend[i,]

            # no significant trend before and after breakpoint
            if (prod(ci_pretrend)<0 && prod(ci_trendrecov)<0) {
              bpnumb_type0 <- bpnumb_type0 + 1
            }
            # interrupted increase, significant increase before and after breakpoint
            if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend>0 && trendrecov>0) {
              bpnumb_type1 <- bpnumb_type1 + 1
            }
            # increase after period of no significant trend
            if (prod(ci_pretrend)<0 && prod(ci_trendrecov)>0 && trendrecov>0) {
              bpnumb_type2 <- bpnumb_type2 + 1
            }
            # period of no significant trend after increase
            if (prod(ci_pretrend)>0 && prod(ci_trendrecov)<0 && pretrend>0) {
              bpnumb_type3 <- bpnumb_type3 + 1
            }
            # interrupted decrease, significant decrease before and after breakpoint
            if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend<0 && trendrecov<0) {
              bpnumb_type4 <- bpnumb_type4 +1
            }
            # decrease after period of no significant trend
            if (prod(ci_pretrend)<0 && prod(ci_trendrecov)>0 && trendrecov<0) {
              bpnumb_type5 <- bpnumb_type5 +1
            }
            # period of no significant trend after decrease
            if (prod(ci_pretrend)>0 && prod(ci_trendrecov)<0 && pretrend<0) {
              bpnumb_type6 <- bpnumb_type6 +1
            }
            # positive reversal
            if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend<0 && trendrecov>0) {
              bpnumb_type7 <- bpnumb_type7 +1
            }
            # negative reversal
            if (prod(ci_pretrend)>0 && prod(ci_trendrecov)>0 && pretrend>0 && trendrecov<0) {
              bpnumb_type8 <- bpnumb_type8 + 1
            }
          }
        }
        # Calculate (drought and flood) resilience indicators ---------------------

        ##If no (significant) breakpoint found set breakpoint position, breakpoint timing,
        #and recovery trend to "NV"; set DBP to 0
        if(bpnumb==0){
          bpnumb_type0 <- 0
          bpnumb_type1 <- 0
          bpnumb_type2 <- 0
          bpnumb_type3 <- 0
          bpnumb_type4 <- 0
          bpnumb_type5 <- 0
          bpnumb_type6 <- 0
          bpnumb_type7 <- 0
          bpnumb_type8 <- 0
          DBP <- 0
          DBP_bpt <- NV
          DBP_CI_lower <- NV
          DBP_CI_upper <- NV
          DBP_trendrecov <- NV
          DBP_pretrend <- NV
          DBP_preNDVI <- NV
          DBP_MagA <- NV
          DBP_MagR <- NV
          DBP_AmpDiff <- NV
          DBP_Type <- NV
          FBP <- 0
          FBP_bpt <- NV
          FBP_CI_lower <- NV
          FBP_CI_upper <- NV
          FBP_trendrecov <- NV
          FBP_pretrend <- NV
          FBP_preNDVI <- NV
          FBP_MagA <- NV
          FBP_MagR <- NV
          FBP_AmpDiff <- NV
          FBP_Type <- NV
          #Set breakpoint position to NA (needed for further calculations)
          DBP_bpd <- NA
          FBP_bpd <- NA
        } else {
          ##If breakpoint, check if one breakpoint ocurred around drought period ("drought breakpoint")
          DBP_bpd <- match(bd[bd >= dr[1] & bd <= dr[2]], bd) #bpd gives you the position of bp.

          ##If drought breakpoint occured set DBP to 1 and calculate breakpoint timing,
          #trend of recovery, trend in segment before BP, preNDVI, Magnitude of change
          #(MagA & MagR) based on mean NDVI of s years before and after BP.
          #If there is more than 1 BP in the specified time interval, the first one is chosen!
          if (length(DBP_bpd) > 1) {
            DBP_bpd <- DBP_bpd[1]
          }

          if (length(DBP_bpd) > 0) {
            DBP <- 1

            #Calculate BP timing
            DBP_bpt <- bd[DBP_bpd]

            #Calculate the confidence interval around the breakpoint timing
            #This is written in a trycatch to avoid an error,
            #in case the confidence interval is outside the period of bpp$time
            tryCatch({
              DBP_ci_dates <- bpp$time[c(ci[["confint"]][DBP_bpd,1],ci[["confint"]][DBP_bpd,3])]
              DBP_CI_lower <- DBP_ci_dates[1]
              DBP_CI_upper <- DBP_ci_dates[2]
            },
            error=function(cond) {
              message("Original error message:")
              message(cond[1])
              return(NA)
            }
            )

            #Calculate trend slopes
            seg2 <- (DBP_bpd+1)
            DBP_trendrecov <- trends[seg2]
            seg1 <- (DBP_bpd)
            DBP_pretrend <- trends[seg1]
            
            #Calculate Confidence Intervals trends slopes
            DBP_ci_trendrecov <- ci_trend[seg2,]
            DBP_ci_pretrend <- ci_trend[seg1,]

            #Calculate preNDVI, Magnitude of change (MagA & MagR) based on mean NDVI
            #of s years before and after BP. Use subset of bpp$response that falls
            #within date vector
            ti <- c(DBP_bpt-s, DBP_bpt+s) # +/- s years around breakpoint
            S1 <- subset(bpp, time >= ti[1] & time <= DBP_bpt)
            S2 <- subset(bpp, time > DBP_bpt & time <= ti[2])

            DBP_preNDVI <- mean(S1$response)
            M2 <- mean(S2$response)

            DBP_MagA <- (M2-DBP_preNDVI) #absolute change magnitude
            DBP_MagR <- DBP_MagA/DBP_preNDVI #relative change magnitude

            #Calculate Difference in mean amplitudes between segments
            mean_ampsseg1 <- mean(amps[which(grepl(seg1,names(amps)))]) #mean amplitude in segment before bp
            mean_ampsseg2 <- mean(amps[which(grepl(seg2,names(amps)))]) #mean amplitude in segment after bp

            DBP_AmpDiff <- (mean_ampsseg2-mean_ampsseg1)/mean_ampsseg1 #relative difference between mean amplitudes

            # Determine the breakpoint typology
            # no significant trend before and after breakpoint
            if (prod(DBP_ci_pretrend)<0 && prod(DBP_ci_trendrecov)<0) {
              DBP_Type <- 0
            }
            # interrupted increase, significant increase before and after breakpoint
            if (prod(DBP_ci_pretrend)>0 && prod(DBP_ci_trendrecov)>0 && DBP_pretrend>0 && DBP_trendrecov>0) {
              DBP_Type <- 1
            }
            # increase after period of no significant trend
            if (prod(DBP_ci_pretrend)<0 && prod(DBP_ci_trendrecov)>0 && DBP_trendrecov>0) {
              DBP_Type <- 2
            }
            # period of no significant trend after increase
            if (prod(DBP_ci_pretrend)>0 && prod(DBP_ci_trendrecov)<0  && pretrend>0) {
              DBP_Type <- 3
            }
            # interrupted decrease, significant decrease before and after breakpoint
            if (prod(DBP_ci_pretrend)>0 && prod(DBP_ci_trendrecov)>0 && DBP_pretrend<0 && DBP_trendrecov<0) {
              DBP_Type <- 4
            }
            # decrease after period of no significant trend
            if (prod(DBP_ci_pretrend)<0 && prod(DBP_ci_trendrecov)>0 && DBP_trendrecov<0) {
              DBP_Type <- 5
            }
            # period of no significant trend after decrease
            if (prod(DBP_ci_pretrend)>0 && prod(DBP_ci_trendrecov)<0 && pretrend<0) {
              DBP_Type <- 6
            }
            # positive reversal
            if (prod(DBP_ci_pretrend)>0 && prod(DBP_ci_trendrecov)>0 && DBP_pretrend<0 && DBP_trendrecov>0) {
              DBP_Type <- 7
            }
            # negative reversal
            if (prod(DBP_ci_pretrend)>0 && prod(DBP_ci_trendrecov)>0 && DBP_pretrend>0 && DBP_trendrecov<0) {
              DBP_Type <- 8
            }
          } else {
            ##If no breakpoint around drought occured set patrameters to NA
            warning('No drought breakpoint found')
            DBP <- 0
            DBP_bpt <- NV
            DBP_CI_lower <- NV
            DBP_CI_upper <- NV
            DBP_trendrecov <- NV
            DBP_pretrend <- NV
            DBP_preNDVI <- NV
            DBP_MagA <- NV
            DBP_MagR <- NV
            DBP_AmpDiff <- NV
            DBP_Type <- NV
          }
          ##If breakpoint, check if one breakpoint ocurred around flood period ("flood breakpoint")
          FBP_bpd <- match(bd[bd >= fr[1] & bd <= fr[2]], bd) #bpd gives you the position of bp.

          ##If flood breakpoint occured set DBP to 1 and calculate breakpoint timing,
          #trend of recovery, trend in segment before BP, preNDVI, Magnitude of change
          #(MagA & MagR) based on mean NDVI of s years before and after BP.
          #If there is more than 1 BP in the specified time interval, the first one is chosen!
          if (length(FBP_bpd) > 1) {
            FBP_bpd <- FBP_bpd[1]
          }

          if (length(FBP_bpd) > 0) {
            FBP <- 1

            #Calculate BP timing
            FBP_bpt <- bd[FBP_bpd]

            #Calculate the confidence interval around the breakpoint timing
            tryCatch({
              FBP_ci_dates <- bpp$time[c(ci[["confint"]][FBP_bpd,1],ci[["confint"]][FBP_bpd,3])]
              FBP_CI_lower <- FBP_ci_dates[1]
              FBP_CI_upper <- FBP_ci_dates[2]
            },
            error=function(cond) {
              message("Original error message:")
              message(cond[1])
              return(NA)
            }
            )

            #Calculate trend slopes
            seg2 <- (FBP_bpd+1)
            FBP_trendrecov <- trends[seg2]
            seg1 <- (FBP_bpd)
            FBP_pretrend <- trends[seg1]
            
            #Calculate Confidence Intervals trends slopes
            FBP_ci_trendrecov <- ci_trend[seg2,]
            FBP_ci_pretrend <- ci_trend[seg1,]

            #Calculate preNDVI, Magnitude of change (MagA & MagR) based on mean NDVI
            #of s years before and after BP. Use subset of bpp$response that falls
            #within date vector
            ti <- c(FBP_bpt-s, FBP_bpt+s) # +/- s years around breakpoint
            S1 <- subset(bpp, time >= ti[1] & time <= FBP_bpt)
            S2 <- subset(bpp, time > FBP_bpt & time <= ti[2])

            FBP_preNDVI <- mean(S1$response)
            M2 <- mean(S2$response)

            FBP_MagA <- (M2-FBP_preNDVI) #absolute change magnitude
            FBP_MagR <- FBP_MagA/FBP_preNDVI #relative change magnitude

            #Calculate Difference in mean amplitudes between segments
            mean_ampsseg1 <- mean(amps[which(grepl(seg1,names(amps)))]) #mean amplitude in segment before bp
            mean_ampsseg2 <- mean(amps[which(grepl(seg2,names(amps)))]) #mean amplitude in segment after bp

            FBP_AmpDiff <- (mean_ampsseg2-mean_ampsseg1)/mean_ampsseg1 #relative difference between mean amplitudes

            # Determine the breakpoint typology
            # no significant trend before and after breakpoint
            if (prod(FBP_ci_pretrend)<0 && prod(FBP_ci_trendrecov)<0) {
              FBP_Type <- 0
            }
            # interrupted increase, significant increase before and after breakpoint
            if (prod(FBP_ci_pretrend)>0 && prod(FBP_ci_trendrecov)>0 && FBP_pretrend>0 && FBP_trendrecov>0) {
              FBP_Type <- 1
            }
            # increase after period of no significant trend
            if (prod(FBP_ci_pretrend)<0 && prod(FBP_ci_trendrecov)>0 && FBP_trendrecov>0) {
              FBP_Type <- 2
            }
            # period of no significant trend after increase
            if (prod(FBP_ci_pretrend)>0 && prod(FBP_ci_trendrecov)<0  && pretrend>0) {
              FBP_Type <- 3
            }
            # interrupted decrease, significant decrease before and after breakpoint
            if (prod(FBP_ci_pretrend)>0 && prod(FBP_ci_trendrecov)>0 && FBP_pretrend<0 && FBP_trendrecov<0) {
              FBP_Type <- 4
            }
            # decrease after period of no significant trend
            if (prod(FBP_ci_pretrend)<0 && prod(FBP_ci_trendrecov)>0 && FBP_trendrecov<0) {
              FBP_Type <- 5
            }
            # period of no significant trend after decrease
            if (prod(FBP_ci_pretrend)>0 && prod(FBP_ci_trendrecov)<0 && pretrend<0) {
              FBP_Type <- 6
            }
            # positive reversal
            if (prod(FBP_ci_pretrend)>0 && prod(FBP_ci_trendrecov)>0 && FBP_pretrend<0 && FBP_trendrecov>0) {
              FBP_Type <- 7
            }
            # negative reversal
            if (prod(FBP_ci_pretrend)>0 && prod(FBP_ci_trendrecov)>0 && FBP_pretrend>0 && FBP_trendrecov<0) {
              FBP_Type <- 8
            }
          } else {
            ##If no breakpoint around flood occured set patrameters to NA
            warning('No flood breakpoint found')
            FBP <- 0
            FBP_bpt <- NV
            FBP_CI_lower <- NV
            FBP_CI_upper <- NV
            FBP_trendrecov <- NV
            FBP_pretrend <- NV
            FBP_preNDVI <- NV
            FBP_MagA <- NV
            FBP_MagR <- NV
            FBP_AmpDiff <- NV
            FBP_Type <- NV
          }
        }
      } else {
        # If bpoints is NA due to error in bpoint calculation:
        warning('Error occurred during computation of breakpoints. All ouptut variables are set to NA')
        # set breakdates to NA
        bd_vector <- vector()
        for (i in 1:max_breaks){
          bd_vector[i] <- assign(namesbd[i], NA)
        }
        # set lower CI around breakdates to NA
        loCIbd_vector <- vector()
        for (i in 1:max_breaks){
          loCIbd_vector[i] <- assign(namesloCIbd[i], NA)
        }
        # set upper CI around breakdates to NA
        upCIbd_vector <- vector()
        for (i in 1:max_breaks){
          upCIbd_vector[i] <- assign(namesupCIbd[i], NA)
        }
        bpnumb <- NA
        bpnumb_type0 <- NA
        bpnumb_type1 <- NA
        bpnumb_type2 <- NA
        bpnumb_type3 <- NA
        bpnumb_type4 <- NA
        bpnumb_type5 <- NA
        bpnumb_type6 <- NA
        bpnumb_type7 <- NA
        bpnumb_type8 <- NA
        MIni <- NA
        MEnd <- NA
        Int <- NA
        Trend_nobp <- NA
        Trend_nobp_lowerCI <- NA
        Trend_nobp_upperCI <- NA
        DBP <- NA
        DBP_bpt <- NA
        DBP_CI_lower <- NA
        DBP_CI_upper <- NA
        DBP_trendrecov <- NA
        DBP_pretrend <- NA
        DBP_preNDVI <- NA
        DBP_MagA <- NA
        DBP_MagR <- NA
        DBP_AmpDiff <- NA
        DBP_Type <- NA
        FBP <- NA
        FBP_bpt <- NA
        FBP_CI_lower <- NA
        FBP_CI_upper <- NA
        FBP_trendrecov <- NA
        FBP_pretrend <- NA
        FBP_preNDVI <- NA
        FBP_MagA <- NA
        FBP_MagR <- NA
        FBP_AmpDiff <- NA
        FBP_Type <- NA
      }
    } else {
      #If all values NA (pixels without NDVI values!) set output variables to NA and give warning
      warning('No observations in time series')
      # set breakdates to NA
      bd_vector <- vector()
      for (i in 1:max_breaks){
        bd_vector[i] <- assign(namesbd[i], NA)
      }
      # set lower CI around breakdates to NA
      loCIbd_vector <- vector()
      for (i in 1:max_breaks){
        loCIbd_vector[i] <- assign(namesloCIbd[i], NA)
      }
      # set upper CI around breakdates to NA
      upCIbd_vector <- vector()
      for (i in 1:max_breaks){
        upCIbd_vector[i] <- assign(namesupCIbd[i], NA)
      }
      bpnumb <- NA
      bpnumb_type0 <- NA
      bpnumb_type1 <- NA
      bpnumb_type2 <- NA
      bpnumb_type3 <- NA
      bpnumb_type4 <- NA
      bpnumb_type5 <- NA
      bpnumb_type6 <- NA
      bpnumb_type7 <- NA
      bpnumb_type8 <- NA
      MIni <- NA
      MEnd <- NA
      Int <- NA
      Trend_nobp <- NA
      Trend_nobp_lowerCI <- NA
      Trend_nobp_upperCI <- NA
      DBP <- NA
      DBP_bpt <- NA
      DBP_CI_lower <- NA
      DBP_CI_upper <- NA
      DBP_trendrecov <- NA
      DBP_pretrend <- NA
      DBP_preNDVI <- NA
      DBP_MagA <- NA
      DBP_MagR <- NA
      DBP_AmpDiff <- NA
      DBP_Type <- NA
      FBP <- NA
      FBP_bpt <- NA
      FBP_CI_lower <- NA
      FBP_CI_upper <- NA
      FBP_trendrecov <- NA
      FBP_pretrend <- NA
      FBP_preNDVI <- NA
      FBP_MagA <- NA
      FBP_MagR <- NA
      FBP_AmpDiff <- NA
      FBP_Type <- NA
    }

    # Save and return output --------------------------------------------------

    #Save all indicators::
    resind <- cbind(bpnumb, bpnumb_type0, bpnumb_type1, bpnumb_type2, bpnumb_type3, 
                    bpnumb_type4, bpnumb_type5, bpnumb_type6, bpnumb_type7, bpnumb_type8, 
                    MIni, MEnd, as.numeric(Int), as.numeric(Trend_nobp), as.numeric(Trend_nobp_lowerCI),
                    as.numeric(Trend_nobp_upperCI),t(bd_vector),t(loCIbd_vector),t(upCIbd_vector),
                    DBP, DBP_bpt, DBP_CI_lower, DBP_CI_upper,
                    as.numeric(DBP_trendrecov), as.numeric(DBP_pretrend), DBP_preNDVI, DBP_MagA, DBP_MagR,
                    DBP_AmpDiff, DBP_Type, FBP, FBP_bpt, FBP_CI_lower,
                    FBP_CI_upper, as.numeric(FBP_trendrecov), as.numeric(FBP_pretrend), FBP_preNDVI, FBP_MagA,
                    FBP_MagR, FBP_AmpDiff, FBP_Type)

    colnames(resind) <- c('BPNumb', 'BPNumbType0', 'BPNumbType1','BPNumbType2','BPNumbType3',
                          'BPNumbType4', 'BPNumbType5','BPNumbType6','BPNumbType7','BPNumbType8',
                          'Initial NDVI', 'End NDVI', 'Intercept', 'Trend_noBP','Trend_noBP_lowerCI',
                          'Trend_noBP_upperCI',namesbd,namesloCIbd,namesupCIbd,
                          'DBP','DBPTime','DBPlowerCI','DBPupperCI',
                          'RecTrendDBP', 'PreTrendDBP', 'DBPPreNDVI', 'DBPMagObsA', 'DBPMagObsR',
                          'DBPAmpDiffR', 'DBPType', 'FBP','FBPTime','FBPlowerCI',
                          'FBPupperCI', 'RecTrendFBP', 'PreTrendFBP', 'FBPPreNDVI', 'FBPMagObsA',
                          'FBPMagObsR', 'FBPAmpDiffR','FBPType')

    return(resind)
  }

  out <- bfastSpatial::mc.calc(x=x, fun=fun, mc.cores=mc.cores)
  names(out) <- c('BPNumb', 'BPNumbType0', 'BPNumbType1','BPNumbType2','BPNumbType3',
                  'BPNumbType4', 'BPNumbType5','BPNumbType6','BPNumbType7','BPNumbType8',
                  'Initial NDVI', 'End NDVI', 'Intercept', 'Trend_noBP','Trend_noBP_lowerCI',
                  'Trend_noBP_upperCI',namesbd,namesloCIbd,namesupCIbd, 
                  'DBP','DBPTime','DBPlowerCI','DBPupperCI',
                  'RecTrendDBP', 'PreTrendDBP', 'DBPPreNDVI', 'DBPMagObsA', 'DBPMagObsR',
                  'DBPAmpDiffR', 'DBPType', 'FBP','FBPTime','FBPlowerCI',
                  'FBPupperCI', 'RecTrendFBP', 'PreTrendFBP', 'FBPPreNDVI', 'FBPMagObsA',
                  'FBPMagObsR', 'FBPAmpDiffR','FBPType')
  
  return(out)
}

# Load NDVI stack ---------------------------------------------------------

setwd(workdirIN)

NDVI <- brick(infile_rscript)

# create dataframe with dates
dates <- as.Date(substr(names(NDVI),5,12),format="%Y%m%d")

# Apply resIndSpatial and save results ------------------------------------

output <- resIndSpatial(NDVI, dates,
                        type='irregular', sc=1, order=2,
                        formula = response ~ (trend + harmon),
                        h=1/10, plevel=0.05,
                        dr=c(1998.748,2002.745), s=3, fr=c(2014.748,2015.745),
                        mc.cores=24)

setwd(workdirOUT)

writeRaster(output,outfile_rscript)
