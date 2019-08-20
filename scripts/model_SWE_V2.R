##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script calculates uses the hydromad model in the R package http://hydromad.catchment.org/
#to estimate snow water equilavent (SWE) from daily mean temperature and daily total precipitation data
#and then to use temperature thresholds that best predicted SWE to partition total precipitation into
#snowfall and rainfall

####################################################################################
#Initial set up                                                                   
####################################################################################

#call libraries

library(plyr)
library(hydromad)
library(matrixStats)
library(zoo)

######################################################################################################
#run  hydromad model for each station and year in Daily_Complete to simulate SWE for all stations
######################################################################################################

#call file names in COMBSWE subdirectory
files.COMP = list.files(path = "Daily_Complete/",full.names = T)

#create blank data.frame for containing SWE and SNOW model validation results
swemods = data.frame(matrix(vector(), nrow = length(files.COMP), ncol = 12))

#add column names to swemods
names(swemods) = c("station", "wyear_lo_swe", "wyear_hi_swe", "nyears_swe", "ndat_swe", "mod_ID_swe", "rsq_swe",
                        "wyear_lo_snow", "wyear_hi_snow", "nyears_snow", "ndat_snow", "rsq_snow")

#start loop through files.COMP

for(j in seq(length(files.COMP))){
  alldat = read.table(files.COMP[j], head = T, sep = ",")

  station = as.character(unique(alldat$station))
  swemods$station[j] = station
  
  
  #create new column for water year, which begins July 1 and ends June 30 of the followig year for SWE modeling
  alldat$wyear = ifelse(alldat$doy < 182, alldat$year - 1, alldat$year)
  
  #determine start and end of SNWD record
  SNWD = alldat[complete.cases(alldat$SNWD), ]
  
  #add earliest and latest years and total number of years for which there are WESD data to swemods
  swemods$wyear_lo_swe[j] = min(SNWD$wyear)
  swemods$wyear_hi_swe[j] = max(SNWD$wyear)
  swemods$nyears_swe[j] = max(SNWD$wyear) - min(SNWD$wyear)
  
  #############################################
  #gap fill temperature and precipitation data#
  #############################################
  
  #gap fill missing temperature and precipitation data because the model will fail without complete data
  #temperature data can be interpolated
  alldat$TMAXgap = ifelse(is.na(alldat$TMAX) == T, na.approx(alldat$TMAX, rule = 2, na.rm = F), alldat$TMAX)
  alldat$TMINgap = ifelse(is.na(alldat$TMIN) == T, na.approx(alldat$TMIN, rule = 2, na.rm = F), alldat$TMIN)
  
  #precipitation data assumed to be zero
  alldat$PRCPgap = ifelse(is.na(alldat$PRCP) == T, 0, alldat$PRCP)
  
  #calcuate daily average temperature
  alldat$TMEAN = (alldat$TMAXgap + alldat$TMINgap) / 2
  
  #calculate the average of TMEAN and TMIN
  alldat$TMINE = (alldat$TMEAN + alldat$TMINgap) / 2
  
  #convert into units that the hydromad model expects (degrees C for temperature, mm for precipitation)
  alldat$TMEANc = alldat$TMEAN / 10
  alldat$TMINEc = alldat$TMINE / 10
  alldat$PRCPm = alldat$PRCPgap / 10
  
  ######################
  #select model subsets#
  ######################
  
  #omit duplicate years to determine when to start / end model subsets. This allows SWE modeling for discrete winters
  start.rows = !duplicated(alldat$wyear)
  end.rows = !duplicated(alldat$wyear, fromLast = TRUE)
  
  #extract original row indexes associated with the start and end of a wyear
  sr.ind <- seq_along(alldat$wyear)[start.rows]
  er.ind <- seq_along(alldat$wyear)[end.rows]
  
  #check that the number of subsets is the same for the start and end
  print(length(sr.ind))
  print(length(er.ind))
  
  #########################
  #select model parameters#
  #########################

  #create vectors with all possible values for temperature thresholds for rainfall (Tmax) and snowmelt (Tmin)
  #and all possible values for melt coefficient (kd)
  
  Tmax = seq(-2, 2, by = 0.5)
  Tmin = seq(-2, 2, by = 0.5)
  kd = seq(1, 6, by = 0.5)
  
  #compute all possible combinations of Tmax and Tmin
  tpar = expand.grid(Tmax, Tmin)
  names(tpar) = c("Tmax", "Tmin")
  
  #select combinations where Tmax is > Tmin (model will not run if Tmax < Tmin
  tpar = tpar[tpar$Tmax > tpar$Tmin, ]
  
  #extract Tmax and Tmin for merging with kd
  Tmax = tpar$Tmax
  Tmin = tpar$Tmin
  
  #compute all possible combinations of Tmax, Tmin, and kd
  pars = data.frame(expand.grid(Tmax, Tmin, kd))
  names(pars) = c("Tmax", "Tmin", "kd")
  
  #select combinations where Tmax is > Tmin
  pars = pars[pars$Tmax > pars$Tmin, ]
  
  #make id column for removing duplicates
  pars$id = paste(pars$Tmax, pars$Tmin, pars$kd, sep = " ")
  pars = subset(pars, !duplicated(pars$id))

  #subset to only include Tmax, Tmin, and kd
  pars = pars[ , c("Tmax", "Tmin", "kd")]
  
  #transpose so that columns become rows (makes it easier to fill WEI)
  pars = t(pars)
  
  #create container to hold results of model
  WEI = data.frame(matrix(nrow = nrow(alldat), ncol = ncol(pars)))
  
  ##########
  #Run model
  ##########
  
  #begin loop to select precipitation and temperature data
  for (h in seq(1,length(sr.ind))) {
  
    #select precipitation subsets for each wyear
    P = ts(alldat$PRCPm[sr.ind[h]:er.ind[h]])
    #select temperature subsets for each wyear
    E = ts(alldat$TMEANc[sr.ind[h]:er.ind[h]])
    #combine P and E into single file
    met.1 = ts.union(P,E)
    
    #begin loop to select parameter sets
    for (i in 1:ncol(pars)) {
      
    #loop across pars columns to choose parameter suite of Tmax, Tmin, and kd
    Tmax = pars[1, i]
    Tmin = pars[2, i]
    kd = pars[3, i]
    
    #run hydromad model
    metmod = snow.sim(met.1, Tmax = Tmax, Tmin = Tmin, cr = 1, cs = 0.88,
                      kd = kd, kf = 0, rcap = 0,
                      d = 0, f = 0, e = 0, LSWE_0 = 0, ISWE_0 = 0, return_state = T)
    
    #write output of model runs for each parameter suite
    WEI[[i]][sr.ind[h]:er.ind[h]] =  metmod[ , 'SWE']
    
    }
  }
  
  #merge WEI with alldat
  alldat.1 = cbind(alldat, WEI)
  
  #change column names of predicted values to something more intuitive and consistent
  pred = data.frame("predno" = c(1:396))
  pred$prednam = ifelse(pred$predno < 10, paste("pred", "000", pred$predno, sep = ""),
                        ifelse(pred$predno < 100, paste("pred", "00", pred$predno, sep = ""),
                               ifelse(pred$predno < 1000, paste("pred", "0", pred$predno, sep = ""),
                                      paste("pred", pred$predno, sep = ""))))
  
  names(alldat.1)[28:423] = pred$prednam
  
  #####################################
  #Compare simulated to observed values
  #####################################
  
  #omit all NAs and zero values (since zeros can occur due to failure to record data, not because snow was not present)
  snoz = alldat.1[complete.cases(alldat.1[ , "SNWD"]) & alldat.1$SNWD != 0, ]
  
  #determine number of records in dataset
  SNWDval = snoz$SNWD
  
  #calculate the r2 between observed SWE and estimated WEI
  func <- function(snoz)
  {
    return(data.frame(unlist(cor(snoz[ , "SNWD"], snoz[28:423]))^2))
  }
  
  WEIcor =  func(snoz)
  #WEIcor.y = ddply(snoz, .(wyear), func)
  
  #########################################################
  #select predicted SWE that was generated from highest r2
  ########################################################
  
  #identify the max correlation to select best predictor of WEI
  WEImod = data.frame("rmax" = as.numeric(max(WEIcor)))
  
  #make into matrix for calculating rowMaxs
  WEIcor.1 =  as.matrix(WEIcor[ , -1])
  
  #calculate max r2 (rmax) and make object dataframe (makes merging with alldat.1 easier)
  WEImod = data.frame("rmax" = as.numeric(rowMaxs(WEIcor.1)))
  
  #match the model name to the max r2
  WEImod$mod_ID = colnames(WEIcor.1)[max.col(WEIcor.1, ties.method="first")]

  #write number of records used to calculate r2, r2 value, and parameter suite in swemods
  swemods$ndat_swe[j] = length(SNWDval)
  swemods$mod_ID_swe[j] = WEImod$mod_ID
  
  #merge alldat.1 and WEImod
  alldat.2 = merge(alldat.1, WEImod, all.x = TRUE, all.y = TRUE)
  
  #select modeled WEI generated from best model fit
  alldat.2$pred = as.numeric(alldat.2[cbind(seq_len(nrow(alldat.2)), match(as.character(alldat.2$mod_ID), colnames(alldat.2)))])
  
  #omit na values to determine correlation between fitted and observed values
  adno = alldat.2[complete.cases(alldat.2$pred) & complete.cases(alldat.2$SNWD), ]
  
  #calculate r2 between fitted and observed vales and write to swemods
  swemods$rsq_swe[j] = cor(adno$pred, adno$SNWD)^2
  
  #omit predicted WEI values generated from data that had been flagged as NA during record completeness screening
  alldat.2$mSWE = ifelse(is.na(alldat.2$PRCPfin) == T | is.na(alldat.2$TMEAN) == T, NA, alldat.2$pred)
  
  #remove  all the columns not needed
  alldat.3 = alldat.2[ , -(28:423)]
  
  ################################################################################
  #create modeled snowfall and rainfall columns using best fit model parameter set
  ################################################################################
  
  #make pars names the same as pred
  pars = t(pars)
  pars = data.frame(pars)
  pars$pred = pred$prednam
  
  #merge aalldat.2 with pars to select model parameters used to generate mSNOW
  alldat.4 = merge(alldat.3, pars, by.x = "mod_ID", by.y = "pred", all.x = T, all.y = F)
  
  #sort data by wyear and doy
  alldat.5 = alldat.4[order(alldat.3$wyear, alldat.3$doy), ]
  
  #calculate modeled snowfall and rainfall
  alldat.5$mLIQ = ifelse(alldat.5$TMEAN > Tmax, alldat.5$PRCPfin, 0)
  alldat.5$mSNOW = ifelse(alldat.5$TMEAN < Tmin, alldat.5$PRCPfin, 0)
  
  #######################################
  #compare modeled with measured snowfall
  #######################################
  
  #omit all NAs and zero values (since zeros can occur due to failure to record data, not because snow was not present)
  snfz = alldat.5[complete.cases(alldat.5$SNOW) & complete.cases(alldat.5$mSNOW) & alldat.5$SNOW != 0, ]
  
  #determine number of records in dataset
  SNOWval = snfz$SNOW
  
  #write start year, end year, number of records used to calculate r2, and r2 value in swemods
  swemods$ndat_snow[j] = length(SNOWval)
  swemods$wyear_lo_snow[j] = min(snfz$wyear)
  swemods$wyear_hi_snow[j] = max(snfz$wyear)
  swemods$nyears_snow[j] = max(snfz$wyear) - min(snfz$wyear)
  
  swemods$rsq_snow[j] = cor(snfz$SNOW, snfz$mSNOW)^2
  
  #remove columns used in snow modeling no longer needed in downstream analysis
  alldat.6 = alldat.5[ , c("year", "month", "ym", "doy", "DATE", "station", 
                           "PRCPfin", "SNOWfin", "SNWDfin", "TMAXfin", "TMINfin", 
                           "mod_ID", "mSWE", "mSNOW", "mLIQ")]
  
  #create output table name
  out_name = paste("Contosta_Projects/Daily_modSWE/", 
                   sub(":", "_", station), 
                   "_", 
                   "Daily_modSWE", 
                   ".csv", 
                   sep="")
  
  write.csv(alldat.6, out_name)
  
  #end loop

}

write.table(swemods, "swe_mods_output.csv", col.names = T, row.names = F, sep = ",")

